/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "middleFieldForm.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "linear.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"
#include "OFstream.H"
#include "fvm.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(middleFieldForm, 0);
    addToRunTimeSelectionTable(functionObject, middleFieldForm, dictionary);
}

namespace
{

//- End points of the segment shared by two faces.  Conformal faces share
//  exactly two points; across a refinement interface more are shared, but they
//  are collinear, so the two furthest apart bound the segment.
bool sharedSegment
(
    const face& a,
    const face& b,
    const pointField& pts,
    point& p0,
    point& p1
)
{
    labelList shared(a.size());
    label nShared = 0;

    for (const label pointi : a)
    {
        if (b.found(pointi)) shared[nShared++] = pointi;
    }

    if (nShared < 2) return false;

    if (nShared == 2)
    {
        p0 = pts[shared[0]];
        p1 = pts[shared[1]];
        return true;
    }

    scalar dMaxSqr = -1;
    for (label i = 0; i < nShared; ++i)
    {
        for (label j = i + 1; j < nShared; ++j)
        {
            const scalar dSqr = magSqr(pts[shared[i]] - pts[shared[j]]);
            if (dSqr > dMaxSqr)
            {
                dMaxSqr = dSqr;
                p0 = pts[shared[i]];
                p1 = pts[shared[j]];
            }
        }
    }

    return true;
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::middleFieldForm::middleFieldForm
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_("U"),
    UDName_("UD"),
    PhiDName_("PhiD"),
    snGradNormal_(true),
    UsName_("Us"),
    zetaName_("zeta"),
    faceZoneName_(word::null),
    faceZoneID_(-1),
    cvPoint_(Zero),
    freeSurfacePatchName_(word::null),
    freeSurfacePatchID_(-1),
    rhoRef_(1000),
    gMag_(9.81),
    CofR_(Zero),
    psiBar_(false),
    psiBarSolver_("PhiS"),
    psiBarNCorr_(5),
    psiBarStartTime_(0),
    psiBarRadius_(-1),
    psiBarSplit_(true),
    psiBarRotRef_(Zero),
    omegaE_(0),
    psiTAcc_(0),
    psiTPrev_(0),
    psiStarted_(false),
    psiLastSolved_(-GREAT)
{
    read(dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::middleFieldForm::reset()
{
    surfaceForce_ = Zero;    surfaceMoment_ = Zero;
    elevationForce_ = Zero;  elevationMoment_ = Zero;
    stripForce_ = Zero;      stripMoment_ = Zero;
}


void Foam::functionObjects::middleFieldForm::surfaceIntegral
(
    const surfaceVectorField& Uf,
    const surfaceVectorField& UDf,
    const surfaceScalarField& snGradPhiD
)
{
    const faceZone& fz = mesh_.faceZones()[faceZoneID_];
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // faceAreas()/faceCentres() are sized nFaces, unlike Sf()/Cf() whose
    // internal fields stop at nInternalFaces.  Control-surface faces can land
    // on a processor patch after decomposition, so the primitiveMesh versions
    // are the ones that stay in bounds.
    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();

    const auto faceValue = [&](const auto& fld, const label facei)
    {
        if (facei < mesh_.nInternalFaces())
        {
            return fld[facei];
        }
        const label patchi = pbm.whichPatch(facei);
        return fld.boundaryField()[patchi][pbm[patchi].whichFace(facei)];
    };

    for (const label facei : fz)
    {
        if (facei >= mesh_.nInternalFaces())
        {
            const label patchi = pbm.whichPatch(facei);

            // A face on a processor patch is seen by both ranks; keep the
            // owner side only
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const auto& ppp =
                    refCast<const processorPolyPatch>(pbm[patchi]);
                if (!ppp.owner()) continue;
            }
        }

        vector u(faceValue(Uf, facei));

        // Orient outward, away from the body.  snGrad is taken along the mesh
        // normal, so flipping the area vector must flip it too.
        vector Sf(faceAreas[facei]);
        scalar sgn = 1;
        if ((Sf & (faceCentres[facei] - cvPoint_)) < 0)
        {
            Sf = -Sf;
            sgn = -1;
        }

        if (snGradNormal_)
        {
            // Swap the interpolated disturbance normal component for the
            // compact one.  The incident part of u is untouched.
            const vector nHat(Sf/mag(Sf));
            const scalar unCompact = sgn*faceValue(snGradPhiD, facei);
            const scalar unInterp = (faceValue(UDf, facei) & nHat);

            u += (unCompact - unInterp)*nHat;
        }

        const vector f(rhoRef_*(0.5*magSqr(u)*Sf - (u & Sf)*u));

        surfaceForce_ += f;
        surfaceMoment_ += (faceCentres[facei] - CofR_) ^ f;
    }
}


void Foam::functionObjects::middleFieldForm::waterlineIntegral()
{
    const faceZone& fz = mesh_.faceZones()[faceZoneID_];
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

    const label fsStart = fsPatch.start();
    const label fsEnd = fsStart + fsPatch.size();

    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const cellList& cells = mesh_.cells();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();

    const auto& zeta = mesh_.lookupObject<volScalarField>(zetaName_);
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& Us = mesh_.lookupObject<volVectorField>(UsName_);

    const scalarField& zetaP = zeta.boundaryField()[freeSurfacePatchID_];
    const vectorField& Up = U.boundaryField()[freeSurfacePatchID_];
    const vectorField& Wp = Us.boundaryField()[freeSurfacePatchID_];

    for (const label cvFacei : fz)
    {
        // Outward horizontal normal of the control surface
        vector n(faceAreas[cvFacei]);
        if ((n & (faceCentres[cvFacei] - cvPoint_)) < 0) n = -n;

        n.z() = 0;
        const scalar nMag = mag(n);
        if (nMag < SMALL) continue;   // horizontal face: no waterline here
        n /= nMag;

        const face& cvFace = faces[cvFacei];

        // Cells this rank owns either side of the face.  A face on a processor
        // patch has no local neighbour; the adjoining rank supplies that half.
        FixedList<label, 2> sideCells(-1);
        sideCells[0] = owner[cvFacei];
        if (cvFacei < mesh_.nInternalFaces())
        {
            sideCells[1] = neighbour[cvFacei];
        }

        for (const label celli : sideCells)
        {
            if (celli < 0) continue;

            for (const label facej : cells[celli])
            {
                if (facej < fsStart || facej >= fsEnd) continue;

                point p0(Zero), p1(Zero);
                if (!sharedSegment(faces[facej], cvFace, pts, p0, p1)) continue;

                const scalar L = mag(p1 - p0);
                if (L < SMALL) continue;

                const label fsi = facej - fsStart;
                const point mid(0.5*(p0 + p1));

                const scalar z = zetaP[fsi];
                const vector& u = Up[fsi];
                const vector& W = Wp[fsi];

                // 1/2: this cell is one of the two sides of the face
                const vector fElev(-0.5*rhoRef_*gMag_*sqr(z)*n*(0.5*L));
                const vector fStrip
                (
                    -rhoRef_*z*((u & n)*W + (W & n)*u)*(0.5*L)
                );


                elevationForce_ += fElev;
                stripForce_ += fStrip;

                elevationMoment_ += (mid - CofR_) ^ fElev;
                stripMoment_ += (mid - CofR_) ^ fStrip;
            }
        }
    }
}


void Foam::functionObjects::middleFieldForm::createFiles()
{
    if (!Pstream::master() || forceFilePtr_) return;

    forceFilePtr_ = createFile("force");
    momentFilePtr_ = createFile("moment");

    for (auto* os : {forceFilePtr_.get(), momentFilePtr_.get()})
    {
        writeHeader(*os, "Mean wave loads, midfield formulation");
        writeHeaderValue(*os, "CofR", CofR_);
        writeHeaderValue(*os, "rhoInf", rhoRef_);
        writeCommented(*os, "Time");
        writeTabbed(*os, "total_x\ttotal_y\ttotal_z");
        writeTabbed(*os, "surface_x\tsurface_y\tsurface_z");
        writeTabbed(*os, "elevation_x\televation_y\televation_z");
        writeTabbed(*os, "strip_x\tstrip_y\tstrip_z");
        *os << endl;
    }
}


void Foam::functionObjects::middleFieldForm::writeFiles()
{
    if (!Pstream::master()) return;

    createFiles();

    writeCurrentTime(forceFilePtr_());
    forceFilePtr_() << tab << force() << tab << surfaceForce_
                    << tab << elevationForce_ << tab << stripForce_ << endl;

    writeCurrentTime(momentFilePtr_());
    momentFilePtr_() << tab << moment() << tab << surfaceMoment_
                     << tab << elevationMoment_ << tab << stripMoment_ << endl;
}


// * * * * * * * * * * * * * * psi-bar terms  * * * * * * * * * * * * * * * //

void Foam::functionObjects::middleFieldForm::readPsiBar(const dictionary& dict)
{
    psiBar_ = dict.getOrDefault<Switch>("psiBar", Switch(false));
    if (!psiBar_) return;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    psiBarHullPatches_ =
        dict.getOrDefault<wordRes>("psiBarHullPatches", wordRes({"sphere"}));
    psiBarFarFieldPatches_ = dict.getOrDefault<wordRes>
    (
        "psiBarFarFieldPatches", wordRes({"inlet", "outlet", "front", "back"})
    );
    psiBarHullIDs_ = pbm.patchSet(psiBarHullPatches_).sortedToc();
    psiBarFarFieldIDs_ = pbm.patchSet(psiBarFarFieldPatches_).sortedToc();

    if (psiBarHullIDs_.empty() || psiBarFarFieldIDs_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "psiBar: no hull patch matches " << psiBarHullPatches_
            << " or no far-field patch matches " << psiBarFarFieldPatches_
            << exit(FatalIOError);
    }

    psiBarSolver_ = dict.getOrDefault<word>("psiBarSolver", "PhiS");
    psiBarNCorr_ = dict.getOrDefault<label>("psiBarCorrectors", 5);
    psiBarRadius_ = dict.getOrDefault<scalar>("psiBarRadius", -1);
    psiBarSplit_ = dict.getOrDefault<Switch>("psiBarSplit", Switch(true));

    // Rotation centre of the body motion: linBodyMotion takes the rotations
    // about (xG, 0, 0) of bodyMotionProperties
    psiBarRotRef_ = Zero;
    {
        IOobject io
        (
            "bodyMotionProperties",
            mesh_.time().constant(),
            mesh_.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        if (io.typeHeaderOk<IOdictionary>(true))
        {
            const IOdictionary bodyDict(io);
            psiBarRotRef_ = point(bodyDict.getOrDefault<scalar>("xG", 0), 0, 0);
        }
    }
    dict.readIfPresent("psiBarRotationCentre", psiBarRotRef_);

    // Encounter frequency (constant/waveConditions, as the solver), for the
    // default averaging window and the period count in the output
    omegaE_ = 0;
    {
        IOobject io
        (
            "waveConditions",
            mesh_.time().constant(),
            mesh_.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        if (io.typeHeaderOk<IOdictionary>(true))
        {
            const IOdictionary waveDict(io);
            const scalar k =
                constant::mathematical::twoPi/waveDict.get<scalar>("waveLength");
            const scalar depth = waveDict.get<scalar>("waterDepth");
            const scalar U0 = waveDict.get<scalar>("currentSpeed");
            const scalar head = waveDict.get<scalar>("headingAngle");
            omegaE_ =
                Foam::sqrt(gMag_*k*Foam::tanh(k*depth)) + k*U0*Foam::cos(head);
        }
    }

    // Averaging window: the last psiBarPeriods encounter periods before
    // endTime (the window meanLoads.py fits), unless psiBarStartTime is given
    const scalar nPer = dict.getOrDefault<scalar>("psiBarPeriods", 6);
    psiBarStartTime_ = 0;
    if (mag(omegaE_) > SMALL)
    {
        psiBarStartTime_ = max
        (
            mesh_.time().endTime().value()
          - nPer*constant::mathematical::twoPi/mag(omegaE_),
            scalar(0)
        );
    }
    dict.readIfPresent("psiBarStartTime", psiBarStartTime_);

    label nHull = 0;
    for (const label patchi : psiBarHullIDs_) nHull += pbm[patchi].size();
    hullForcing_.setSize(nHull, Zero);
    fsForcing_.setSize(pbm[freeSurfacePatchID_].size(), Zero);
    psiTAcc_ = 0;
    psiStarted_ = false;
    psiLastSolved_ = -GREAT;

    Info<< "    psiBar          : on" << nl
        << "      hull patches  " << psiBarHullPatches_ << nl
        << "      far field     " << psiBarFarFieldPatches_ << nl
        << "      rotation ref  " << psiBarRotRef_ << nl
        << "      omega_e       " << omegaE_ << nl
        << "      average from  t = " << psiBarStartTime_ << nl
        << "      fs radius     " << psiBarRadius_ << endl;
}


void Foam::functionObjects::middleFieldForm::accumulatePsiBar()
{
    const scalar t = mesh_.time().value();
    if (t < psiBarStartTime_ - 0.5*mesh_.time().deltaTValue()) return;

    if (!psiStarted_)
    {
        psiTPrev_ = t;
        psiStarted_ = true;
        return;                         // no interval to weight yet
    }

    const scalar dt = t - psiTPrev_;
    psiTPrev_ = t;
    if (dt <= 0) return;

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& zeta = mesh_.lookupObject<volScalarField>(zetaName_);
    const tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    // Body state published by linBodyMotion (global axes); absent for a
    // restrained body, where zero is the right value
    const auto bodyVec = [&](const word& nm)
    {
        const auto* ptr = mesh_.findObject<uniformDimensionedVectorField>(nm);
        return ptr ? ptr->value() : vector(Zero);
    };
    const vector xi(bodyVec("bodyDisp"));
    const vector theta(bodyVec("bodyRot"));
    const vector xiDot(bodyVec("bodyVel"));
    const vector omega(bodyVec("bodyOmega"));

    // Hull: d(psibar)/dn = -[(X.grad)grad(phi)].n - (grad(phi) - Xdot).(theta x n)
    // with n the patch normal (out of the fluid); the condition holds for
    // either orientation used consistently, and fixedGradient takes the
    // derivative along this same normal.
    label off = 0;
    for (const label patchi : psiBarHullIDs_)
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        const vectorField nf(fp.nf());
        const vectorField& Cf = fp.Cf();
        const vectorField& Ub = U.boundaryField()[patchi];
        const tensorField& gUb = gradU.boundaryField()[patchi];

        forAll(nf, i)
        {
            const vector r(Cf[i] - psiBarRotRef_);
            const vector X(xi + (theta ^ r));
            const vector Xdot(xiDot + (omega ^ r));

            // (X & gradU) = (X.grad)u, since gradU_ji = d_j u_i
            const scalar h =
              - ((X & gUb[i]) & nf[i])
              - ((Ub[i] - Xdot) & (theta ^ nf[i]));

            hullForcing_[off + i] += h*dt;
        }
        off += nf.size();
    }

    // Free surface: d(psibar)/dz = -zeta d2phi/dz2, d2phi/dz2 = -(du/dx + dv/dy)
    {
        const scalarField& zP = zeta.boundaryField()[freeSurfacePatchID_];
        const tensorField& gUf = gradU.boundaryField()[freeSurfacePatchID_];
        const vectorField& Cf = mesh_.boundary()[freeSurfacePatchID_].Cf();

        forAll(zP, i)
        {
            if (psiBarRadius_ > 0)
            {
                vector d(Cf[i] - CofR_);
                d.z() = 0;
                if (mag(d) > psiBarRadius_) continue;
            }

            const scalar phizz = -(gUf[i].xx() + gUf[i].yy());
            fsForcing_[i] += -zP[i]*phizz*dt;
        }
    }

    psiTAcc_ += dt;
}


void Foam::functionObjects::middleFieldForm::solvePsiBar
(
    const bool useFs,
    const bool useHull,
    vector& F,
    vector& M,
    scalar& qHull,
    scalar& qFs
)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!psiBarPtr_)
    {
        // Constraint patches (processor, empty, ...) keep their own type
        wordList bcTypes(pbm.size(), zeroGradientFvPatchScalarField::typeName);
        for (const label patchi : psiBarHullIDs_)
        {
            bcTypes[patchi] = fixedGradientFvPatchScalarField::typeName;
        }
        bcTypes[freeSurfacePatchID_] = fixedGradientFvPatchScalarField::typeName;
        for (const label patchi : psiBarFarFieldIDs_)
        {
            bcTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }

        psiBarPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "psiBar",
                    mesh_.time().timeName(),
                    mesh_.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedScalar(dimArea/dimTime, Zero),
                bcTypes
            )
        );
    }

    volScalarField& psi = psiBarPtr_();
    const scalar norm = 1.0/max(psiTAcc_, SMALL);

    qHull = 0;
    qFs = 0;

    label off = 0;
    for (const label patchi : psiBarHullIDs_)
    {
        auto& pf = refCast<fixedGradientFvPatchScalarField>
        (
            psi.boundaryFieldRef()[patchi]
        );
        const scalarField magSf(mag(mesh_.boundary()[patchi].Sf()));
        forAll(pf, i)
        {
            pf.gradient()[i] = useHull ? norm*hullForcing_[off + i] : 0;
            qHull += pf.gradient()[i]*magSf[i];
        }
        off += pf.size();
    }
    {
        auto& pf = refCast<fixedGradientFvPatchScalarField>
        (
            psi.boundaryFieldRef()[freeSurfacePatchID_]
        );
        const scalarField magSf(mag(mesh_.boundary()[freeSurfacePatchID_].Sf()));
        forAll(pf, i)
        {
            pf.gradient()[i] = useFs ? norm*fsForcing_[i] : 0;
            qFs += pf.gradient()[i]*magSf[i];
        }
    }
    for (const label patchi : psiBarFarFieldIDs_)
    {
        psi.boundaryFieldRef()[patchi] == 0.0;
    }
    reduce(qHull, sumOp<scalar>());
    reduce(qFs, sumOp<scalar>());

    psi.primitiveFieldRef() = 0;
    psi.correctBoundaryConditions();

    scalar res = 0;
    for (label corr = 0; corr < psiBarNCorr_; ++corr)
    {
        fvScalarMatrix psiEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), psi)
         == dimensionedScalar(psi.dimensions()/dimArea, Zero)
        );
        res = psiEqn.solve(mesh_.solverDict(psiBarSolver_)).initialResidual();
        if (res < 1e-7) break;
    }
    psi.correctBoundaryConditions();

    // Cross momentum flux through the control surface, as surfaceIntegral():
    // outward normals, owner side of processor faces, compact normal
    // derivative.  f = rho [ (W.g) Sf - (W.Sf) g - (g.Sf) W ], g = grad(psibar)
    const auto& Us = mesh_.lookupObject<volVectorField>(UsName_);
    const tmp<surfaceVectorField> tWf(fvc::interpolate(Us));
    const tmp<surfaceVectorField> tgf(fvc::interpolate(fvc::grad(psi)));
    const tmp<surfaceScalarField> tsn(fvc::snGrad(psi));

    const faceZone& fz = mesh_.faceZones()[faceZoneID_];
    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();

    const auto faceValue = [&](const auto& fld, const label facei)
    {
        if (facei < mesh_.nInternalFaces())
        {
            return fld[facei];
        }
        const label patchi = pbm.whichPatch(facei);
        return fld.boundaryField()[patchi][pbm[patchi].whichFace(facei)];
    };

    F = Zero;
    M = Zero;

    for (const label facei : fz)
    {
        if (facei >= mesh_.nInternalFaces())
        {
            const label patchi = pbm.whichPatch(facei);
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const auto& ppp = refCast<const processorPolyPatch>(pbm[patchi]);
                if (!ppp.owner()) continue;
            }
        }

        vector Sf(faceAreas[facei]);
        scalar sgn = 1;
        if ((Sf & (faceCentres[facei] - cvPoint_)) < 0)
        {
            Sf = -Sf;
            sgn = -1;
        }
        const vector nHat(Sf/mag(Sf));

        vector g(faceValue(tgf(), facei));
        g += (sgn*faceValue(tsn(), facei) - (g & nHat))*nHat;
        const vector W(faceValue(tWf(), facei));

        const vector f(rhoRef_*((W & g)*Sf - (W & Sf)*g - (g & Sf)*W));

        F += f;
        M += (faceCentres[facei] - CofR_) ^ f;
    }

    reduce(F, sumOp<vector>());
    reduce(M, sumOp<vector>());

    Info<< "    psiBar solve (" << (useFs ? "fs" : "") << (useFs && useHull ? "+" : "")
        << (useHull ? "hull" : "") << "): residual " << res
        << "  F " << F << "  Mz " << M.z() << endl;
}


void Foam::functionObjects::middleFieldForm::psiBarOutput()
{
    if (psiTAcc_ <= SMALL) return;

    vector Ffs(Zero), Mfs(Zero), Fh(Zero), Mh(Zero), F(Zero), M(Zero);
    scalar qH = 0, qF = 0;

    if (psiBarSplit_)
    {
        solvePsiBar(true, false, Ffs, Mfs, qH, qF);
        solvePsiBar(false, true, Fh, Mh, qH, qF);
    }
    solvePsiBar(true, true, F, M, qH, qF);      // last: the field kept is the total

    psiLastSolved_ = mesh_.time().value();

    const scalar nPer =
        mag(omegaE_) > SMALL
      ? psiTAcc_*mag(omegaE_)/constant::mathematical::twoPi
      : 0;

    if (mesh_.time().writeTime())
    {
        psiBarPtr_->instance() = mesh_.time().timeName();
        psiBarPtr_->write();
    }

    if (Pstream::master())
    {
        if (!psiFilePtr_)
        {
            psiFilePtr_ = createFile("psiBar");
            writeHeader(*psiFilePtr_, "psi-bar terms of the mean loads (add to force/moment)");
            writeHeaderValue(*psiFilePtr_, "CofR", CofR_);
            writeHeaderValue(*psiFilePtr_, "averaged from", psiBarStartTime_);
            writeCommented(*psiFilePtr_, "Time");
            writeTabbed(*psiFilePtr_, "periods");
            writeTabbed(*psiFilePtr_, "Fpsi_x\tFpsi_y\tFpsi_z\tMpsi_x\tMpsi_y\tMpsi_z");
            writeTabbed(*psiFilePtr_, "Fpsi_fs_x\tFpsi_fs_y\tMpsi_fs_z");
            writeTabbed(*psiFilePtr_, "Fpsi_hull_x\tFpsi_hull_y\tMpsi_hull_z");
            writeTabbed(*psiFilePtr_, "Q_hull\tQ_fs");
            *psiFilePtr_ << endl;
        }
        writeCurrentTime(*psiFilePtr_);
        *psiFilePtr_
            << tab << nPer
            << tab << F.x() << tab << F.y() << tab << F.z()
            << tab << M.x() << tab << M.y() << tab << M.z()
            << tab << Ffs.x() << tab << Ffs.y() << tab << Mfs.z()
            << tab << Fh.x() << tab << Fh.y() << tab << Mh.z()
            << tab << qH << tab << qF << endl;
    }

    Log << type() << ' ' << name() << " psiBar (" << nPer << " periods):" << nl
        << "    F_psi " << F << "  (should be small: G&P eq.15)" << nl
        << "    Mz_psi " << M.z() << "  (free surface " << Mfs.z()
        << ", hull " << Mh.z() << ")" << nl
        << "    net flux Q: hull " << qH << ", free surface " << qF << endl;

    setResult("psiBarForce", F);
    setResult("psiBarMoment", M);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::middleFieldForm::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    dict.readIfPresent("U", UName_);
    dict.readIfPresent("UD", UDName_);
    dict.readIfPresent("PhiD", PhiDName_);
    snGradNormal_ = dict.getOrDefault<Switch>("snGradNormal", Switch(true));
    dict.readIfPresent("Us", UsName_);
    dict.readIfPresent("zeta", zetaName_);

    dict.readEntry("faceZone", faceZoneName_);
    faceZoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
    if (faceZoneID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "faceZone " << faceZoneName_ << " not found" << exit(FatalIOError);
    }

    dict.readEntry("freeSurfacePatch", freeSurfacePatchName_);
    freeSurfacePatchID_ =
        mesh_.boundaryMesh().findPatchID(freeSurfacePatchName_);
    if (freeSurfacePatchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "freeSurfacePatch " << freeSurfacePatchName_ << " not found"
            << exit(FatalIOError);
    }

    dict.readEntry("cvPoint", cvPoint_);
    dict.readEntry("CofR", CofR_);
    dict.readEntry("rhoInf", rhoRef_);
    gMag_ = dict.getOrDefault<scalar>("gMag", 9.81);

    Info<< type() << ' ' << name() << ':' << nl
        << "    control surface : " << faceZoneName_ << nl
        << "    free surface    : " << freeSurfacePatchName_ << nl
        << "    rhoInf          : " << rhoRef_ << endl;

    readPsiBar(dict);

    return true;
}


bool Foam::functionObjects::middleFieldForm::execute()
{
    reset();

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& UD = mesh_.lookupObject<volVectorField>(UDName_);
    const auto& PhiD = mesh_.lookupObject<volScalarField>(PhiDName_);

    // fvc::interpolate honours interpolationSchemes, so the face value of the
    // velocity can be given a skewness correction from fvSchemes; the old
    // linearInterpolate() hard-coded uncorrected geometric weights.
    const tmp<surfaceVectorField> tUf(fvc::interpolate(U));
    const tmp<surfaceVectorField> tUDf(fvc::interpolate(UD));
    const tmp<surfaceScalarField> tSnGradPhiD(fvc::snGrad(PhiD));

    surfaceIntegral(tUf(), tUDf(), tSnGradPhiD());
    waterlineIntegral();

    reduce(surfaceForce_, sumOp<vector>());
    reduce(surfaceMoment_, sumOp<vector>());
    reduce(elevationForce_, sumOp<vector>());
    reduce(elevationMoment_, sumOp<vector>());
    reduce(stripForce_, sumOp<vector>());
    reduce(stripMoment_, sumOp<vector>());


    Log << type() << ' ' << name() << " write:" << nl
        << "    total     " << force() << nl
        << "    surface   " << surfaceForce_ << nl
        << "    elevation " << elevationForce_ << nl
        << "    strip     " << stripForce_ << endl;

    setResult("force", force());
    setResult("moment", moment());

    if (psiBar_) accumulatePsiBar();

    return true;
}


bool Foam::functionObjects::middleFieldForm::write()
{
    if (writeToFile())
    {
        writeFiles();
    }

    if (psiBar_ && mesh_.time().writeTime())
    {
        psiBarOutput();
    }

    return true;
}


bool Foam::functionObjects::middleFieldForm::end()
{
    if (psiBar_ && psiLastSolved_ != mesh_.time().value())
    {
        psiBarOutput();
    }

    return true;
}


Foam::vector Foam::functionObjects::middleFieldForm::force() const
{
    return surfaceForce_ + elevationForce_ + stripForce_;
}


Foam::vector Foam::functionObjects::middleFieldForm::moment() const
{
    return surfaceMoment_ + elevationMoment_ + stripMoment_;
}


// ************************************************************************* //
