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
#include "calculatedFvPatchFields.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "indirectPrimitivePatch.H"
#include "syncTools.H"

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
    psiLastSolved_(-GREAT),
    stokesInst_(0),
    stokesForceInst_(Zero),
    stokesMomentInst_(Zero),
    stokesSum_(0),
    stokesForceSum_(Zero),
    stokesMomentSum_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::middleFieldForm::reset()
{
    surfaceForce_ = Zero;    surfaceMoment_ = Zero;
    elevationForce_ = Zero;  elevationMoment_ = Zero;
    stripForce_ = Zero;      stripMoment_ = Zero;
    stokesInst_ = 0;
    stokesForceInst_ = Zero; stokesMomentInst_ = Zero;
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

                // Mass-balance diagnostics: Stokes transport through the strip
                // and the part of fStrip it carries
                const vector fStokes(-rhoRef_*z*(u & n)*W*(0.5*L));
                stokesInst_ += z*(u & n)*(0.5*L);
                stokesForceInst_ += fStokes;
                stokesMomentInst_ += (mid - CofR_) ^ fStokes;
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

    // The dynamic free-surface and gradient-form hull forcing have been
    // removed (see the class description)
    for (const word key : {"psiBarForcing", "psiBarHullForcing"})
    {
        if (dict.found(key))
        {
            WarningInFunction
                << key << " is obsolete and ignored: the forcing is always"
                << " kinematic on the free surface and in divergence form on"
                << " the hull" << endl;
        }
    }

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
    hullP_.setSize(nHull, Zero);
    hullR_.setSize(nHull, Zero);
    fsQ_.setSize(pbm[freeSurfacePatchID_].size(), Zero);
    psiTAcc_ = 0;
    psiStarted_ = false;
    psiLastSolved_ = -GREAT;
    stokesSum_ = 0;
    stokesForceSum_ = Zero;
    stokesMomentSum_ = Zero;

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

    // Hull: X_n v_t and -X.(n x omega) for the divergence form, with
    // X = xi + theta x r, v = u - Xdot and n the patch normal (out of the
    // fluid; the condition holds for either orientation used consistently,
    // and fixedGradient takes the derivative along this same normal).
    // v.n = 0 by the first-order body condition; it is projected so that
    // only the tangential part enters.
    label off = 0;
    for (const label patchi : psiBarHullIDs_)
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        const vectorField nf(fp.nf());
        const vectorField& Cf = fp.Cf();
        const vectorField& Ub = U.boundaryField()[patchi];

        forAll(nf, i)
        {
            const vector r(Cf[i] - psiBarRotRef_);
            const vector X(xi + (theta ^ r));
            const vector Xdot(xiDot + (omega ^ r));
            const vector v(Ub[i] - Xdot);
            const vector vt(v - (v & nf[i])*nf[i]);

            hullP_[off + i] += (X & nf[i])*vt*dt;
            hullR_[off + i] += -(X & (nf[i] ^ omega))*dt;
        }
        off += nf.size();
    }

    // Free surface: the Stokes transport zeta u_h, on all faces (the radius
    // cut goes on its divergence)
    {
        const scalarField& zP = zeta.boundaryField()[freeSurfacePatchID_];
        const vectorField& Uf = U.boundaryField()[freeSurfacePatchID_];

        forAll(zP, i)
        {
            fsQ_[i] += zP[i]*vector(Uf[i].x(), Uf[i].y(), 0)*dt;
        }
    }

    // Strip Stokes transport, reduced in execute()
    stokesSum_ += stokesInst_*dt;
    stokesForceSum_ += stokesForceInst_*dt;
    stokesMomentSum_ += stokesMomentInst_*dt;

    psiTAcc_ += dt;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::middleFieldForm::fsKinematicForcing(scalar& sWL) const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];
    const label fsStart = fsPatch.start();
    const labelList& fc = fsPatch.faceCells();

    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const cellList& cells = mesh_.cells();
    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();

    // Mean Stokes transport <zeta u_h> per free-surface face
    const vectorField q(fsQ_/max(psiTAcc_, SMALL));

    // ... put in the top-layer cells and interpolated to their side faces,
    // which gives one value per edge, the same seen from either side, and
    // across processor boundaries
    volVectorField Qv
    (
        IOobject
        (
            "psiBarStokesTransport",
            mesh_.time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimArea/dimTime, Zero),
        calculatedFvPatchVectorField::typeName
    );
    forAll(fc, i)
    {
        Qv.primitiveFieldRef()[fc[i]] = q[i];
    }
    Qv.correctBoundaryConditions();
    const tmp<surfaceVectorField> tQf(linearInterpolate(Qv));
    const surfaceVectorField& Qf = tQf();

    boolList isHull(pbm.size(), false);
    for (const label patchi : psiBarHullIDs_) isHull[patchi] = true;

    auto tdiv = tmp<scalarField>::New(fsPatch.size(), Zero);
    scalarField& div = tdiv.ref();
    sWL = 0;

    forAll(fc, i)
    {
        const label fsFacei = fsStart + i;
        const face& fsFace = faces[fsFacei];
        const point& c = faceCentres[fsFacei];
        scalar flux = 0;

        // The edges of the face are the segments it shares with the other
        // faces of its cell (as in waterlineIntegral)
        for (const label facej : cells[fc[i]])
        {
            if (facej == fsFacei) continue;

            point p0(Zero), p1(Zero);
            if (!sharedSegment(faces[facej], fsFace, pts, p0, p1)) continue;

            // In-plane edge normal, out of the face, times the edge length
            vector m(p1.y() - p0.y(), p0.x() - p1.x(), 0);
            if ((m & (0.5*(p0 + p1) - c)) < 0) m = -m;
            if (mag(m) < SMALL) continue;

            vector qe(q[i]);
            bool onHull = false;
            if (facej < mesh_.nInternalFaces())
            {
                qe = Qf[facej];
            }
            else
            {
                const label patchi = pbm.whichPatch(facej);
                const polyPatch& pp = pbm[patchi];
                if (pp.coupled())
                {
                    qe = Qf.boundaryField()[patchi][pp.whichFace(facej)];
                }
                else if
                (
                    isA<symmetryPlanePolyPatch>(pp)
                 || isA<symmetryPolyPatch>(pp)
                )
                {
                    continue;   // no transport through a symmetry plane
                }
                else if (isHull[patchi])
                {
                    onHull = true;
                }
            }

            const scalar fe = (qe & m);
            flux += fe;
            if (onHull) sWL += fe;
        }

        div[i] = flux/max(mag(faceAreas[fsFacei]), VSMALL);

        if (psiBarRadius_ > 0)
        {
            vector d(c - CofR_);
            d.z() = 0;
            if (mag(d) > psiBarRadius_) div[i] = 0;
        }
    }

    reduce(sWL, sumOp<scalar>());

    return tdiv;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::middleFieldForm::hullDivergenceForcing() const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& faceAreas = mesh_.faceAreas();
    const scalar norm = 1.0/max(psiTAcc_, SMALL);

    // All hull faces, in the order of hullP_ / hullR_, as one patch.  The
    // edges are taken from its own topology (not from the wall cells, which
    // at a convex corner may touch the hull along an edge only).
    labelList hullFaces(hullP_.size());
    {
        label k = 0;
        for (const label patchi : psiBarHullIDs_)
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp, i) hullFaces[k++] = pp.start() + i;
        }
    }
    const indirectPrimitivePatch hp
    (
        IndirectList<face>(mesh_.faces(), hullFaces),
        mesh_.points()
    );
    const labelList& mp = hp.meshPoints();
    const faceList& lf = hp.localFaces();
    const pointField& lp = hp.localPoints();

    // Point values, area weighted over the faces round the point and summed
    // across processors: an edge takes the mean of its two end points, so the
    // two faces either side of it (on any processor) see exactly the same
    // value, normal and length -- the fluxes cancel to round-off and the net
    // flux is that through the waterline plus the rotation term (checked:
    // equal to 6 digits on the KVLCC2 mesh).
    vectorField pPoint(hp.nPoints(), Zero);
    vectorField nPoint(hp.nPoints(), Zero);
    scalarField wPoint(hp.nPoints(), Zero);
    forAll(lf, fi)
    {
        const vector& S = faceAreas[hullFaces[fi]];
        const scalar magS = mag(S);
        for (const label pi : lf[fi])
        {
            pPoint[pi] += norm*hullP_[fi]*magS;
            nPoint[pi] += S;
            wPoint[pi] += magS;
        }
    }
    syncTools::syncPointList(mesh_, mp, pPoint, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh_, mp, nPoint, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh_, mp, wPoint, plusEqOp<scalar>(), scalar(0));
    forAll(pPoint, pi)
    {
        pPoint[pi] /= max(wPoint[pi], VSMALL);
        nPoint[pi] /= max(mag(nPoint[pi]), VSMALL);
    }

    auto tdiv = tmp<scalarField>::New(hullP_.size(), Zero);
    scalarField& div = tdiv.ref();

    forAll(lf, fi)
    {
        const face& f = lf[fi];
        const scalar magS = mag(faceAreas[hullFaces[fi]]);
        scalar flux = 0;

        forAll(f, k)
        {
            const label a = f[k];
            const label b = f.nextLabel(k);

            const vector e(lp[b] - lp[a]);
            const vector pe(0.5*(pPoint[a] + pPoint[b]));
            const vector ne(nPoint[a] + nPoint[b]);

            // Edge conormal: in the surface, normal to the edge, times the
            // edge length.  e x n points out of the face for the point order
            // of a boundary face (right-handed about its outward normal),
            // and the face on the other side runs the edge the other way, so
            // the two get exactly opposite conormals.
            vector m(e ^ ne);
            const scalar magM = mag(m);
            if (magM < VSMALL) continue;
            m *= mag(e)/magM;

            flux += (pe & m);
        }

        div[fi] = flux/magS + norm*hullR_[fi];
    }

    return tdiv;
}


void Foam::functionObjects::middleFieldForm::solvePsiBar
(
    const scalarField& fsGrad,
    const scalarField& hullGrad,
    const bool useFs,
    const bool useHull,
    vector& F,
    vector& M,
    scalar& qHull,
    scalar& qFs,
    scalar& phiC,
    vector& Fmass,
    vector& Mmass
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
            pf.gradient()[i] = useHull ? hullGrad[off + i] : 0;
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
            pf.gradient()[i] = useFs ? fsGrad[i] : 0;
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
    phiC = 0;
    Fmass = Zero;
    Mmass = Zero;

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

        // Mass flux of psibar out of the control volume, and the part of f
        // it carries
        const vector fm(-rhoRef_*(g & Sf)*W);
        phiC += (g & Sf);
        Fmass += fm;
        Mmass += (faceCentres[facei] - CofR_) ^ fm;
    }

    reduce(F, sumOp<vector>());
    reduce(M, sumOp<vector>());
    reduce(phiC, sumOp<scalar>());
    reduce(Fmass, sumOp<vector>());
    reduce(Mmass, sumOp<vector>());

    Info<< "    psiBar solve (" << (useFs ? "fs" : "") << (useFs && useHull ? "+" : "")
        << (useHull ? "hull" : "") << "): residual " << res
        << "  F " << F << "  Mz " << M.z() << endl;
}


void Foam::functionObjects::middleFieldForm::psiBarOutput()
{
    if (psiTAcc_ <= SMALL) return;

    const scalar norm = 1.0/psiTAcc_;

    // Boundary forcing; sWL: Stokes transport into the hull at the waterline
    scalar sWL = 0;
    const tmp<scalarField> tfsGrad(fsKinematicForcing(sWL));
    const tmp<scalarField> thullGrad(hullDivergenceForcing());
    const scalarField& fsGrad = tfsGrad();
    const scalarField& hullGrad = thullGrad();

    // Strip Stokes transport through the control surface and its momentum
    const scalar sC = norm*stokesSum_;
    const vector Fst(norm*stokesForceSum_);
    const vector Mst(norm*stokesMomentSum_);

    vector Ffs(Zero), Mfs(Zero), Fh(Zero), Mh(Zero), F(Zero), M(Zero);
    scalar qH = 0, qF = 0, phiC = 0;
    vector Fm(Zero), Mm(Zero);

    if (psiBarSplit_)
    {
        solvePsiBar(fsGrad, hullGrad, true, false, Ffs, Mfs, qH, qF, phiC, Fm, Mm);
        solvePsiBar(fsGrad, hullGrad, false, true, Fh, Mh, qH, qF, phiC, Fm, Mm);
    }
    // last: the field kept, and phiC, Fm, Mm, are those of the total
    solvePsiBar(fsGrad, hullGrad, true, true, F, M, qH, qF, phiC, Fm, Mm);

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
            writeHeader
            (
                *psiFilePtr_,
                "mass balance: Phi_C + S_C = 0 (control volume),"
                " Q_hull + S_WL = 0 (hull)"
            );
            writeCommented(*psiFilePtr_, "Time");
            writeTabbed(*psiFilePtr_, "periods");
            writeTabbed(*psiFilePtr_, "Fpsi_x\tFpsi_y\tFpsi_z\tMpsi_x\tMpsi_y\tMpsi_z");
            writeTabbed(*psiFilePtr_, "Fpsi_fs_x\tFpsi_fs_y\tMpsi_fs_z");
            writeTabbed(*psiFilePtr_, "Fpsi_hull_x\tFpsi_hull_y\tMpsi_hull_z");
            writeTabbed(*psiFilePtr_, "Q_hull\tQ_fs\tS_WL\tS_C\tPhi_C");
            writeTabbed(*psiFilePtr_, "Fstokes_x\tFstokes_y\tMstokes_z");
            writeTabbed(*psiFilePtr_, "Fmass_x\tFmass_y\tMmass_z");
            *psiFilePtr_ << endl;
        }
        writeCurrentTime(*psiFilePtr_);
        *psiFilePtr_
            << tab << nPer
            << tab << F.x() << tab << F.y() << tab << F.z()
            << tab << M.x() << tab << M.y() << tab << M.z()
            << tab << Ffs.x() << tab << Ffs.y() << tab << Mfs.z()
            << tab << Fh.x() << tab << Fh.y() << tab << Mh.z()
            << tab << qH << tab << qF << tab << sWL << tab << sC << tab << phiC
            << tab << Fst.x() << tab << Fst.y() << tab << Mst.z()
            << tab << Fm.x() << tab << Fm.y() << tab << Mm.z() << endl;
    }

    Log << type() << ' ' << name() << " psiBar (" << nPer << " periods):" << nl
        << "    F_psi " << F << "  (should be small: G&P eq.15)" << nl
        << "    Mz_psi " << M.z() << "  (free surface " << Mfs.z()
        << ", hull " << Mh.z() << ")" << nl
        << "    net flux Q: hull " << qH << ", free surface " << qF << nl
        << "    mass, control volume: Phi_C + S_C = " << phiC << " + " << sC
        << " = " << phiC + sC << nl
        << "    mass, hull:          Q_hull + S_WL = " << qH << " + " << sWL
        << " = " << qH + sWL << nl
        << "    momentum: Fstokes + Fmass = " << Fst << " + " << Fm
        << " = " << Fst + Fm << endl;

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

    if (psiBar_)
    {
        reduce(stokesInst_, sumOp<scalar>());
        reduce(stokesForceInst_, sumOp<vector>());
        reduce(stokesMomentInst_, sumOp<vector>());
        accumulatePsiBar();
    }

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
