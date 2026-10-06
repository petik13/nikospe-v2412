/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "middleFieldFormRot.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "linear.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "processorPolyPatch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(middleFieldFormRot, 0);
    addToRunTimeSelectionTable(functionObject, middleFieldFormRot, dictionary);
}

namespace
{

//- End points of the segment shared by two faces (as middleFieldForm)
bool sharedSegmentRot
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

Foam::functionObjects::middleFieldFormRot::middleFieldFormRot
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_("U"),
    UsName_("Us"),
    VSName_("VS"),
    zetaName_("zeta"),
    UDName_("UD"),
    PhiDName_("PhiD"),
    UIName_("UI"),
    zetaIName_("zetaI"),
    snGradNormal_(true),
    cornerTerm_(true),
    cornerReported_(false),
    faceZoneName_(word::null),
    faceZoneID_(-1),
    cvPoint_(Zero),
    freeSurfacePatchName_(word::null),
    freeSurfacePatchID_(-1),
    hullPatchName_("sphere"),
    hullPatchID_(-1),
    rhoRef_(1000),
    gMag_(9.81),
    CofR_(Zero),
    OmegaDict_(Zero),
    cvBox_(),
    fsInside_(),
    fsInsideBuilt_(false),
    P_(Zero), PI_(Zero), Hz_(0), Q_(0),
    Pc_(Zero), Hc_(0), Qc_(0),
    Pold_(Zero), HzOld_(0), tOld_(0), lastTimeIndex_(-1), haveOld_(false)
{
    reset();
    storageForce_ = Zero;
    storageMoment_ = Zero;
    read(dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::middleFieldFormRot::reset()
{
    surfaceForce_ = Zero;    surfaceMoment_ = Zero;
    elevationForce_ = Zero;  elevationMoment_ = Zero;
    stripForce_ = Zero;      stripMoment_ = Zero;
    coriolisForce_ = Zero;   coriolisMoment_ = Zero;
    centripetalForce_ = Zero; centripetalMoment_ = Zero;
}


Foam::vector Foam::functionObjects::middleFieldFormRot::Omega() const
{
    const auto* ptr = mesh_.findObject<uniformDimensionedVectorField>("frameOmega");
    return ptr ? ptr->value() : OmegaDict_;
}


void Foam::functionObjects::middleFieldFormRot::buildInside()
{
    if (fsInsideBuilt_) return;
    fsInsideBuilt_ = true;

    // Horizontal extent of the control volume from its faceZone.  The zone
    // faces lie ON the box boundary, so the free-surface faces of the cells
    // inside have their centres strictly inside it.
    const faceZone& fz = mesh_.faceZones()[faceZoneID_];
    const vectorField& fc = mesh_.faceCentres();

    point bbMin(GREAT, GREAT, GREAT);
    point bbMax(-GREAT, -GREAT, -GREAT);
    for (const label facei : fz)
    {
        bbMin = min(bbMin, fc[facei]);
        bbMax = max(bbMax, fc[facei]);
    }
    reduce(bbMin, minOp<point>());
    reduce(bbMax, maxOp<point>());
    cvBox_ = boundBox(bbMin, bbMax);

    const fvPatch& fsp = mesh_.boundary()[freeSurfacePatchID_];
    const vectorField& Cf = fsp.Cf();
    const scalarField& magSf = fsp.magSf();

    DynamicList<label> inside(fsp.size());
    scalar area = 0;

    forAll(Cf, i)
    {
        const point& c = Cf[i];
        if
        (
            c.x() > bbMin.x() && c.x() < bbMax.x()
         && c.y() > bbMin.y() && c.y() < bbMax.y()
        )
        {
            inside.append(i);
            area += magSf[i];
        }
    }
    fsInside_.transfer(inside);

    reduce(area, sumOp<scalar>());
    const label nIn = returnReduce(fsInside_.size(), sumOp<label>());

    Info<< type() << ' ' << name() << ": control volume x "
        << bbMin.x() << " .. " << bbMax.x() << ", y "
        << bbMin.y() << " .. " << bbMax.y() << nl
        << "    free surface inside: " << nIn << " faces, area " << area
        << " m^2 (box " << (bbMax.x() - bbMin.x())*(bbMax.y() - bbMin.y())
        << " m^2, the difference is the waterplane)" << endl;
}


void Foam::functionObjects::middleFieldFormRot::surfaceIntegral
(
    const surfaceVectorField& Uf,
    const surfaceVectorField& UDf,
    const surfaceScalarField& snGradPhiD
)
{
    // As middleFieldForm::surfaceIntegral
    const faceZone& fz = mesh_.faceZones()[faceZoneID_];
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

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

            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const auto& ppp =
                    refCast<const processorPolyPatch>(pbm[patchi]);
                if (!ppp.owner()) continue;
            }
        }

        vector u(faceValue(Uf, facei));

        vector Sf(faceAreas[facei]);
        scalar sgn = 1;
        if ((Sf & (faceCentres[facei] - cvPoint_)) < 0)
        {
            Sf = -Sf;
            sgn = -1;
        }

        if (snGradNormal_)
        {
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


void Foam::functionObjects::middleFieldFormRot::waterlineIntegral()
{
    // As middleFieldForm::waterlineIntegral; W from Us is the relative flow,
    // which in the rotating frame contains -Omega x r
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
        vector n(faceAreas[cvFacei]);
        if ((n & (faceCentres[cvFacei] - cvPoint_)) < 0) n = -n;

        n.z() = 0;
        const scalar nMag = mag(n);
        if (nMag < SMALL) continue;
        n /= nMag;

        const face& cvFace = faces[cvFacei];

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
                if (!sharedSegmentRot(faces[facej], cvFace, pts, p0, p1)) continue;

                const scalar L = mag(p1 - p0);
                if (L < SMALL) continue;

                const label fsi = facej - fsStart;
                const point mid(0.5*(p0 + p1));

                const scalar z = zetaP[fsi];
                const vector& u = Up[fsi];
                const vector& W = Wp[fsi];

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


void Foam::functionObjects::middleFieldFormRot::momentumIntegrals()
{
    buildInside();

    P_ = Zero;
    PI_ = Zero;
    Hz_ = 0;
    Q_ = 0;

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& zeta = mesh_.lookupObject<volScalarField>(zetaName_);
    const auto* UIptr = mesh_.findObject<volVectorField>(UIName_);
    const auto* zetaIptr = mesh_.findObject<volScalarField>(zetaIName_);

    // --- Free-surface strip: rho int_AF zeta grad_h(phi) dA ----------------
    {
        const label fsi = freeSurfacePatchID_;
        const fvPatch& fsp = mesh_.boundary()[fsi];
        const vectorField& Cf = fsp.Cf();
        const scalarField& magSf = fsp.magSf();
        const scalarField& zp = zeta.boundaryField()[fsi];
        const vectorField& up = U.boundaryField()[fsi];

        for (const label i : fsInside_)
        {
            const vector uh(up[i].x(), up[i].y(), 0);
            const vector r(Cf[i] - CofR_);
            const scalar zA = rhoRef_*zp[i]*magSf[i];

            P_  += zA*uh;
            Hz_ += zA*(r.x()*uh.y() - r.y()*uh.x());
            Q_  += zA*(r.x()*uh.x() + r.y()*uh.y());
        }

        // Incident self-part, for monitoring the cancellation with the strip
        if (UIptr && zetaIptr)
        {
            const scalarField& zIp = zetaIptr->boundaryField()[fsi];
            const vectorField& uIp = UIptr->boundaryField()[fsi];

            for (const label i : fsInside_)
            {
                PI_ += rhoRef_*zIp[i]*magSf[i]*vector(uIp[i].x(), uIp[i].y(), 0);
            }
        }
    }

    // --- Hull: fluid gained where the hull moves away, rho (xi.n) grad(phi) --
    // n is the patch normal, out of the fluid (into the body)
    if (hullPatchID_ >= 0)
    {
        const auto* dispPtr =
            mesh_.findObject<uniformDimensionedVectorField>("bodyDisp");
        const auto* rotPtr =
            mesh_.findObject<uniformDimensionedVectorField>("bodyRot");

        if (dispPtr && rotPtr)
        {
            const vector xi(dispPtr->value());
            const vector th(rotPtr->value());

            const fvPatch& hp = mesh_.boundary()[hullPatchID_];
            const vectorField& Cf = hp.Cf();
            const vectorField nf(hp.nf());
            const scalarField& magSf = hp.magSf();
            const vectorField& uh = U.boundaryField()[hullPatchID_];

            forAll(Cf, i)
            {
                const vector r(Cf[i] - CofR_);
                const vector S(xi + (th ^ r));
                const scalar sA = rhoRef_*(S & nf[i])*magSf[i];
                const vector& u = uh[i];

                P_  += sA*vector(u.x(), u.y(), 0);
                Hz_ += sA*(r.x()*u.y() - r.y()*u.x());
                Q_  += sA*(r.x()*u.x() + r.y()*u.y());
            }
        }
    }

    // --- Waterline corner: adds its local parts to P_, Hz_, Q_ --------------
    cornerIntegrals(Omega());

    reduce(P_, sumOp<vector>());
    reduce(PI_, sumOp<vector>());
    reduce(Hz_, sumOp<scalar>());
    reduce(Q_, sumOp<scalar>());
}


void Foam::functionObjects::middleFieldFormRot::cornerIntegrals
(
    const vector& Om
)
{
    // The fluid between z = 0 and the free surface that the hull leaves or
    // takes.  Per waterline edge of length L its mass is
    //     m = rho zeta (xi.n_h) L
    // (wall-sided at the waterline), moving with the relative steady flow W.
    // Its parts of P, H, Q are added to P_, Hz_, Q_ (local sums, reduced by
    // the caller); the centripetal term is its part of
    // -rho Omega x int_D V_S dV.
    Pc_ = Zero;
    Hc_ = 0;
    Qc_ = 0;
    centripetalForce_ = Zero;
    centripetalMoment_ = Zero;

    if (!cornerTerm_ || hullPatchID_ < 0) return;

    const auto* dispPtr =
        mesh_.findObject<uniformDimensionedVectorField>("bodyDisp");
    const auto* rotPtr =
        mesh_.findObject<uniformDimensionedVectorField>("bodyRot");
    if (!dispPtr || !rotPtr) return;

    const vector xi(dispPtr->value());
    const vector th(rotPtr->value());

    const auto& zeta = mesh_.lookupObject<volScalarField>(zetaName_);
    const auto& Us = mesh_.lookupObject<volVectorField>(UsName_);
    const auto* VSptr = mesh_.findObject<volVectorField>(VSName_);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& hullPatch = pbm[hullPatchID_];
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];
    const label fsStart = fsPatch.start();
    const label fsEnd = fsStart + fsPatch.size();

    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const cellList& cells = mesh_.cells();
    const labelList& hullCells = hullPatch.faceCells();
    const vectorField nf(mesh_.boundary()[hullPatchID_].nf());

    const scalarField& zetaFs = zeta.boundaryField()[freeSurfacePatchID_];
    const vectorField& WFs = Us.boundaryField()[freeSurfacePatchID_];
    const vectorField& WHull = Us.boundaryField()[hullPatchID_];

    vector mVS(Zero);       // sum m V_S,h
    scalar mrVS = 0;        // sum m r_h.V_S
    scalar lengthSum = 0;
    label nEdges = 0;

    forAll(hullCells, i)
    {
        // Horizontal unit normal, out of the fluid (into the body)
        vector nh(nf[i].x(), nf[i].y(), 0);
        const scalar nhMag = mag(nh);
        if (nhMag < SMALL) continue;
        nh /= nhMag;

        const face& hullFace = faces[hullPatch.start() + i];

        for (const label facej : cells[hullCells[i]])
        {
            if (facej < fsStart || facej >= fsEnd) continue;

            point p0(Zero), p1(Zero);
            if (!sharedSegmentRot(faces[facej], hullFace, pts, p0, p1)) continue;

            const scalar L = mag(p1 - p0);
            if (L < SMALL) continue;

            const label fsi = facej - fsStart;
            const vector r(0.5*(p0 + p1) - CofR_);
            const vector S(xi + (th ^ r));

            const scalar m = rhoRef_*zetaFs[fsi]*(S & nh)*L;

            vector W(0.5*(WHull[i] + WFs[fsi]));
            W.z() = 0;

            Pc_ += m*W;
            Hc_ += m*(r.x()*W.y() - r.y()*W.x());
            Qc_ += m*(r.x()*W.x() + r.y()*W.y());

            if (VSptr)
            {
                vector V
                (
                    0.5*
                    (
                        VSptr->boundaryField()[hullPatchID_][i]
                      + VSptr->boundaryField()[freeSurfacePatchID_][fsi]
                    )
                );
                V.z() = 0;

                mVS += m*V;
                mrVS += m*(r.x()*V.x() + r.y()*V.y());
            }

            lengthSum += L;
            ++nEdges;
        }
    }

    P_ += Pc_;
    Hz_ += Hc_;
    Q_ += Qc_;

    reduce(Pc_, sumOp<vector>());
    reduce(Hc_, sumOp<scalar>());
    reduce(Qc_, sumOp<scalar>());
    reduce(mVS, sumOp<vector>());
    reduce(mrVS, sumOp<scalar>());

    centripetalForce_ = -(Om ^ mVS);
    centripetalMoment_ = vector(0, 0, -Om.z()*mrVS);

    if (!cornerReported_)
    {
        cornerReported_ = true;
        reduce(lengthSum, sumOp<scalar>());
        reduce(nEdges, sumOp<label>());

        Info<< type() << ' ' << name() << ": waterline corner, "
            << nEdges << " edges, length " << lengthSum << " m";
        if (!VSptr)
        {
            Info<< " (no field " << VSName_
                << ": the centripetal term is left out)";
        }
        Info<< endl;
    }
}


void Foam::functionObjects::middleFieldFormRot::createFiles()
{
    if (!Pstream::master() || forceFilePtr_) return;

    forceFilePtr_ = createFile("force");
    momentFilePtr_ = createFile("moment");
    momentumFilePtr_ = createFile("momentum");

    for (auto* os : {forceFilePtr_.get(), momentFilePtr_.get()})
    {
        writeHeader(*os, "Mean wave loads, midfield formulation, rotating frame");
        writeHeaderValue(*os, "CofR", CofR_);
        writeHeaderValue(*os, "rhoInf", rhoRef_);
        writeHeader
        (
            *os,
            "total = chen + coriolis + storage + centripetal; chen = surface"
            " + elevation + strip (middleFieldForm); the waterline corner is in"
            " coriolis, storage (through P) and centripetal"
        );
        writeHeaderValue(*os, "cornerTerm", word(cornerTerm_ ? "on" : "off"));
        writeCommented(*os, "Time");
        writeTabbed(*os, "total_x\ttotal_y\ttotal_z");
        writeTabbed(*os, "chen_x\tchen_y\tchen_z");
        writeTabbed(*os, "surface_x\tsurface_y\tsurface_z");
        writeTabbed(*os, "elevation_x\televation_y\televation_z");
        writeTabbed(*os, "strip_x\tstrip_y\tstrip_z");
        writeTabbed(*os, "coriolis_x\tcoriolis_y\tcoriolis_z");
        writeTabbed(*os, "storage_x\tstorage_y\tstorage_z");
        writeTabbed(*os, "centripetal_x\tcentripetal_y\tcentripetal_z");
        *os << endl;
    }

    writeHeader(momentumFilePtr_(), "Relative momentum in the control volume");
    writeHeader
    (
        momentumFilePtr_(),
        "P = rho int_AF zeta u_h dA + rho int_SB (xi.n) u_h dS"
        " + int_GH m_c W_h dl;"
        " P_I = rho int_AF zetaI uI_h dA (incident self-part);"
        " Pc, Hc, Qc = waterline-corner parts (included in P, Hz, Q)"
    );
    writeCommented(momentumFilePtr_(), "Time");
    writeTabbed
    (
        momentumFilePtr_(),
        "P_x\tP_y\tP_z\tPI_x\tPI_y\tPI_z\tHz\tQ\tOmega_z"
        "\tPc_x\tPc_y\tPc_z\tHc\tQc"
    );
    momentumFilePtr_() << endl;
}


void Foam::functionObjects::middleFieldFormRot::writeFiles()
{
    if (!Pstream::master()) return;

    createFiles();

    writeCurrentTime(forceFilePtr_());
    forceFilePtr_() << tab << force() << tab << chenForce()
                    << tab << surfaceForce_ << tab << elevationForce_
                    << tab << stripForce_ << tab << coriolisForce_
                    << tab << storageForce_ << tab << centripetalForce_
                    << endl;

    writeCurrentTime(momentFilePtr_());
    momentFilePtr_() << tab << moment() << tab << chenMoment()
                     << tab << surfaceMoment_ << tab << elevationMoment_
                     << tab << stripMoment_ << tab << coriolisMoment_
                     << tab << storageMoment_ << tab << centripetalMoment_
                     << endl;

    writeCurrentTime(momentumFilePtr_());
    momentumFilePtr_() << tab << P_.x() << tab << P_.y() << tab << P_.z()
                       << tab << PI_.x() << tab << PI_.y() << tab << PI_.z()
                       << tab << Hz_ << tab << Q_
                       << tab << Omega().z()
                       << tab << Pc_.x() << tab << Pc_.y() << tab << Pc_.z()
                       << tab << Hc_ << tab << Qc_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::middleFieldFormRot::read(const dictionary& dict)
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
    dict.readIfPresent("VS", VSName_);
    cornerTerm_ = dict.getOrDefault<Switch>("cornerTerm", Switch(true));
    cornerReported_ = false;
    dict.readIfPresent("zeta", zetaName_);
    dict.readIfPresent("UI", UIName_);
    dict.readIfPresent("zetaI", zetaIName_);

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

    hullPatchName_ = dict.getOrDefault<word>("hullPatch", "sphere");
    hullPatchID_ = mesh_.boundaryMesh().findPatchID(hullPatchName_);
    if (hullPatchID_ < 0)
    {
        WarningInFunction
            << "hullPatch " << hullPatchName_ << " not found: the hull part"
               " of the control-volume momentum is left out" << endl;
    }

    dict.readEntry("cvPoint", cvPoint_);
    dict.readEntry("CofR", CofR_);
    dict.readEntry("rhoInf", rhoRef_);
    gMag_ = dict.getOrDefault<scalar>("gMag", 9.81);
    OmegaDict_ = dict.getOrDefault<vector>("Omega", Zero);

    fsInsideBuilt_ = false;

    Info<< type() << ' ' << name() << ':' << nl
        << "    control surface : " << faceZoneName_ << nl
        << "    free surface    : " << freeSurfacePatchName_ << nl
        << "    hull            : " << hullPatchName_ << nl
        << "    rhoInf          : " << rhoRef_ << nl
        << "    waterline corner: " << (cornerTerm_ ? "on" : "off") << nl
        << "    Omega           : frameOmega if registered, else "
        << OmegaDict_ << nl
        << "    CofR            : " << CofR_
        << "  (must be the rotation centre of the frame)" << endl;

    return true;
}


bool Foam::functionObjects::middleFieldFormRot::execute()
{
    reset();

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& UD = mesh_.lookupObject<volVectorField>(UDName_);
    const auto& PhiD = mesh_.lookupObject<volScalarField>(PhiDName_);

    const tmp<surfaceVectorField> tUf(fvc::interpolate(U));
    const tmp<surfaceVectorField> tUDf(fvc::interpolate(UD));
    const tmp<surfaceScalarField> tSnGradPhiD(fvc::snGrad(PhiD));

    // --- Chen (middleFieldForm) terms ---------------------------------------
    surfaceIntegral(tUf(), tUDf(), tSnGradPhiD());
    waterlineIntegral();

    reduce(surfaceForce_, sumOp<vector>());
    reduce(surfaceMoment_, sumOp<vector>());
    reduce(elevationForce_, sumOp<vector>());
    reduce(elevationMoment_, sumOp<vector>());
    reduce(stripForce_, sumOp<vector>());
    reduce(stripMoment_, sumOp<vector>());

    // --- Rotating-frame terms -----------------------------------------------
    momentumIntegrals();

    const vector Om(Omega());

    // Coriolis: -2 Omega x P, and -2 Omega_z Q for yaw
    coriolisForce_ = -2.0*(Om ^ P_);
    coriolisMoment_ = vector(0, 0, -2.0*Om.z()*Q_);

    // Storage: backward difference, once per time step (telescopes over a
    // period, so the period mean is exact)
    const label timeIndex = mesh_.time().timeIndex();
    const scalar t = mesh_.time().value();

    if (timeIndex != lastTimeIndex_)
    {
        if (haveOld_ && t > tOld_ + VSMALL)
        {
            storageForce_ = -(P_ - Pold_)/(t - tOld_);
            storageMoment_ = vector(0, 0, -(Hz_ - HzOld_)/(t - tOld_));
        }
        else
        {
            storageForce_ = Zero;
            storageMoment_ = Zero;
        }

        Pold_ = P_;
        HzOld_ = Hz_;
        tOld_ = t;
        lastTimeIndex_ = timeIndex;
        haveOld_ = true;
    }

    Log << type() << ' ' << name() << " write:" << nl
        << "    total     " << force() << nl
        << "    chen      " << chenForce() << nl
        << "    coriolis  " << coriolisForce_ << nl
        << "    storage   " << storageForce_ << nl
        << "    centripetal " << centripetalForce_ << nl
        << "    P         " << P_ << "   P_I " << PI_
        << "   P_corner " << Pc_ << endl;

    setResult("force", force());
    setResult("moment", moment());
    setResult("chenForce", chenForce());
    setResult("chenMoment", chenMoment());

    return true;
}


bool Foam::functionObjects::middleFieldFormRot::write()
{
    if (writeToFile())
    {
        writeFiles();
    }

    return true;
}


Foam::vector Foam::functionObjects::middleFieldFormRot::chenForce() const
{
    return surfaceForce_ + elevationForce_ + stripForce_;
}


Foam::vector Foam::functionObjects::middleFieldFormRot::chenMoment() const
{
    return surfaceMoment_ + elevationMoment_ + stripMoment_;
}


Foam::vector Foam::functionObjects::middleFieldFormRot::force() const
{
    return chenForce() + coriolisForce_ + storageForce_ + centripetalForce_;
}


Foam::vector Foam::functionObjects::middleFieldFormRot::moment() const
{
    return chenMoment() + coriolisMoment_ + storageMoment_ + centripetalMoment_;
}


// ************************************************************************* //
