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
    CofR_(Zero)
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

    return true;
}


bool Foam::functionObjects::middleFieldForm::write()
{
    if (writeToFile())
    {
        writeFiles();
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
