/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meanWaveLoads.H"
#include "fvcGrad.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meanWaveLoads, 0);
    addToRunTimeSelectionTable(functionObject, meanWaveLoads, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::meanWaveLoads::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    point origin(Zero);

    // With objectRegistry for access to indirect (global) coordinate systems
    coordSysPtr_ = coordinateSystem::NewIfPresent(obr_, dict);

    if (coordSysPtr_)
    {
        // Report ...
    }
    else if (dict.readIfPresent("CofR", origin))
    {
        const vector e3
        (
            e3Name.empty() ? vector(0, 0, 1) : dict.get<vector>(e3Name)
        );
        const vector e1
        (
            e1Name.empty() ? vector(1, 0, 0) : dict.get<vector>(e1Name)
        );

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        // No 'coordinateSystem' or 'CofR'
        // - enforce a cartesian system

        coordSysPtr_.reset(new coordSystem::cartesian(dict));
    }
}


Foam::volVectorField& Foam::functionObjects::meanWaveLoads::force()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("force"));

    if (!ptr)
    {
        ptr = new volVectorField // create force volVectorField if it doesnt exist
        (
            IOobject
            (
                scopedName("force"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimForce, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


Foam::volVectorField& Foam::functionObjects::meanWaveLoads::moment()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("moment"));

    if (!ptr)
    {
        ptr = new volVectorField // create moment volVectorField 
        (
            IOobject
            (
                scopedName("moment"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimForce*dimLength, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


void Foam::functionObjects::meanWaveLoads::initialise()
{
    if (initialised_)
    {
        return;
    }

    // if // This check was removed cause it failed. Probably because Phi is not present at the start?
    // (
    //     !foundObject<volScalarField>(PhiName_)
    // )
    // {
    //     FatalErrorInFunction
    //         << "Could not find Phi: " << PhiName_
    //         << " in database" << exit(FatalError);
    // }

    if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
    {
        FatalErrorInFunction
            << "Could not find rho:" << rhoName_ << " in database"
            << exit(FatalError);
    }
    initialised_ = true;
}


void Foam::functionObjects::meanWaveLoads::reset()
{
    sumPatchForcesP_ = Zero;
    sumPatchMomentsP_ = Zero;

    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    auto& force = this->force(); // initialize force and moment fields
    auto& moment = this->moment();

    constexpr bool updateAccessTime = false;
    for (const label patchi : patchIDs_)
    {
        force.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
        moment.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
    }
    
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::meanWaveLoads::rho() const // tmp rho volScalarField
{
    if (rhoName_ == "rhoInf")
    {
        return volScalarField::New
        (
            "rho",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_) // set to rhoRef_ since incompressible simulation
        );
    }

    return (lookupObject<volScalarField>(rhoName_));
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::meanWaveLoads::rho(const label patchi) const // rho for each patch
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<scalarField>::New
        (
            mesh_.boundary()[patchi].size(),
            rhoRef_
        );
    }

    const auto& rho = lookupObject<volScalarField>(rhoName_);
    return rho.boundaryField()[patchi];
}


Foam::scalar Foam::functionObjects::meanWaveLoads::rho(const volScalarField& p) const // return rhoRef. Used in calcForcesMoments
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::meanWaveLoads::addToPatchFields // sum all patch contributions. Also asign force to each patch
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fP
)
{
    constexpr bool updateAccessTime = false;

    sumPatchForcesP_ += sum(fP);
    force().boundaryFieldRef(updateAccessTime)[patchi] += fP;

    const vectorField mP(Md^fP);

    sumPatchMomentsP_ += sum(mP);
    moment().boundaryFieldRef(updateAccessTime)[patchi] += mP;
}


void Foam::functionObjects::meanWaveLoads::createIntegratedDataFiles()
{
    if (!forceFilePtr_)
    {
        forceFilePtr_ = newFileAtStartTime("force");
        writeIntegratedDataFileHeader("Force", forceFilePtr_());
    }

    if (!momentFilePtr_)
    {
        momentFilePtr_ = newFileAtStartTime("moment");
        writeIntegratedDataFileHeader("Moment", momentFilePtr_());
    }
}


void Foam::functionObjects::meanWaveLoads::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    const auto vecDesc = [](const word& root)->string
    {
        return root + "_x " + root + "_y " + root + "_z";
    };
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, vecDesc("total"));
    writeTabbed(os, vecDesc("pressure"));


    os  << endl;
}


void Foam::functionObjects::meanWaveLoads::writeIntegratedDataFiles()
{
    const auto& coordSys = coordSysPtr_();

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchForcesP_),
        coordSys.localVector(sumInternalForces_),
        forceFilePtr_()
    );

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchMomentsP_),
        coordSys.localVector(sumInternalMoments_),
        momentFilePtr_()
    );
}


void Foam::functionObjects::meanWaveLoads::writeIntegratedDataFile
(
    const vector& pres,
    const vector& internal,
    OFstream& os
) const
{
    writeCurrentTime(os);

    writeValue(os, pres + internal);
    writeValue(os, pres);


    os  << endl;
}


void Foam::functionObjects::meanWaveLoads::logIntegratedData
(
    const string& descriptor,
    const vector& pres,
    const vector& internal
) const
{
    if (!log)
    {
        return;
    }

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << (pres + internal) << nl // internal should be zero
        << "        Pressure : " << pres << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::meanWaveLoads::meanWaveLoads
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchMomentsP_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    rhoRef_(VGREAT),
    PhiName_("Phi"),
    rhoName_("rho"),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::meanWaveLoads::meanWaveLoads
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchMomentsP_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    rhoRef_(VGREAT),
    PhiName_("Phi"),
    rhoName_("rho"),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::meanWaveLoads::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << ' ' << name() << ':' << endl;

    // Can also use pbm.indices(), but no warnings...
    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();

    // Optional field name entries
    if (dict.readIfPresent<word>("Phi", PhiName_))
    {
        Info<< "    Phi: " << PhiName_ << endl;
    }
    if (dict.readIfPresent<word>("rho", rhoName_))
    {
        Info<< "    rho: " << rhoName_ << endl;
    }
    if (dict.readIfPresent<word>("faceZone", faceZoneName_))
    {
        faceZoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
        if (faceZoneID_ < 0)
        {
            FatalErrorInFunction
                << "faceZone " << faceZoneName_ << " not found" << exit(FatalError);
        }
    }

    dict.readEntry("zeta", zetaName_);
    dict.readEntry("freeSurfacePatch", freeSurfacePatchName_);
    dict.readEntry("cvPoint", cvPoint_);

    freeSurfacePatchID_ =
        mesh_.boundaryMesh().findPatchID(freeSurfacePatchName_);

    if (freeSurfacePatchID_ < 0)
    {
        FatalErrorInFunction
            << "freeSurfacePatch " << freeSurfacePatchName_
            << " not found in boundaryMesh" << exit(FatalError);
    }

    // Reference density needed for incompressible calculations
    if (rhoName_ == "rhoInf")
    {
        rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
        Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
    }

    dict.readEntry("gMag", gMag_);


    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }


    return true;
}


// void Foam::functionObjects::meanWaveLoads::calcForcesMoments()
// {
    // initialise();

    // reset();

    // const point& origin = coordSysPtr_->origin();


    // const auto& Phi = lookupObject<volScalarField>(PhiName_);
    // tmp<volVectorField> tgradPhi = fvc::grad(Phi);
    // const volVectorField& gradPhi = tgradPhi();          // materialize from tmp !! otherwise cant calculate gradient
    // const auto& gradPhib = gradPhi.boundaryField();      // fvPatchVectorField list

    // const auto& Sfb = mesh_.Sf().boundaryField();
    // const auto& Cb = mesh_.C().boundaryField();

    // // const auto& U = lookupObject<volVectorField>(UName_);
    // // tmp<volTensorField> tgradU = fvc::grad(U);
    // // const volTensorField& gradU = tgradU();
    // // const auto& gradUb = gradU.boundaryField();

    // const scalar rhoRef = rho(Phi); // create rho field as Phi?

//     for (const label patchi : patchIDs_)
//     {
//         const vectorField Md(Cb[patchi] - origin);

//         const auto& dphidnSb = gradPhib[patchi] & Sfb[patchi];
//         const vectorField fP(rhoRef*(0.5*Sfb[patchi]*(gradPhib[patchi] & gradPhib[patchi]) - gradPhib[patchi]*dphidnSb));

//         addToPatchFields(patchi, Md, fP);
//     }


//     reduce(sumPatchForcesP_, sumOp<vector>()); // sums per-rank forces into a single global total force vector
//     reduce(sumPatchMomentsP_, sumOp<vector>());
//     reduce(sumInternalForces_, sumOp<vector>());
//     reduce(sumInternalMoments_, sumOp<vector>());
// }

void Foam::functionObjects::meanWaveLoads::calcForcesMoments()
{
    initialise();
    reset();

    const point& origin = coordSysPtr_->origin();

    // Potential and its gradient
    const auto& Phi = lookupObject<volScalarField>(PhiName_);
    // tmp<volVectorField> tgradPhi = fvc::grad(Phi);
    // const volVectorField& gradPhi = tgradPhi();
    // const volVectorField& gradPhi = lookupObject<volScalarField>("U");


    const auto& Sfb = mesh_.Sf().boundaryField();
    const auto& Cb  = mesh_.C().boundaryField();

    const scalar rhoRef = rho(Phi);

    const auto& U = lookupObject<volVectorField>("U");  // or UName_ from dict
    tmp<surfaceVectorField> tUf = fvc::interpolate(U);
    const surfaceVectorField& Uf = tUf();

    const auto& gradPhib = U.boundaryField();   // fvPatchVectorField list


    // ---------------------------------------------------------------------
    // 1. Contributions from selected patches (e.g. hull, if you still want)
    // ---------------------------------------------------------------------
    // for (const label patchi : patchIDs_)
    // {
    //     const vectorField Md(Cb[patchi] - origin);

    //     const auto& gradPhip = gradPhib[patchi];
    //     const auto& Sfp      = Sfb[patchi];

    //     const scalarField dphidnSb(gradPhip & Sfp);

    //     const vectorField fP
    //     (
    //         rhoRef
    //        *(
    //             0.5*Sfp*(gradPhip & gradPhip)
    //           - gradPhip*dphidnSb
    //         )
    //     );

    //     addToPatchFields(patchi, Md, fP);
    // }

    // ---------------------------------------------------------------------
    // 2. Contributions from internal control surface (faceZone)
    // ---------------------------------------------------------------------
    if (faceZoneID_ >= 0)
    {
        const faceZone& fz = mesh_.faceZones()[faceZoneID_];
        const labelList& faces = fz;  // face indices

        const fvMesh& fvm = mesh_;
        const vectorField& Sf = fvm.Sf();
        const vectorField& Cf = fvm.Cf();

        // const labelUList& owner = mesh_.owner();
        // const labelUList& nei   = mesh_.neighbour();

        const point& cvP = cvPoint_;

        forAll(faces, i)
        {
            const label facei = faces[i];

            // Approximate gradPhi on the face
            vector gradPhi_f(Zero);

            if (facei < mesh_.nInternalFaces())
            {
                // const label ownCell = owner[facei];
                // const label neiCell = nei[facei];

                // gradPhi_f = 0.5*(gradPhi[ownCell] + gradPhi[neiCell]);
                gradPhi_f = Uf[facei];
            }
            else
            {
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                const label patchi = pbm.whichPatch(facei);
                const label localFacei = pbm[patchi].whichFace(facei);

                gradPhi_f = gradPhib[patchi][localFacei];
            }

            // Face centre and area vector
            const vector& Cf_f = Cf[facei];
            vector Sf_f = Sf[facei];              // copy, we will flip it

            // Vector from cvPoint to face centre
            const vector d = Cf_f - cvP;

            // If Sf_f points *away* from cvPoint, flip it to point inward
            if ((Sf_f & d) > 0)
            {
                Sf_f = -Sf_f;
            }

            const vector Md = Cf_f - origin;

            const scalar dphidnSb = (gradPhi_f & Sf_f);

            const vector fP
            (
                rhoRef
               *(
                    0.5*Sf_f*(gradPhi_f & gradPhi_f)
                  - gradPhi_f*dphidnSb
                )
            //      rhoRef
            //    *(
            //         -gradPhi_f*dphidnSb
            //     )
            );

            sumPatchForcesP_  += fP;
            sumPatchMomentsP_ += Md ^ fP;
        }
    }

    // ---------------------------------------------------------------------
    // 3. Line integral term: -0.5*rho*g ∫_C n*zeta^2 dl
    //     C = intersection of CV vertical sides (faceZone) with free surface
    // ---------------------------------------------------------------------

    if (faceZoneID_ >= 0)
    {
        const auto& zeta =
            mesh_.lookupObject<volVectorField>(zetaName_);

        const fvMesh& fvm = mesh_;
        const polyMesh& pm = fvm;

        const faceZone& fz = pm.faceZones()[faceZoneID_];

        const edgeList& edges = pm.edges();
        const pointField& pts = pm.points();
        const vectorField& Sf = fvm.Sf();
        const vectorField& Cf = fvm.Cf();

        const labelListList& edgeFaces = pm.edgeFaces();

        const label nFaces = pm.nFaces();
        boolList isCVFace(nFaces, false);
        forAll(fz, i)
        {
            isCVFace[fz[i]] = true;
        }

        const polyBoundaryMesh& pbm = pm.boundaryMesh();
        const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

        // Identify free-surface faces
        boolList isFreeFace(nFaces, false);
        forAll(fsPatch, i)
        {
            const label facei = fsPatch.start() + i;
            isFreeFace[facei] = true;
        }

        // Free-surface elevation on free-surface patch
        const fvPatchVectorField& zetaPatch =
            zeta.boundaryField()[freeSurfacePatchID_];

        const scalar gMag = gMag_; // provided in controlDict

        vector Fzeta(Zero);
        vector Mzeta(Zero);
        // const point& origin = coordSysPtr_->origin();
        const point& cvP = cvPoint_; // control volume point

        forAll(edges, edgei)
        {
            const labelList& eFaces = edgeFaces[edgei]; // faces for this edge

            const label f0 = eFaces[0];
            const label f1 = eFaces[1];

            const bool cv0 = isCVFace[f0];
            const bool cv1 = isCVFace[f1];
            const bool fs0 = isFreeFace[f0];
            const bool fs1 = isFreeFace[f1];

            if (!((cv0 && fs1) || (cv1 && fs0))) continue; // if one face is CV and other is free-surface

            // find which face is on cV and which on free-surface
            const label cvFace = cv0 ? f0 : f1; //(condition) ? (value_if_true) : (value_if_false)
            const label fsFace = fs0 ? f0 : f1;

            // Edge geometry
            const edge& e = edges[edgei];
            const point edgeMid = 0.5*(pts[e[0]] + pts[e[1]]);
            const scalar L = mag(pts[e[1]] - pts[e[0]]);

            // Start from CV face area vector
            vector n = Sf[cvFace];
            const vector& Cf_f = Cf[cvFace];

            // Flip to point inward with respect to cvPoint
            const vector d = Cf_f - cvP;
            if ((n & d) > 0)
            {
                n = -n;
            }

            // Make n horizontal and unit
            n.z() = 0.0;
            const scalar nMagH = mag(n);
            if (nMagH > SMALL)
            {
                n /= nMagH;
            }

            // zeta at this edge: z-component on free-surface face
            const label fsLocalFace = fsFace - fsPatch.start();
            const scalar zetaVal =
                zetaPatch[fsLocalFace].z();

            const scalar coeff = -0.5*rhoRef*gMag*sqr(zetaVal);

            const vector fEdge = coeff*n*L;
            Fzeta += fEdge;

            const vector Md = edgeMid - origin;
            Mzeta += Md ^ fEdge;
        }

        sumPatchForcesP_  += Fzeta;
        sumPatchMomentsP_ += Mzeta;
    }

    // ---------------------------------------------------------------------
    // 3. Parallel reduction
    // ---------------------------------------------------------------------
    reduce(sumPatchForcesP_,   sumOp<vector>());
    reduce(sumPatchMomentsP_,  sumOp<vector>());
    reduce(sumInternalForces_, sumOp<vector>());
    reduce(sumInternalMoments_,sumOp<vector>());
}


Foam::vector Foam::functionObjects::meanWaveLoads::forceEff() const
{
    return sumPatchForcesP_ + sumInternalForces_;
}


Foam::vector Foam::functionObjects::meanWaveLoads::momentEff() const
{
    return sumPatchMomentsP_ + sumInternalMoments_;
}


bool Foam::functionObjects::meanWaveLoads::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    const auto localFp(coordSys.localVector(sumPatchForcesP_));
    const auto localFi(coordSys.localVector(sumInternalForces_));

    logIntegratedData("meanWaveLoads", localFp, localFi);

    const auto localMp(coordSys.localVector(sumPatchMomentsP_));
    const auto localMi(coordSys.localVector(sumInternalMoments_));

    logIntegratedData("moments", localMp, localMi);

    setResult("pressureForce", localFp);
    setResult("internalForce", localFi);
    setResult("pressureMoment", localMp);
    setResult("internalMoment", localMi);

    return true;
}


bool Foam::functionObjects::meanWaveLoads::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment files." << endl;

        createIntegratedDataFiles();
        writeIntegratedDataFiles();
    }

    if (writeFields_)
    {
        Log << "    writing force and moment fields." << endl;

        force().write();
        moment().write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
