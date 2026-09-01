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
#include "processorPolyPatch.H"

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


Foam::tmp<Foam::volScalarField> Foam::functionObjects::meanWaveLoads::rho() const
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
Foam::functionObjects::meanWaveLoads::rho(const label patchi) const
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


Foam::scalar Foam::functionObjects::meanWaveLoads::rho(const volScalarField& p) const
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


void Foam::functionObjects::meanWaveLoads::addToPatchFields
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


void Foam::functionObjects::meanWaveLoads::calcForcesMoments()
{
    initialise();
    reset();

    const point& origin = coordSysPtr_->origin();

    const auto& Phi = lookupObject<volScalarField>(PhiName_);

    const auto& Sfb = mesh_.Sf().boundaryField();
    const auto& Cb  = mesh_.C().boundaryField();

    const scalar rhoRef = rho(Phi);

    const auto& U = lookupObject<volVectorField>("U");  // or UName_ from dict
    const auto& Ucur = lookupObject<volVectorField>("Ucur");  // for forward speed effect in line integral term
    
    tmp<surfaceVectorField> tUf = fvc::interpolate(U);
    const surfaceVectorField& Uf = tUf();
    const auto& Ufb = Uf.boundaryField();

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

        const point& cvP = cvPoint_;

        forAll(faces, i)
        {
            const label facei = faces[i];

            // Avoid double counting on processor patches
            if (Pstream::parRun() && facei >= mesh_.nInternalFaces())
            {
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                const label patchi = pbm.whichPatch(facei);

                if (isA<processorPolyPatch>(pbm[patchi]))
                {
                    const processorPolyPatch& ppp =
                        refCast<const processorPolyPatch>(pbm[patchi]);

                    if (!ppp.owner())
                    {
                        continue; // prevent double counting
                    }
                }
            }
            
            // Approximate velocity vector on the face
            vector gradPhi_f(Zero);

            if (facei < mesh_.nInternalFaces())
            {
                gradPhi_f = Uf[facei];
            }
            else
            {
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                const label patchi = pbm.whichPatch(facei);
                const label localFacei = pbm[patchi].whichFace(facei);

                gradPhi_f = Ufb[patchi][localFacei];
            }

            // Face centre and area vector
            const vector& Cf_f = Cf[facei];
            vector Sf_f = Sf[facei];              // copy, we will flip it

            // Vector from cvPoint to face centre
            const vector d = Cf_f - cvP;

            // If Sf_f points *away* from cvPoint, flip it to point inward
            if ((Sf_f & d) < 0)
            {
                Sf_f = -Sf_f;
            }

            // Reverted back to using Uf & Sf_f
            const scalar dphidnSb = (gradPhi_f & Sf_f);

            const vector Md = Cf_f - origin;

            const vector fP
            (
                rhoRef
               *(
                    0.5*Sf_f*(gradPhi_f & gradPhi_f)
                  - gradPhi_f*dphidnSb 
                )
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
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

        // Make a quick lookup for “is face in CV zone”
        const faceZone& fz = mesh_.faceZones()[faceZoneID_];
        boolList isCVFace(mesh_.nFaces(), false);
        forAll(fz, i) { isCVFace[fz[i]] = true; } // mark CV faces

        // Mark free-surface faces
        boolList isFSFace(mesh_.nFaces(), false);
        forAll(fsPatch, i)
        {
            isFSFace[fsPatch.start() + i] = true;
        }

        const edgeList& edges = mesh_.edges();
        const pointField& pts = mesh_.points();
        const labelList& fsMeshEdges = fsPatch.meshEdges();

        const fvPatchVectorField& zetaPatch = zeta.boundaryField()[freeSurfacePatchID_];

        // For forward speed terms
        const fvPatchVectorField& U_p = U.boundaryField()[freeSurfacePatchID_];
        const fvPatchVectorField& Ucur_p = Ucur.boundaryField()[freeSurfacePatchID_];

        vector Fzeta(Zero);
        vector Mzeta(Zero);

        forAll(fsMeshEdges, ei)
        {
            const label edgei = fsMeshEdges[ei];
            const labelList& eFaces = mesh_.edgeFaces()[edgei];

            // 1. Find the CV face connected to this edge
            label cvFace = -1;
            forAll(eFaces, k)
            {
                if (isCVFace[eFaces[k]])
                {
                    cvFace = eFaces[k];
                    break;
                }
            }
            if (cvFace < 0) continue; // Edge doesn't touch the control volume

            // 2. MPI Ownership Check (Moved to the very top to prevent double counting)
            if (Pstream::parRun() && cvFace >= mesh_.nInternalFaces())
            {
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                const label patchi = pbm.whichPatch(cvFace);

                if (isA<processorPolyPatch>(pbm[patchi]))
                {
                    const processorPolyPatch& ppp = refCast<const processorPolyPatch>(pbm[patchi]);
                    if (!ppp.owner())
                    {
                        continue; // Let the 'owner' processor handle the force for this shared edge.
                    }
                }
            }

            // 3. Establish Edge Geometry
            const edge& e = edges[edgei];
            const point mid = 0.5*(pts[e[0]] + pts[e[1]]);
            const scalar L = mag(pts[e[1]] - pts[e[0]]);

            // 4. Calculate Inward Horizontal Normal
            vector n = mesh_.Sf()[cvFace];
            const vector d = mesh_.Cf()[cvFace] - cvPoint_;
            if ((n & d) < 0) n = -n; // Flip inward wrt cvPoint_
            
            n.z() = 0; // Make horizontal
            const scalar nm = mag(n);
            if (nm < SMALL) continue;
            n /= nm;

            // 5. Gather Stencil for WLS (crawl via edge points)
            DynamicList<label> stencilFaces(8);
            for (const label pointi : {e[0], e[1]})
            {
                const labelList& pFaces = mesh_.pointFaces()[pointi];
                forAll(pFaces, pfi)
                {
                    const label f = pFaces[pfi];
                    if (isFSFace[f] && !stencilFaces.found(f))
                    {
                        stencilFaces.append(f);
                    }
                }
            }

            scalar zetaZ = 0.0;
            vector Un_Uc = vector::zero;
            vector U_Ucn = vector::zero;

            // 6. Local Weighted Least Squares (WLS) extrapolation to edge
            if (stencilFaces.size() >= 3)
            {
                scalar sumW = 0, sumWdx = 0, sumWdy = 0;
                scalar sumWdx2 = 0, sumWdy2 = 0, sumWdxdy = 0;
                scalar sumWz = 0, sumWzdx = 0, sumWzdy = 0;
                
                // Track vector components for U terms
                vector sumW_UnUc = vector::zero, sumW_UnUc_dx = vector::zero, sumW_UnUc_dy = vector::zero;
                vector sumW_UUcn = vector::zero, sumW_UUcn_dx = vector::zero, sumW_UUcn_dy = vector::zero;

                forAll(stencilFaces, sfi)
                {
                    const label f = stencilFaces[sfi];
                    const label localF = f - fsPatch.start();
                    const vector df = mesh_.Cf()[f] - mid;
                    
                    const scalar w = 1.0 / (mag(df) + SMALL); // IDW weighting matrix
                    
                    const scalar z = zetaPatch[localF].z();
                    const vector un_uc = (U_p[localF] & n) * Ucur_p[localF];
                    const vector u_ucn = U_p[localF] * (Ucur_p[localF] & n);

                    sumW     += w;
                    sumWdx   += w * df.x();
                    sumWdy   += w * df.y();
                    sumWdx2  += w * df.x() * df.x();
                    sumWdy2  += w * df.y() * df.y();
                    sumWdxdy += w * df.x() * df.y();
                    
                    sumWz    += w * z;
                    sumWzdx  += w * z * df.x();
                    sumWzdy  += w * z * df.y();

                    sumW_UnUc    += w * un_uc;
                    sumW_UnUc_dx += w * un_uc * df.x();
                    sumW_UnUc_dy += w * un_uc * df.y();

                    sumW_UUcn    += w * u_ucn;
                    sumW_UUcn_dx += w * u_ucn * df.x();
                    sumW_UUcn_dy += w * u_ucn * df.y();
                }

                tensor M(
                    sumW,   sumWdx,   sumWdy,
                    sumWdx, sumWdx2,  sumWdxdy,
                    sumWdy, sumWdxdy, sumWdy2
                );

                if (mag(det(M)) > SMALL)
                {
                    tensor invM = inv(M);
                    
                    // Solve for zeta (scalar)
                    vector b_z(sumWz, sumWzdx, sumWzdy);
                    vector sol_z = invM & b_z;
                    zetaZ = sol_z.x(); // a0 coefficient is the extrapolated value at mid

                    // Solve for velocity terms (component by component)
                    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
                    {
                        vector b_UnUc(sumW_UnUc[cmpt], sumW_UnUc_dx[cmpt], sumW_UnUc_dy[cmpt]);
                        vector sol_UnUc = invM & b_UnUc;
                        Un_Uc[cmpt] = sol_UnUc.x();

                        vector b_UUcn(sumW_UUcn[cmpt], sumW_UUcn_dx[cmpt], sumW_UUcn_dy[cmpt]);
                        vector sol_UUcn = invM & b_UUcn;
                        U_Ucn[cmpt] = sol_UUcn.x();
                    }
                }
                else
                {
                    // Matrix is singular (collinear points), fallback to simple IDW
                    zetaZ = sumWz / (sumW + SMALL);
                    Un_Uc = sumW_UnUc / (sumW + SMALL);
                    U_Ucn = sumW_UUcn / (sumW + SMALL);
                }
            }
            else if (stencilFaces.size() > 0)
            {
                // Edge case: < 3 points available, fallback to pure IDW
                scalar sumW = 0, sumWz = 0;
                vector sumW_UnUc = vector::zero, sumW_UUcn = vector::zero;
                
                forAll(stencilFaces, sfi)
                {
                    const label f = stencilFaces[sfi];
                    const label localF = f - fsPatch.start();
                    const vector df = mesh_.Cf()[f] - mid;
                    const scalar w = 1.0 / (mag(df) + SMALL);
                    
                    sumW += w;
                    sumWz += w * zetaPatch[localF].z();
                    sumW_UnUc += w * ((U_p[localF] & n) * Ucur_p[localF]);
                    sumW_UUcn += w * (U_p[localF] * (Ucur_p[localF] & n));
                }
                
                zetaZ = sumWz / (sumW + SMALL);
                Un_Uc = sumW_UnUc / (sumW + SMALL);
                U_Ucn = sumW_UUcn / (sumW + SMALL);
            }
            else
            {
                continue; // Highly unlikely missing FS face topology
            }

            // 7. Calculate Edge Force and Moment
            const vector coeff = -0.5 * rhoRef * gMag_ * (sqr(zetaZ) * n + 2 * (Un_Uc + U_Ucn) * zetaZ / gMag_);
            const vector fEdge = coeff * L;

            Fzeta += fEdge;
            Mzeta += (mid - origin) ^ fEdge;
        }

        // Add and parallel-reduce as before
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