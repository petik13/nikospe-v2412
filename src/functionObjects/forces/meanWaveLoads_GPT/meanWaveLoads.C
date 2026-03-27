/*---------------------------------------------------------------------------*\\
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
#include "surfaceInterpolate.H"
#include "syncTools.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "DynamicList.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meanWaveLoads, 0);
    addToRunTimeSelectionTable(functionObject, meanWaveLoads, dictionary);
}
}


void Foam::functionObjects::meanWaveLoads::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    point origin(Zero);

    coordSysPtr_ = coordinateSystem::NewIfPresent(obr_, dict);

    if (coordSysPtr_)
    {
        // already set
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
        coordSysPtr_.reset(new coordSystem::cartesian(dict));
    }
}


Foam::volVectorField& Foam::functionObjects::meanWaveLoads::force()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("force"));

    if (!ptr)
    {
        ptr = new volVectorField
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
        ptr = new volVectorField
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

    if (!foundObject<volScalarField>(PhiName_))
    {
        FatalErrorInFunction
            << "Could not find Phi: " << PhiName_
            << " in database" << exit(FatalError);
    }

    if (!foundObject<volVectorField>(UName_))
    {
        FatalErrorInFunction
            << "Could not find U: " << UName_
            << " in database" << exit(FatalError);
    }

    if (!foundObject<volVectorField>(UcurName_))
    {
        FatalErrorInFunction
            << "Could not find Ucur: " << UcurName_
            << " in database" << exit(FatalError);
    }

    if (!foundObject<volVectorField>(zetaName_))
    {
        FatalErrorInFunction
            << "Could not find zeta: " << zetaName_
            << " in database" << exit(FatalError);
    }

    if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
    {
        FatalErrorInFunction
            << "Could not find rho: " << rhoName_
            << " in database" << exit(FatalError);
    }

    initialised_ = true;
}


void Foam::functionObjects::meanWaveLoads::reset()
{
    sumPatchForcesP_ = Zero;
    sumPatchMomentsP_ = Zero;

    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    auto& force = this->force();
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
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return lookupObject<volScalarField>(rhoName_);
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::meanWaveLoads::rho(const label patchi) const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<scalarField>::New(mesh_.boundary()[patchi].size(), rhoRef_);
    }

    const auto& rhoField = lookupObject<volScalarField>(rhoName_);
    return rhoField.boundaryField()[patchi];
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
    os << endl;
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
    os << endl;
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
        << "        Total    : " << (pres + internal) << nl
        << "        Pressure : " << pres << nl;
}


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
    patchIDs_(),
    rhoRef_(VGREAT),
    PhiName_("Phi"),
    UName_("U"),
    UcurName_("Ucur"),
    rhoName_("rho"),
    zetaName_("zeta"),
    freeSurfacePatchName_(word::null),
    freeSurfacePatchID_(-1),
    gMag_(0.0),
    writeFields_(false),
    initialised_(false),
    boxMin_(point(Zero)),
    boxMax_(point(Zero)),
    excludePatches_(),
    excludePatchIDs_()
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
    patchIDs_(),
    rhoRef_(VGREAT),
    PhiName_("Phi"),
    UName_("U"),
    UcurName_("Ucur"),
    rhoName_("rho"),
    zetaName_("zeta"),
    freeSurfacePatchName_(word::null),
    freeSurfacePatchID_(-1),
    gMag_(0.0),
    writeFields_(false),
    initialised_(false),
    boxMin_(point(Zero)),
    boxMax_(point(Zero)),
    excludePatches_(),
    excludePatchIDs_()
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


bool Foam::functionObjects::meanWaveLoads::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;
    excludePatchIDs_.clear();

    Info<< type() << ' ' << name() << ':' << endl;

    patchIDs_ = pbm.patchSet(dict.getOrDefault<wordRes>("patches", wordRes())).sortedToc();

    dict.readIfPresent<word>("Phi", PhiName_);
    dict.readIfPresent<word>("U", UName_);
    dict.readIfPresent<word>("Ucur", UcurName_);
    dict.readIfPresent<word>("rho", rhoName_);

    dict.readEntry("zeta", zetaName_);
    dict.readEntry("freeSurfacePatch", freeSurfacePatchName_);
    dict.readEntry("gMag", gMag_);
    dict.readEntry("boxMin", boxMin_);
    dict.readEntry("boxMax", boxMax_);

    excludePatches_ = dict.getOrDefault<wordRes>("excludePatches", wordRes());

    freeSurfacePatchID_ = pbm.findPatchID(freeSurfacePatchName_);

    if (freeSurfacePatchID_ < 0)
    {
        FatalErrorInFunction
            << "freeSurfacePatch " << freeSurfacePatchName_
            << " not found in boundaryMesh" << exit(FatalError);
    }

    // always exclude the free-surface patch from the side-surface integral
    excludePatchIDs_.insert(freeSurfacePatchID_);

    const labelList excludeIDs = pbm.patchSet(excludePatches_).sortedToc();

    forAll(excludeIDs, i)
    {
        excludePatchIDs_.insert(excludeIDs[i]);
    }

    if (rhoName_ == "rhoInf")
    {
        rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
        Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
    }

    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written for explicit patch contributions only" << endl;
    }

    return true;
}


void Foam::functionObjects::meanWaveLoads::calcForcesMoments()
{
    initialise();
    reset();

    const point& origin = coordSysPtr_->origin();

    const auto& Phi  = lookupObject<volScalarField>(PhiName_);
    const auto& U    = lookupObject<volVectorField>(UName_);
    const auto& Ucur = lookupObject<volVectorField>(UcurName_);
    const auto& zeta = lookupObject<volVectorField>(zetaName_);

    const scalar rhoRef = rho(Phi);

    tmp<surfaceVectorField> tUf = fvc::interpolate(U);
    const surfaceVectorField& Uf = tUf();
    const auto& Ufb = Uf.boundaryField();

    // Build box-membership on cells directly, without creating a temporary
    // field with copied patch types (which can trigger runtime-coded patch
    // machinery such as codeDict).
    const vectorField& Cc = mesh_.C();
    boolList insideCell(mesh_.nCells(), false);

    forAll(insideCell, celli)
    {
        insideCell[celli] = insideBox(Cc[celli]);
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& nei   = mesh_.neighbour();

    // Boundary-face copy of owner-side inside/outside state.
    // After sync, processor/coupled faces contain the neighbouring side value.
    scalarList nbrInside(mesh_.nFaces() - mesh_.nInternalFaces(), 0.0);

    forAll(nbrInside, bFacei)
    {
        const label facei = mesh_.nInternalFaces() + bFacei;
        nbrInside[bFacei] = insideCell[owner[facei]] ? 1.0 : 0.0;
    }

    syncTools::swapBoundaryFaceList(mesh_, nbrInside);

    const vectorField& Sf = mesh_.Sf();
    const vectorField& Cf = mesh_.Cf();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

    const fvPatchVectorField& zetaPatch = zeta.boundaryField()[freeSurfacePatchID_];
    const fvPatchVectorField& Ufs       = U.boundaryField()[freeSurfacePatchID_];
    const fvPatchVectorField& Ucurfs    = Ucur.boundaryField()[freeSurfacePatchID_];

    boolList isFSFace(mesh_.nFaces(), false);
    forAll(fsPatch, i)
    {
        isFSFace[fsPatch.start() + i] = true;
    }

    DynamicList<label> cvFaces;
    DynamicList<vector> cvSfOriented;

    // ---------------------------------------------------------
    // 1) Surface integral over runtime box boundary
    // ---------------------------------------------------------

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        const bool ownInside = insideCell[owner[facei]];
        const bool neiInside = insideCell[nei[facei]];

        if (ownInside == neiInside)
        {
            continue;
        }

        vector Sf_f = ownInside ? Sf[facei] : -Sf[facei];
        const vector gradPhi_f = Uf[facei];
        const vector Md = Cf[facei] - origin;

        const scalar dphidnSb = (gradPhi_f & Sf_f);

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

        cvFaces.append(facei);
        cvSfOriented.append(Sf_f);
    }


    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, localFacei)
        {
            const label facei = pp.start() + localFacei;

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>(pp);

                if (!ppp.owner())
                {
                    continue;
                }

                const bool ownInside = insideCell[pp.faceCells()[localFacei]];
                const label bFacei = facei - mesh_.nInternalFaces();
                const bool neiInside = (nbrInside[bFacei] > 0.5);

                if (ownInside == neiInside)
                {
                    continue;
                }

                vector Sf_f = ownInside ? Sf[facei] : -Sf[facei];
                const vector gradPhi_f = Ufb[patchi][localFacei];
                const vector Md = Cf[facei] - origin;

                const scalar dphidnSb = (gradPhi_f & Sf_f);

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

                cvFaces.append(facei);
                cvSfOriented.append(Sf_f);

                continue;
            }

            if (excludedPatch(patchi))
            {
                continue;
            }

            const label ownCell = pp.faceCells()[localFacei];
            if (!insideCell[ownCell])
            {
                continue;
            }

            vector Sf_f = Sf[facei];
            const vector gradPhi_f = Ufb[patchi][localFacei];
            const vector Md = Cf[facei] - origin;

            const scalar dphidnSb = (gradPhi_f & Sf_f);

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

            cvFaces.append(facei);
            cvSfOriented.append(Sf_f);
        }
    }

    // ---------------------------------------------------------
    // 2) Free-surface line integral derived from selected CV faces
    // ---------------------------------------------------------
    vector Fzeta(Zero);
    vector Mzeta(Zero);
    HashSet<label> doneEdges;

    forAll(cvFaces, iFace)
    {
        const label facei = cvFaces[iFace];
        const vector& Sf_f = cvSfOriented[iFace];

        vector n = Sf_f;
        n.z() = 0.0;

        const scalar nm = mag(n);
        if (nm < SMALL)
        {
            continue;
        }
        n /= nm;

        const labelList& fEdges = mesh_.faceEdges()[facei];

        forAll(fEdges, ie)
        {
            const label edgei = fEdges[ie];

            if (doneEdges.found(edgei))
            {
                continue;
            }

            const labelList& eFaces = mesh_.edgeFaces()[edgei];

            scalar zetaZ = 0.0;
            label nFS = 0;
            vector Un_Uc = vector::zero;
            vector U_Ucn = vector::zero;

            forAll(eFaces, k)
            {
                const label fsFace = eFaces[k];

                if (!isFSFace[fsFace])
                {
                    continue;
                }

                const label fsLocalFace = fsFace - fsPatch.start();

                zetaZ += zetaPatch[fsLocalFace].z();
                Un_Uc += (Ufs[fsLocalFace] & n) * Ucurfs[fsLocalFace];
                U_Ucn += Ufs[fsLocalFace] * (Ucurfs[fsLocalFace] & n);
                ++nFS;
            }

            if (nFS == 0)
            {
                continue;
            }

            doneEdges.insert(edgei);

            zetaZ /= nFS;
            Un_Uc /= nFS;
            U_Ucn /= nFS;

            const edge& e = mesh_.edges()[edgei];
            const point mid = 0.5*(mesh_.points()[e[0]] + mesh_.points()[e[1]]);
            const scalar L = mag(mesh_.points()[e[1]] - mesh_.points()[e[0]]);

            const vector coeff =
                -0.5*rhoRef*gMag_
               *(
                    sqr(zetaZ)*n
                  + 2.0*(Un_Uc + U_Ucn)*zetaZ/gMag_
                );

            const vector fEdge = coeff * L;

            Fzeta += fEdge;
            Mzeta += (mid - origin) ^ fEdge;
        }
    }

    sumPatchForcesP_  += Fzeta;
    sumPatchMomentsP_ += Mzeta;

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
