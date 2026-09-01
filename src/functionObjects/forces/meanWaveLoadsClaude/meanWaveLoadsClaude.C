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

#include "meanWaveLoadsClaude.H"
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
    defineTypeNameAndDebug(meanWaveLoadsClaude, 0);
    addToRunTimeSelectionTable(functionObject, meanWaveLoadsClaude, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::meanWaveLoadsClaude::setCoordinateSystem
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


Foam::volVectorField& Foam::functionObjects::meanWaveLoadsClaude::force()
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


Foam::volVectorField& Foam::functionObjects::meanWaveLoadsClaude::moment()
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


void Foam::functionObjects::meanWaveLoadsClaude::initialise()
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


void Foam::functionObjects::meanWaveLoadsClaude::reset()
{
    sumPatchForcesP_ = Zero;
    sumPatchMomentsP_ = Zero;

    sumSurfaceForces_ = Zero;
    sumLineForces_ = Zero;
    sumSurfaceMoments_ = Zero;
    sumLineMoments_ = Zero;

    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    sumSteadyElevForces_ = Zero;
    sumSteadyElevMoments_ = Zero;

    sumNearFieldForces_ = Zero;
    sumNearFieldMoments_ = Zero;
    sumNearHullForces_ = Zero;
    sumNearWLForces_ = Zero;
    sumLineZetaForces_ = Zero;
    sumLineStripForces_ = Zero;

    auto& force = this->force(); // initialize force and moment fields
    auto& moment = this->moment();

    constexpr bool updateAccessTime = false;
    for (const label patchi : patchIDs_)
    {
        force.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
        moment.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
    }
    
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::meanWaveLoadsClaude::rho() const // tmp rho volScalarField
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
Foam::functionObjects::meanWaveLoadsClaude::rho(const label patchi) const // rho for each patch
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


Foam::scalar Foam::functionObjects::meanWaveLoadsClaude::rho(const volScalarField& p) const // return rhoRef. Used in calcForcesMoments
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


void Foam::functionObjects::meanWaveLoadsClaude::addToPatchFields // sum all patch contributions. Also asign force to each patch
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


void Foam::functionObjects::meanWaveLoadsClaude::createIntegratedDataFiles()
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


void Foam::functionObjects::meanWaveLoadsClaude::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    const auto vecDesc = [](const word& root)->string
    {
        return root + "_x\t" + root + "_y\t" + root + "_z";
    };
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    
    // Create new headers for the columns
    writeTabbed(os, vecDesc("total"));
    writeTabbed(os, vecDesc("surfaceCV"));
    writeTabbed(os, vecDesc("lineFS"));

    os  << endl;
}


void Foam::functionObjects::meanWaveLoadsClaude::writeIntegratedDataFiles()
{
    const auto& coordSys = coordSysPtr_();

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchForcesP_),   // total
        coordSys.localVector(sumSurfaceForces_),  // surface
        coordSys.localVector(sumLineForces_),     // line
        forceFilePtr_()
    );

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchMomentsP_),
        coordSys.localVector(sumSurfaceMoments_),
        coordSys.localVector(sumLineMoments_),
        momentFilePtr_()
    );
}


void Foam::functionObjects::meanWaveLoadsClaude::writeIntegratedDataFile
(
    const vector& total,
    const vector& surface,
    const vector& line,
    OFstream& os
) const
{
    writeCurrentTime(os);

    writeValue(os, total);
    writeValue(os, surface);
    writeValue(os, line);

    os  << endl;
}


void Foam::functionObjects::meanWaveLoadsClaude::logIntegratedData
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

Foam::functionObjects::meanWaveLoadsClaude::meanWaveLoadsClaude
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
    initialised_(false),
    analyticIncident_(true),
    steadyElevationTerms_(true),
    nearField_(true),
    bodyPatchName_("sphere"),
    bodyPatchID_(-1),
    sumNearFieldForces_(Zero),
    sumNearFieldMoments_(Zero),
    sumNearHullForces_(Zero),
    sumNearWLForces_(Zero),
    sumLineZetaForces_(Zero),
    sumLineStripForces_(Zero),
    UDName_("UD"),
    zetaDName_("zetaD"),
    UIName_("UI"),
    waveAmp_(0),
    waveK_(0),
    waveW_(0),
    waveOmegaE_(0),
    waveDepth_(0),
    waveRampTime_(0),
    waveU0_(0),
    waveParamsOK_(false),
    incidentChecked_(false),
    sumSteadyElevForces_(Zero),
    sumSteadyElevMoments_(Zero)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::meanWaveLoadsClaude::meanWaveLoadsClaude
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
    initialised_(false),
    analyticIncident_(true),
    steadyElevationTerms_(true),
    nearField_(true),
    bodyPatchName_("sphere"),
    bodyPatchID_(-1),
    sumNearFieldForces_(Zero),
    sumNearFieldMoments_(Zero),
    sumNearHullForces_(Zero),
    sumNearWLForces_(Zero),
    sumLineZetaForces_(Zero),
    sumLineStripForces_(Zero),
    UDName_("UD"),
    zetaDName_("zetaD"),
    UIName_("UI"),
    waveAmp_(0),
    waveK_(0),
    waveW_(0),
    waveOmegaE_(0),
    waveDepth_(0),
    waveRampTime_(0),
    waveU0_(0),
    waveParamsOK_(false),
    incidentChecked_(false),
    sumSteadyElevForces_(Zero),
    sumSteadyElevMoments_(Zero)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::meanWaveLoadsClaude::read(const dictionary& dict)
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

    // Analytic incident-wave evaluation (default on)
    analyticIncident_ = dict.getOrDefault<Switch>("analyticIncident", Switch(true));
    dict.readIfPresent<word>("UD", UDName_);
    dict.readIfPresent<word>("zetaD", zetaDName_);
    dict.readIfPresent<word>("UI", UIName_);
    Info<< "    analyticIncident: " << analyticIncident_ << endl;

    nearField_ = dict.getOrDefault<Switch>("nearField", Switch(true));
    dict.readIfPresent<word>("bodyPatch", bodyPatchName_);
    if (nearField_)
    {
        bodyPatchID_ = pbm.findPatchID(bodyPatchName_);
        if (bodyPatchID_ < 0)
        {
            WarningInFunction
                << "bodyPatch " << bodyPatchName_ << " not found: near-field"
                   " formulation disabled." << endl;
            nearField_ = false;
        }
    }
    Info<< "    nearField: " << nearField_
        << "  (bodyPatch " << bodyPatchName_ << ")" << endl;

    steadyElevationTerms_ =
        dict.getOrDefault<Switch>("steadyElevationTerms", Switch(true));
    Info<< "    steadyElevationTerms: " << steadyElevationTerms_ << endl;

    // Always read: waveU0_ is needed by the O(Fr^2) group even when the
    // analytic incident evaluation is switched off.
    incidentChecked_ = false;
    readWaveConditions();
    if (!waveParamsOK_)
    {
        analyticIncident_ = false;
        steadyElevationTerms_ = false;
    }


    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }


    return true;
}



namespace
{

//- End points of the segment shared by two faces.
//  Conformal faces share exactly two points.  Across a refinement interface
//  more points can be shared; they are collinear, so the two furthest apart
//  bound the shared segment.  Returns false if fewer than two are shared.
bool sharedSegment
(
    const Foam::face& a,
    const Foam::face& b,
    const Foam::pointField& pts,
    Foam::point& p0,
    Foam::point& p1
)
{
    Foam::labelList shared(a.size());
    Foam::label nShared = 0;

    for (const Foam::label pointi : a)
    {
        if (b.found(pointi))
        {
            shared[nShared++] = pointi;
        }
    }

    if (nShared < 2)
    {
        return false;
    }

    if (nShared == 2)
    {
        p0 = pts[shared[0]];
        p1 = pts[shared[1]];
        return true;
    }

    Foam::scalar dMaxSqr = -1;
    for (Foam::label i = 0; i < nShared; ++i)
    {
        for (Foam::label j = i + 1; j < nShared; ++j)
        {
            const Foam::scalar dSqr =
                Foam::magSqr(pts[shared[i]] - pts[shared[j]]);

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

//- Gravity used by the wave kinematics.  Hard-coded to match shipFlow.C and
//  waveCurrentPotential3DFvPatchScalarField.C, which both use 9.81 when they
//  build PhiI/UI/zetaI.  Using anything else here would put the analytic
//  incident field out of step with the one the solver actually imposed.
constexpr Foam::scalar gWave_ = 9.81;

} // End anonymous namespace


void Foam::functionObjects::meanWaveLoadsClaude::readWaveConditions()
{
    waveParamsOK_ = false;

    IOobject wcIO
    (
        "waveCurConditions",
        time_.constant(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!wcIO.typeHeaderOk<IOdictionary>(true))
    {
        WarningInFunction
            << "constant/waveCurConditions not found: analytic incident"
               " evaluation disabled, falling back to the discrete fields."
            << endl;
        return;
    }

    IOdictionary wc(wcIO);

    const scalar steepness  = wc.getOrDefault<scalar>("steepness", 0.0);
    const scalar wavelength = wc.getOrDefault<scalar>("wavelength", 0.0);
    const scalar U0         = wc.getOrDefault<scalar>("currentspeed", 0.0);
    const scalar headAng    = wc.getOrDefault<scalar>("head_ang", 0.0);
    const scalar rampperiod = wc.getOrDefault<scalar>("rampperiod", 0.0);
    waveDepth_ = wc.getOrDefault<scalar>("waterdepth", 0.0);

    if (wavelength <= SMALL || waveDepth_ <= SMALL)
    {
        WarningInFunction
            << "Invalid wavelength/waterdepth in waveCurConditions:"
               " analytic incident evaluation disabled." << endl;
        return;
    }

    waveU0_  = U0;
    waveAmp_ = 0.5*steepness*wavelength;
    waveK_   = 2.0*constant::mathematical::pi/wavelength;
    waveW_   = Foam::sqrt(gWave_*waveK_*Foam::tanh(waveK_*waveDepth_));
    waveOmegaE_ = waveW_ + waveK_*U0*Foam::cos(headAng);
    waveRampTime_ = rampperiod*2.0*constant::mathematical::pi/waveW_;

    waveParamsOK_ = true;

    Info<< "    incident wave: amp = " << waveAmp_
        << "  k = " << waveK_
        << "  w = " << waveW_
        << "  w_e = " << waveOmegaE_
        << "  h = " << waveDepth_
        << "  rampTime = " << waveRampTime_ << endl;
}


Foam::scalar Foam::functionObjects::meanWaveLoadsClaude::incidentZeta
(
    const point& p,
    const scalar t
) const
{
    if (!waveParamsOK_) return 0.0;

    const scalar ramp =
        waveRampTime_ > SMALL
      ? 0.5*(1.0 - Foam::cos(constant::mathematical::pi*min(1.0, t/waveRampTime_)))
      : 1.0;

    return waveAmp_*Foam::cos(waveK_*p.x() - waveOmegaE_*t)*ramp;
}


Foam::vector Foam::functionObjects::meanWaveLoadsClaude::incidentU
(
    const point& p,
    const scalar t
) const
{
    if (!waveParamsOK_) return vector::zero;

    const scalar ramp =
        waveRampTime_ > SMALL
      ? 0.5*(1.0 - Foam::cos(constant::mathematical::pi*min(1.0, t/waveRampTime_)))
      : 1.0;

    const scalar phase   = waveK_*p.x() - waveOmegaE_*t;
    const scalar chDenom = Foam::cosh(waveK_*waveDepth_);
    const scalar coshT   = Foam::cosh(waveK_*(waveDepth_ + p.z()))/chDenom;
    const scalar sinhT   = Foam::sinh(waveK_*(waveDepth_ + p.z()))/chDenom;
    const scalar Ua      = waveAmp_*waveK_*gWave_/waveW_*ramp;

    return vector(Ua*coshT*Foam::cos(phase), 0.0, Ua*sinhT*Foam::sin(phase));
}


void Foam::functionObjects::meanWaveLoadsClaude::checkIncident()
{
    if (incidentChecked_ || !analyticIncident_ || !waveParamsOK_) return;
    incidentChecked_ = true;

    if (!mesh_.foundObject<volVectorField>(UIName_))
    {
        WarningInFunction
            << "Field " << UIName_ << " not found: cannot verify that the"
               " analytic incident wave matches the one the solver imposed."
            << endl;
        return;
    }

    const auto& UI = mesh_.lookupObject<volVectorField>(UIName_);
    const scalar t = time_.value();
    const vectorField& C = mesh_.C();

    scalar maxErr = 0.0;
    scalar maxRef = SMALL;
    forAll(C, celli)
    {
        maxErr = max(maxErr, mag(incidentU(C[celli], t) - UI[celli]));
        maxRef = max(maxRef, mag(UI[celli]));
    }
    reduce(maxErr, maxOp<scalar>());
    reduce(maxRef, maxOp<scalar>());

    const scalar relErr = maxErr/maxRef;

    Info<< "    incident consistency: max|UI_analytic - " << UIName_ << "| = "
        << maxErr << "  (" << 100*relErr << "% of max|UI|)" << endl;

    if (relErr > 1e-6)
    {
        WarningInFunction
            << "Analytic incident wave differs from the solver's " << UIName_
            << " field by " << 100*relErr << "%.  The wave parameters or the"
               " formulae are out of step with the solver; analytic incident"
               " evaluation will bias the loads." << endl;
    }
}


void Foam::functionObjects::meanWaveLoadsClaude::calcNearField
(
    const volVectorField& U,
    const scalar tNow
)
{
    // ---------------------------------------------------------------------
    // Near-field (direct pressure integration), Chen (2022) Eq.(29), for a
    // FIXED wall-sided body with psi_bar dropped:
    //
    //     F_N = rho * int_hull 1/2 <|u|^2> n dS  -  1/2 rho g int_WL <zeta^2> n dl
    //
    // Chen's n is "positively into the fluid", i.e. OUT of the body, which is
    // -Sf/|Sf| for the hull patch (OpenFOAM's Sf points out of the fluid).
    //
    // Sign check on the waterline term: run-up is larger on the upstream face,
    // where n_out_body is -x, so -1/2 rho g zeta^2 n gives +x -- the drift
    // pushes downstream, as it must.  The hull term is a suction where |u| is
    // large and partly cancels it; that near-cancellation is the whole reason
    // this is worth comparing against the midfield route.
    //
    // The waterline is found from cells owning BOTH a hull face and a
    // free-surface face; each cell belongs to exactly one rank, so this is
    // parallel-exact with no communication and no double counting.
    // ---------------------------------------------------------------------
    const point& origin = coordSysPtr_->origin();
    const scalar rhoRef = rho(lookupObject<volScalarField>(PhiName_));

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& bp = pbm[bodyPatchID_];
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();
    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const labelList& faceOwner = mesh_.faceOwner();

    const fvPatchVectorField& Ub = U.boundaryField()[bodyPatchID_];

    // ---- hull term:  rho * sum 1/2 |u|^2 * (-Sf) --------------------------
    vector Fh(Zero), Mh(Zero);
    forAll(bp, i)
    {
        const label facei = bp.start() + i;

        vector u = Ub[i];
        if (analyticIncident_)
        {
            u += incidentU(faceCentres[facei], tNow);
        }

        const vector fH(-rhoRef*0.5*magSqr(u)*faceAreas[facei]);
        Fh += fH;
        Mh += (faceCentres[facei] - origin) ^ fH;
    }

    // ---- waterline term:  -1/2 rho g sum <zeta^2> L n_out ----------------
    const auto& zeta =
        analyticIncident_
      ? mesh_.lookupObject<volVectorField>(zetaDName_)
      : mesh_.lookupObject<volVectorField>(zetaName_);
    const fvPatchVectorField& zetaPatch =
        zeta.boundaryField()[freeSurfacePatchID_];

    const label fsStart = fsPatch.start();
    const label fsEnd = fsStart + fsPatch.size();
    const label bStart = bp.start();
    const cellList& cells = mesh_.cells();

    vector Fw(Zero), Mw(Zero);
    forAll(bp, i)
    {
        const label bodyFacei = bStart + i;
        const label celli = faceOwner[bodyFacei];

        for (const label facej : cells[celli])
        {
            if (facej < fsStart || facej >= fsEnd) continue;   // free-surface face?

            point p0(Zero), p1(Zero);
            if (!sharedSegment(faces[facej], faces[bodyFacei], pts, p0, p1)) continue;

            const scalar L = mag(p1 - p0);
            if (L < SMALL) continue;

            const point mid(0.5*(p0 + p1));

            // outward horizontal normal OF THE BODY  =  -Sf/|Sf|, z removed
            vector n = -faceAreas[bodyFacei];
            n.z() = 0;
            const scalar nm = mag(n);
            if (nm < SMALL) continue;
            n /= nm;

            scalar zetaZ = zetaPatch[facej - fsStart].z();
            if (analyticIncident_)
            {
                zetaZ += incidentZeta(mid, tNow);
            }

            const vector fW(-0.5*rhoRef*gMag_*sqr(zetaZ)*L*n);
            Fw += fW;
            Mw += (mid - origin) ^ fW;
        }
    }

    sumNearHullForces_ += Fh;
    sumNearWLForces_   += Fw;
    sumNearFieldForces_  += Fh + Fw;
    sumNearFieldMoments_ += Mh + Mw;
}


void Foam::functionObjects::meanWaveLoadsClaude::writeNearField()
{
    if (!Pstream::master()) return;

    if (!nearFilePtr_)
    {
        fileName dir("./postProcessing/" + name());
        mkDir(dir);
        nearFilePtr_.reset(new OFstream(dir/"forceNear.dat"));
        *nearFilePtr_
            << "# Near-field (Chen Eq.29) drift force, fixed wall-sided body" << nl
            << "# CofR : " << coordSysPtr_->origin() << nl
            << "# Time\ttotal_x\ttotal_y\ttotal_z"
               "\thull_x\thull_y\thull_z\twaterline_x\twaterline_y\twaterline_z"
            << endl;
    }

    *nearFilePtr_
        << time_.timeOutputValue() << '\t'
        << sumNearFieldForces_.x() << '\t' << sumNearFieldForces_.y() << '\t'
        << sumNearFieldForces_.z() << '\t'
        << sumNearHullForces_.x() << '\t' << sumNearHullForces_.y() << '\t'
        << sumNearHullForces_.z() << '\t'
        << sumNearWLForces_.x() << '\t' << sumNearWLForces_.y() << '\t'
        << sumNearWLForces_.z() << endl;
}


void Foam::functionObjects::meanWaveLoadsClaude::calcForcesMoments()
{
    initialise();
    reset();

    const point& origin = coordSysPtr_->origin();

    const auto& Phi = lookupObject<volScalarField>(PhiName_);

    const scalar rhoRef = rho(Phi);

    const scalar tNow = time_.value();

    // With analyticIncident_ the discrete fields supply only the SCATTERED
    // part and the incident part is added in closed form at the exact
    // evaluation point.  Otherwise fall back to the total discrete fields,
    // which reproduces the original behaviour.
    if (analyticIncident_)
    {
        if
        (
            !mesh_.foundObject<volVectorField>(UDName_)
         || !mesh_.foundObject<volVectorField>(zetaDName_)
        )
        {
            WarningInFunction
                << "Fields " << UDName_ << '/' << zetaDName_
                << " not found: disabling analytic incident evaluation."
                << endl;
            analyticIncident_ = false;
        }
    }
    checkIncident();

    const auto& U =
        analyticIncident_
      ? lookupObject<volVectorField>(UDName_)
      : lookupObject<volVectorField>("U");

    const auto& Ucur = lookupObject<volVectorField>("Ucur");  // for forward speed effect in line integral term
    tmp<surfaceVectorField> tUf = fvc::interpolate(U);
    const surfaceVectorField& Uf = tUf();
    const auto& Ufb = Uf.boundaryField();

    // const auto& gradPhib = U.boundaryField();   // fvPatchVectorField list


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

        // NOTE: primitiveMesh::faceAreas()/faceCentres() are sized nFaces, so
        // they are valid for the CV faces that end up on a processor patch
        // after decomposition.  mesh_.Sf()/mesh_.Cf() only hold nInternalFaces
        // values in their internal field, so indexing those with a boundary
        // face label reads past the end of the array (silently, in an Opt
        // build).  In serial the whole faceZone is internal and it never
        // shows; in parallel it corrupts every CV face on a rank boundary.
        const vectorField& Sf = mesh_.faceAreas();
        const vectorField& Cf = mesh_.faceCentres();

        // const labelUList& owner = mesh_.owner();
        // const labelUList& nei   = mesh_.neighbour();

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

                gradPhi_f = Ufb[patchi][localFacei];
                // Pout << "Face " << facei << " is a boundary face on patch " << patchi
                //     << " with local index " << localFacei
                //     << " and gradPhi_f = " << gradPhi_f << endl;
            }

            // Face centre and area vector
            const vector& Cf_f = Cf[facei];

            // Incident part evaluated exactly at the face centre rather than
            // interpolated from cell centres.  Interpolation attenuates it by
            // cos(k*dx/2); since the incident field is what has to cancel
            // between this integral and the waterline one, that attenuation
            // lands straight on the (much smaller) difference.
            if (analyticIncident_)
            {
                gradPhi_f += incidentU(Cf_f, tNow);
            }

            vector Sf_f = Sf[facei];              // copy, we will flip it

            // Vector from cvPoint to face centre
            const vector d = Cf_f - cvP;

            // If Sf_f points *away* from cvPoint, flip it to point inward
            if ((Sf_f & d) < 0)
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

            sumSurfaceForces_  += fP;
            sumSurfaceMoments_ += Md ^ fP;
        }
    }

    // ---------------------------------------------------------------------
    // 3. Line integral term along the waterline C of the control volume
    //      C = intersection of the CV vertical sides with the free surface
    //
    //    Driven from the CV faces rather than from free-surface mesh edges.
    //    Each CV side face contributes the segment it shares with the free
    //    surface, taken with weight 0.5 from each of the two cells either
    //    side of it.  Every cell belongs to exactly one rank, so an internal
    //    CV face collects 0.5 + 0.5 locally while a CV face lying on a
    //    processor patch collects 0.5 from each rank.  That is exact for any
    //    decomposition, needs no communication and cannot double count.
    //
    //    The previous edge-driven version keyed off mesh_.edgeFaces(), which
    //    only sees faces local to the rank, so at a partition boundary the
    //    zeta/U average was one-sided and the "first CV face in eFaces" pick
    //    was not reproducible between serial and parallel runs.
    //
    //    zeta and u are taken per segment, so the quadratic term accumulates
    //    <zeta^2> (not <zeta>^2) and the O(U) term accumulates the product
    //    <zeta*u> -- both without any pre-averaging.
    // ---------------------------------------------------------------------

    if (faceZoneID_ >= 0)
    {
        const auto& zeta =
            analyticIncident_
          ? mesh_.lookupObject<volVectorField>(zetaDName_)
          : mesh_.lookupObject<volVectorField>(zetaName_);

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const polyPatch& fsPatch = pbm[freeSurfacePatchID_];
        const label fsStart = fsPatch.start();
        const label fsEnd = fsStart + fsPatch.size();

        const faceZone& fz = mesh_.faceZones()[faceZoneID_];

        const faceList& faces = mesh_.faces();
        const pointField& pts = mesh_.points();
        const cellList& cells = mesh_.cells();
        const labelList& faceOwner = mesh_.faceOwner();
        const labelList& faceNeighbour = mesh_.faceNeighbour();

        // Sized nFaces, so valid for CV faces sitting on a processor patch
        const vectorField& faceAreas = mesh_.faceAreas();
        const vectorField& faceCentres = mesh_.faceCentres();

        const fvPatchVectorField& zetaPatch =
            zeta.boundaryField()[freeSurfacePatchID_];
        const fvPatchVectorField& U_p = U.boundaryField()[freeSurfacePatchID_];
        const fvPatchVectorField& Ucur_p =
            Ucur.boundaryField()[freeSurfacePatchID_];

        vector Fzeta(Zero);
        vector Mzeta(Zero);

        forAll(fz, i)
        {
            const label cvFacei = fz[i];

            vector n = faceAreas[cvFacei];

            // Orient outward wrt cvPoint_
            if ((n & (faceCentres[cvFacei] - cvPoint_)) < 0) n = -n;

            // Make horizontal and unit
            n.z() = 0;
            const scalar nm = mag(n);
            if (nm < SMALL) continue;   // horizontal CV face: no line contribution
            n /= nm;

            const face& cvFace = faces[cvFacei];

            // Cells this rank owns either side of the CV face.  A face on a
            // processor patch has no local neighbour - the adjoining rank
            // supplies that half.
            label sideCells[2] = {faceOwner[cvFacei], -1};
            if (cvFacei < mesh_.nInternalFaces())
            {
                sideCells[1] = faceNeighbour[cvFacei];
            }

            for (int side = 0; side < 2; ++side)
            {
                const label celli = sideCells[side];
                if (celli < 0) continue;

                for (const label facej : cells[celli])
                {
                    // Only a free-surface face of this cell can meet C
                    if (facej < fsStart || facej >= fsEnd) continue;

                    point p0(Zero);
                    point p1(Zero);
                    if (!sharedSegment(faces[facej], cvFace, pts, p0, p1)) continue;

                    const scalar L = mag(p1 - p0);
                    if (L < SMALL) continue;

                    const point mid(0.5*(p0 + p1));
                    const label fsLocalFace = facej - fsStart;

                    // Free-surface state on this segment.  The scattered part
                    // comes from the discrete solution on the adjacent face;
                    // the incident part is closed form, so evaluate it at the
                    // segment midpoint instead of inheriting the face-centre
                    // value.
                    scalar zetaZ = zetaPatch[fsLocalFace].z();
                    vector uFS = U_p[fsLocalFace];

                    if (analyticIncident_)
                    {
                        zetaZ += incidentZeta(mid, tNow);
                        uFS   += incidentU(mid, tNow);
                    }

                    const vector Un_Uc((uFS & n)*Ucur_p[fsLocalFace]);
                    const vector U_Ucn(uFS*(Ucur_p[fsLocalFace] & n));

                    //   -1/2 rho g zeta^2 n            : wave elevation
                    //   -rho zeta [ (u.n)W + u(W.n) ]  : strip flux, O(U)
                    // Kept apart so the O(U) part can be sized on its own.
                    const vector coeffZ(-0.5*rhoRef*gMag_*sqr(zetaZ)*n);
                    const vector coeffS(-rhoRef*(Un_Uc + U_Ucn)*zetaZ);

                    // 0.5: this cell is one of the two sides of the CV face
                    const vector fEdgeZ(0.5*coeffZ*L);
                    const vector fEdgeS(0.5*coeffS*L);
                    const vector fEdge(fEdgeZ + fEdgeS);

                    sumLineZetaForces_  += fEdgeZ;
                    sumLineStripForces_ += fEdgeS;

                    Fzeta += fEdge;
                    Mzeta += (mid - origin) ^ fEdge;

                    // --------------------------------------------------------
                    // Chen (2022), Eq.(46): the O(Fr^2) waterline group
                    //
                    //   F^W_M = closed_int { eta_0 [ 2 eta_2 n + phi_n grad_h phi ]
                    //                      + eta_2 (w.n) grad_h phi^s } dl
                    //
                    // with (his Eq.18)
                    //   eta_0 = -1/2 (w.w - 1)  ->  (U_inf^2 - |W|^2)/(2g)
                    //   eta_2 = -1/2 grad(phi).grad(phi) - Fr w.grad(psi_bar)
                    //        ->  -|u|^2/(2g)      (psi_bar dropped, as in Chen
                    //                              section V and refs [5],[7])
                    //
                    // Chen's n points INTO the fluid, i.e. towards the body on
                    // the control surface, so every term flips sign against the
                    // outward n used here.  Each of the three carries two powers
                    // of the steady flow: the group is O(eps^2 U^2) and even in
                    // U, and it vanishes identically at U = 0.
                    //
                    // eta_2 is formed instantaneously; since eta_0 and W are
                    // steady, time-averaging the product afterwards recovers
                    // <eta_2> correctly, so no second-order BVP is needed.
                    // --------------------------------------------------------
                    if (steadyElevationTerms_)
                    {
                        const vector& W = Ucur_p[fsLocalFace];

                        const scalar etaS =
                            (sqr(waveU0_) - magSqr(W))/(2*gMag_);   // eta_0
                        const scalar etaBar = -magSqr(uFS)/(2*gMag_); // eta_2

                        vector Wh(W);    Wh.z() = 0;
                        vector uh(uFS);  uh.z() = 0;

                        const vector fW
                        (
                           -rhoRef
                           *(
                                2*gMag_*etaS*etaBar*n
                              + etaS*(uFS & n)*uh
                              + etaBar*(W & n)*Wh
                            )
                        );

                        const vector fEdgeW(0.5*fW*L);

                        sumSteadyElevForces_  += fEdgeW;
                        sumSteadyElevMoments_ += (mid - origin) ^ fEdgeW;
                    }
                }
            }
        }

        sumLineForces_  += Fzeta;
        sumLineMoments_ += Mzeta;
    }

    // ---------------------------------------------------------------------
    // 4. Parallel reduction and Total Summation
    // ---------------------------------------------------------------------
    // Reduce the individual components across all processors
    reduce(sumSurfaceForces_,  sumOp<vector>());
    reduce(sumSurfaceMoments_, sumOp<vector>());
    reduce(sumLineForces_,     sumOp<vector>());
    reduce(sumLineMoments_,    sumOp<vector>());

    reduce(sumSteadyElevForces_,  sumOp<vector>());
    reduce(sumSteadyElevMoments_, sumOp<vector>());

    reduce(sumLineZetaForces_,  sumOp<vector>());
    reduce(sumLineStripForces_, sumOp<vector>());

    if (steadyElevationTerms_)
    {
        Log << "    Chen Eq.(46) O(Fr^2) waterline group: "
            << sumSteadyElevForces_ << endl;

        sumLineForces_  += sumSteadyElevForces_;
        sumLineMoments_ += sumSteadyElevMoments_;
    }

    reduce(sumInternalForces_, sumOp<vector>());
    reduce(sumInternalMoments_,sumOp<vector>());

    // Reconstruct the totals for standard OpenFOAM logging
    sumPatchForcesP_  = sumSurfaceForces_ + sumLineForces_;
    sumPatchMomentsP_ = sumSurfaceMoments_ + sumLineMoments_;

    // ---------------------------------------------------------------------
    // 5. Independent near-field route over the hull (Chen Eq.29)
    // ---------------------------------------------------------------------
    if (nearField_)
    {
        calcNearField(U, tNow);

        reduce(sumNearHullForces_,    sumOp<vector>());
        reduce(sumNearWLForces_,      sumOp<vector>());
        reduce(sumNearFieldForces_,   sumOp<vector>());
        reduce(sumNearFieldMoments_,  sumOp<vector>());

        Log << "    near-field  total " << sumNearFieldForces_ << nl
            << "                hull  " << sumNearHullForces_
            << "  waterline " << sumNearWLForces_ << nl
            << "    midfield    total " << sumPatchForcesP_ << endl;

        writeNearField();
    }

    Log << "    midfield split: surfaceCV " << sumSurfaceForces_
        << "  lineZeta " << sumLineZetaForces_
        << "  lineStrip(O(U)) " << sumLineStripForces_
        << "  steadyElev(O(U^2)) " << sumSteadyElevForces_ << endl;

    writeMidBreakdown();
}


void Foam::functionObjects::meanWaveLoadsClaude::writeMidBreakdown()
{
    if (!Pstream::master()) return;

    if (!midFilePtr_)
    {
        fileName dir("./postProcessing/" + name());
        mkDir(dir);
        midFilePtr_.reset(new OFstream(dir/"forceMid.dat"));
        *midFilePtr_
            << "# Midfield (Chen Eq.43) drift force, split by contribution" << nl
            << "# CofR : " << coordSysPtr_->origin() << nl
            << "# Time\ttotal_x\ttotal_y\ttotal_z"
               "\tsurfaceCV_x\tsurfaceCV_y\tsurfaceCV_z"
               "\tlineZeta_x\tlineZeta_y\tlineZeta_z"
               "\tlineStrip_x\tlineStrip_y\tlineStrip_z"
               "\tsteadyElev_x\tsteadyElev_y\tsteadyElev_z" << nl
            << "# lineStrip is O(U); steadyElev is O(U^2); both vanish at U = 0"
            << endl;
    }

    const vector& a = sumPatchForcesP_;
    const vector& b = sumSurfaceForces_;
    const vector& c = sumLineZetaForces_;
    const vector& d = sumLineStripForces_;
    const vector& e = sumSteadyElevForces_;

    *midFilePtr_
        << time_.timeOutputValue() << '\t'
        << a.x() << '\t' << a.y() << '\t' << a.z() << '\t'
        << b.x() << '\t' << b.y() << '\t' << b.z() << '\t'
        << c.x() << '\t' << c.y() << '\t' << c.z() << '\t'
        << d.x() << '\t' << d.y() << '\t' << d.z() << '\t'
        << e.x() << '\t' << e.y() << '\t' << e.z() << endl;
}


Foam::vector Foam::functionObjects::meanWaveLoadsClaude::forceEff() const
{
    return sumPatchForcesP_ + sumInternalForces_;
}


Foam::vector Foam::functionObjects::meanWaveLoadsClaude::momentEff() const
{
    return sumPatchMomentsP_ + sumInternalMoments_;
}


bool Foam::functionObjects::meanWaveLoadsClaude::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    const auto localFp(coordSys.localVector(sumPatchForcesP_));
    const auto localFi(coordSys.localVector(sumInternalForces_));

    logIntegratedData("meanWaveLoadsClaude", localFp, localFi);

    const auto localMp(coordSys.localVector(sumPatchMomentsP_));
    const auto localMi(coordSys.localVector(sumInternalMoments_));

    logIntegratedData("moments", localMp, localMi);

    setResult("pressureForce", localFp);
    setResult("internalForce", localFi);
    setResult("pressureMoment", localMp);
    setResult("internalMoment", localMi);

    return true;
}


bool Foam::functionObjects::meanWaveLoadsClaude::write()
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
