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

#include "linBodyMotionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

using vector6 = Foam::FixedList<Foam::scalar, 6>;
using tensor6 = Foam::FixedList<vector6, 6>;

inline vector6 zero6()
{
    vector6 v;
    for (int i = 0; i < 6; ++i) v[i] = 0;
    return v;
}

inline tensor6 zero66()
{
    tensor6 A;
    for (int i = 0; i < 6; ++i) A[i] = zero6();
    return A;
}

inline vector6 operator+(const vector6& a, const vector6& b)
{
    vector6 r;
    for (int i = 0; i < 6; ++i) r[i] = a[i] + b[i];
    return r;
}

inline vector6 operator-(const vector6& a, const vector6& b)
{
    vector6 r;
    for (int i = 0; i < 6; ++i) r[i] = a[i] - b[i];
    return r;
}

inline vector6 operator*(const Foam::scalar s, const vector6& a)
{
    vector6 r;
    for (int i = 0; i < 6; ++i) r[i] = s*a[i];
    return r;
}

inline vector6 mult6(const tensor6& A, const vector6& x)
{
    vector6 r = zero6();
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            r[i] += A[i][j]*x[j];
    return r;
}

//- Gaussian elimination with partial pivoting
vector6 solve6(tensor6 A, vector6 b)
{
    for (int k = 0; k < 6; ++k)
    {
        int piv = k;
        Foam::scalar amax = Foam::mag(A[k][k]);
        for (int i = k + 1; i < 6; ++i)
        {
            const Foam::scalar v = Foam::mag(A[i][k]);
            if (v > amax) { amax = v; piv = i; }
        }

        if (amax < Foam::VSMALL)
        {
            FatalErrorInFunction
                << "Singular 6x6 system at row " << k
                << ": check the mass, damping and restoring matrices."
                << Foam::abort(Foam::FatalError);
        }

        if (piv != k)
        {
            for (int j = k; j < 6; ++j) std::swap(A[k][j], A[piv][j]);
            std::swap(b[k], b[piv]);
        }

        for (int i = k + 1; i < 6; ++i)
        {
            const Foam::scalar f = A[i][k]/A[k][k];
            A[i][k] = 0;
            for (int j = k + 1; j < 6; ++j) A[i][j] -= f*A[k][j];
            b[i] -= f*b[k];
        }
    }

    vector6 x = zero6();
    for (int i = 5; i >= 0; --i)
    {
        Foam::scalar s = b[i];
        for (int j = i + 1; j < 6; ++j) s -= A[i][j]*x[j];
        x[i] = s/A[i][i];
    }
    return x;
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linBodyMotionFvPatchScalarField::
linBodyMotionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(p, iF),
    PhiIName_("PhiI"),
    UsName_("Us"),
    gradUsName_("gradUs"),
    pSName_("pS"),
    pName_("p"),
    fixBody_(false),
    mass_(1), xG_(0), zG_(0),
    I44_(1), I55_(1), I66_(1), I46_(0),
    C33_(0), C44_(0), C55_(0), C35_(0),
    rollDampingFraction_(0),
    mooringK11_(0), mooringK22_(0), mooringK66_(0),
    mooringDampingFraction_(0),
    heading_(0),
    addedMass_(scalar(0)),
    addedMassFromTimeStep_(true),
    accelRelaxation_(1),
    accelDamping_(1),
    gamma_(0.5), beta_(0.25),
    waveAmp_(0), waveNumber_(0), omega_(0), omegaE_(0), depth_(0), rampTime_(0),
    Uinf_(Zero),
    rho_(1000),
    steadyRestoring_(true),
    Ks_(zero66()), KsValid_(false),
    q_(zero6()), qDot_(zero6()), qDdot_(zero6()),
    qNew_(zero6()), qDotNew_(zero6()), qDdotNew_(zero6()),
    qDdotPrevIter_(zero6()),
    lastTimeIndex_(-1)
{
    gradient() = Zero;
    readProperties();
}


Foam::linBodyMotionFvPatchScalarField::
linBodyMotionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict, IOobjectOption::NO_READ),
    PhiIName_(dict.getOrDefault<word>("PhiI", "PhiI")),
    UsName_(dict.getOrDefault<word>("Us", "Us")),
    gradUsName_(dict.getOrDefault<word>("gradUs", "gradUs")),
    pSName_(dict.getOrDefault<word>("pS", "pS")),
    pName_(dict.getOrDefault<word>("p", "p")),
    fixBody_(false),
    mass_(1), xG_(0), zG_(0),
    I44_(1), I55_(1), I66_(1), I46_(0),
    C33_(0), C44_(0), C55_(0), C35_(0),
    rollDampingFraction_(0),
    mooringK11_(0), mooringK22_(0), mooringK66_(0),
    mooringDampingFraction_(0),
    heading_(0),
    addedMass_(scalar(0)),
    addedMassFromTimeStep_(true),
    accelRelaxation_(1),
    accelDamping_(1),
    gamma_(0.5), beta_(0.25),
    waveAmp_(0), waveNumber_(0), omega_(0), omegaE_(0), depth_(0), rampTime_(0),
    Uinf_(Zero),
    rho_(1000),
    steadyRestoring_(true),
    Ks_(zero66()), KsValid_(false),
    q_(zero6()), qDot_(zero6()), qDdot_(zero6()),
    qNew_(zero6()), qDotNew_(zero6()), qDdotNew_(zero6()),
    qDdotPrevIter_(zero6()),
    lastTimeIndex_(-1)
{
    gradient() = Zero;
    readProperties();
}


Foam::linBodyMotionFvPatchScalarField::
linBodyMotionFvPatchScalarField
(
    const linBodyMotionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    PhiIName_(ptf.PhiIName_),
    UsName_(ptf.UsName_),
    gradUsName_(ptf.gradUsName_),
    pSName_(ptf.pSName_),
    pName_(ptf.pName_),
    fixBody_(ptf.fixBody_),
    mass_(ptf.mass_), xG_(ptf.xG_), zG_(ptf.zG_),
    I44_(ptf.I44_), I55_(ptf.I55_), I66_(ptf.I66_), I46_(ptf.I46_),
    C33_(ptf.C33_), C44_(ptf.C44_), C55_(ptf.C55_), C35_(ptf.C35_),
    rollDampingFraction_(ptf.rollDampingFraction_),
    mooringK11_(ptf.mooringK11_), mooringK22_(ptf.mooringK22_),
    mooringK66_(ptf.mooringK66_),
    mooringDampingFraction_(ptf.mooringDampingFraction_),
    heading_(ptf.heading_),
    addedMass_(ptf.addedMass_),
    addedMassFromTimeStep_(ptf.addedMassFromTimeStep_),
    accelRelaxation_(ptf.accelRelaxation_),
    accelDamping_(ptf.accelDamping_),
    gamma_(ptf.gamma_), beta_(ptf.beta_),
    waveAmp_(ptf.waveAmp_), waveNumber_(ptf.waveNumber_),
    omega_(ptf.omega_), omegaE_(ptf.omegaE_),
    depth_(ptf.depth_), rampTime_(ptf.rampTime_),
    Uinf_(ptf.Uinf_),
    rho_(ptf.rho_),
    steadyRestoring_(ptf.steadyRestoring_),
    Ks_(ptf.Ks_), KsValid_(ptf.KsValid_),
    q_(ptf.q_), qDot_(ptf.qDot_), qDdot_(ptf.qDdot_),
    qNew_(ptf.qNew_), qDotNew_(ptf.qDotNew_), qDdotNew_(ptf.qDdotNew_),
    qDdotPrevIter_(ptf.qDdotPrevIter_),
    lastTimeIndex_(ptf.lastTimeIndex_)
{}


Foam::linBodyMotionFvPatchScalarField::
linBodyMotionFvPatchScalarField
(
    const linBodyMotionFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchField<scalar>(ptf),
    PhiIName_(ptf.PhiIName_),
    UsName_(ptf.UsName_),
    gradUsName_(ptf.gradUsName_),
    pSName_(ptf.pSName_),
    pName_(ptf.pName_),
    fixBody_(ptf.fixBody_),
    mass_(ptf.mass_), xG_(ptf.xG_), zG_(ptf.zG_),
    I44_(ptf.I44_), I55_(ptf.I55_), I66_(ptf.I66_), I46_(ptf.I46_),
    C33_(ptf.C33_), C44_(ptf.C44_), C55_(ptf.C55_), C35_(ptf.C35_),
    rollDampingFraction_(ptf.rollDampingFraction_),
    mooringK11_(ptf.mooringK11_), mooringK22_(ptf.mooringK22_),
    mooringK66_(ptf.mooringK66_),
    mooringDampingFraction_(ptf.mooringDampingFraction_),
    heading_(ptf.heading_),
    addedMass_(ptf.addedMass_),
    addedMassFromTimeStep_(ptf.addedMassFromTimeStep_),
    accelRelaxation_(ptf.accelRelaxation_),
    accelDamping_(ptf.accelDamping_),
    gamma_(ptf.gamma_), beta_(ptf.beta_),
    waveAmp_(ptf.waveAmp_), waveNumber_(ptf.waveNumber_),
    omega_(ptf.omega_), omegaE_(ptf.omegaE_),
    depth_(ptf.depth_), rampTime_(ptf.rampTime_),
    Uinf_(ptf.Uinf_),
    rho_(ptf.rho_),
    steadyRestoring_(ptf.steadyRestoring_),
    Ks_(ptf.Ks_), KsValid_(ptf.KsValid_),
    q_(ptf.q_), qDot_(ptf.qDot_), qDdot_(ptf.qDdot_),
    qNew_(ptf.qNew_), qDotNew_(ptf.qDotNew_), qDdotNew_(ptf.qDdotNew_),
    qDdotPrevIter_(ptf.qDdotPrevIter_),
    lastTimeIndex_(ptf.lastTimeIndex_)
{}


Foam::linBodyMotionFvPatchScalarField::
linBodyMotionFvPatchScalarField
(
    const linBodyMotionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(ptf, iF),
    PhiIName_(ptf.PhiIName_),
    UsName_(ptf.UsName_),
    gradUsName_(ptf.gradUsName_),
    pSName_(ptf.pSName_),
    pName_(ptf.pName_),
    fixBody_(ptf.fixBody_),
    mass_(ptf.mass_), xG_(ptf.xG_), zG_(ptf.zG_),
    I44_(ptf.I44_), I55_(ptf.I55_), I66_(ptf.I66_), I46_(ptf.I46_),
    C33_(ptf.C33_), C44_(ptf.C44_), C55_(ptf.C55_), C35_(ptf.C35_),
    rollDampingFraction_(ptf.rollDampingFraction_),
    mooringK11_(ptf.mooringK11_), mooringK22_(ptf.mooringK22_),
    mooringK66_(ptf.mooringK66_),
    mooringDampingFraction_(ptf.mooringDampingFraction_),
    heading_(ptf.heading_),
    addedMass_(ptf.addedMass_),
    addedMassFromTimeStep_(ptf.addedMassFromTimeStep_),
    accelRelaxation_(ptf.accelRelaxation_),
    accelDamping_(ptf.accelDamping_),
    gamma_(ptf.gamma_), beta_(ptf.beta_),
    waveAmp_(ptf.waveAmp_), waveNumber_(ptf.waveNumber_),
    omega_(ptf.omega_), omegaE_(ptf.omegaE_),
    depth_(ptf.depth_), rampTime_(ptf.rampTime_),
    Uinf_(ptf.Uinf_),
    rho_(ptf.rho_),
    steadyRestoring_(ptf.steadyRestoring_),
    Ks_(ptf.Ks_), KsValid_(ptf.KsValid_),
    q_(ptf.q_), qDot_(ptf.qDot_), qDdot_(ptf.qDdot_),
    qNew_(ptf.qNew_), qDotNew_(ptf.qDotNew_), qDdotNew_(ptf.qDdotNew_),
    qDdotPrevIter_(ptf.qDdotPrevIter_),
    lastTimeIndex_(ptf.lastTimeIndex_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linBodyMotionFvPatchScalarField::readProperties()
{
    const IOdictionary bodyDict
    (
        IOobject
        (
            "bodyMotionProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fixBody_ = bodyDict.getOrDefault<Switch>("fixBody", Switch(false));

    mass_ = bodyDict.get<scalar>("mass");
    xG_   = bodyDict.getOrDefault<scalar>("xG", 0);
    zG_   = bodyDict.getOrDefault<scalar>("zG", 0);

    I44_ = bodyDict.getOrDefault<scalar>("I44", 1);
    I55_ = bodyDict.getOrDefault<scalar>("I55", 1);
    I66_ = bodyDict.getOrDefault<scalar>("I66", 1);
    I46_ = bodyDict.getOrDefault<scalar>("I46", 0);

    C33_ = bodyDict.getOrDefault<scalar>("C33", 0);
    C44_ = bodyDict.getOrDefault<scalar>("C44", 0);
    C55_ = bodyDict.getOrDefault<scalar>("C55", 0);
    C35_ = bodyDict.getOrDefault<scalar>("C35", 0);

    rollDampingFraction_ =
        bodyDict.getOrDefault<scalar>("rollDampingFraction", 0.02);

    mooringK11_ = bodyDict.getOrDefault<scalar>("mooringK11", 0);
    mooringK22_ = bodyDict.getOrDefault<scalar>("mooringK22", 0);
    mooringK66_ = bodyDict.getOrDefault<scalar>("mooringK66", 0);
    mooringDampingFraction_ =
        bodyDict.getOrDefault<scalar>("mooringDampingFraction", 0);

    heading_ = bodyDict.getOrDefault<scalar>("heading", 0);
    rho_     = bodyDict.getOrDefault<scalar>("rhoInf", 1000);

    // Per-DOF artificial added mass; a single scalar is accepted and applied
    // to all six.  Defaults reproduce the values that damped the coupling in
    // practice: light on surge, full inertia on the rest.
    addedMass_ = scalar(0);
    if (bodyDict.found("addedMass"))
    {
        ITstream& is = bodyDict.lookup("addedMass");
        if (is.peek().isNumber())
        {
            addedMass_ = bodyDict.get<scalar>("addedMass");
        }
        else
        {
            const scalarList am(bodyDict.get<scalarList>("addedMass"));
            if (am.size() != 6)
            {
                FatalIOErrorInFunction(bodyDict)
                    << "addedMass needs 6 entries (surge sway heave roll pitch"
                       " yaw) or a single scalar, got " << am.size()
                    << exit(FatalIOError);
            }
            forAll(addedMass_, i) addedMass_[i] = am[i];
        }
    }

    // Forward-speed restoring carried by the steady dynamic pressure: the
    // transfer of p_S onto the displaced hull, the rotation of the normal in
    // it, and the movement of the lever arm.  This is the force-side partner
    // of the m-terms; without it the formulation is not consistent at forward
    // speed.  It is identically zero for a uniform basis flow, so switching
    // it off only matters for a genuine double-body W.
    steadyRestoring_ =
        bodyDict.getOrDefault<Switch>("steadyRestoring", Switch(true));

    addedMassFromTimeStep_ =
        bodyDict.getOrDefault<Switch>("addedMassFromTimeStep", Switch(true));
    accelRelaxation_ = bodyDict.getOrDefault<scalar>("accelRelaxation", 1);
    accelDamping_    = bodyDict.getOrDefault<scalar>("accelDamping", 1);
    gamma_ = bodyDict.getOrDefault<scalar>("newmarkGamma", 0.5);
    beta_  = bodyDict.getOrDefault<scalar>("newmarkBeta", 0.25);

    const IOdictionary waveDict
    (
        IOobject
        (
            "waveConditions",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar g = mag(meshObjects::gravity::New(db().time()).value());

    const scalar steepness = waveDict.get<scalar>("steepness");
    const scalar waveLength = waveDict.get<scalar>("waveLength");
    const scalar U0 = waveDict.get<scalar>("currentSpeed");
    const scalar headAngle = waveDict.get<scalar>("headingAngle");
    const scalar rampPeriods = waveDict.get<scalar>("rampPeriods");
    depth_ = waveDict.get<scalar>("waterDepth");

    waveAmp_    = 0.5*steepness*waveLength;
    waveNumber_ = constant::mathematical::twoPi/waveLength;
    omega_      = Foam::sqrt(g*waveNumber_*Foam::tanh(waveNumber_*depth_));
    Uinf_       = U0*vector(Foam::cos(headAngle), Foam::sin(headAngle), 0);
    omegaE_     = omega_ + waveNumber_*Uinf_.x();
    rampTime_   = rampPeriods*constant::mathematical::twoPi/omega_;
}


Foam::scalar Foam::linBodyMotionFvPatchScalarField::ramp(const scalar t) const
{
    if (rampTime_ <= SMALL || t >= rampTime_)
    {
        return 1.0;
    }

    return 0.5*(1.0 - Foam::cos(constant::mathematical::pi*t/rampTime_));
}


Foam::tmp<Foam::vectorField>
Foam::linBodyMotionFvPatchScalarField::incidentVelocity
(
    const vectorField& pts,
    const scalar t
) const
{
    // u_I = grad(phi_I),  phi_I = A g/w0 cosh k(h+z)/cosh kh sin(kx - we t)
    auto tuI = tmp<vectorField>::New(pts.size(), Zero);
    vectorField& uI = tuI.ref();

    const scalar g = mag(meshObjects::gravity::New(db().time()).value());
    const scalar r = ramp(t);
    const scalar Ua = waveAmp_*waveNumber_*g/omega_*r;
    const scalar chDen = Foam::cosh(waveNumber_*depth_);

    forAll(pts, i)
    {
        const scalar theta = waveNumber_*pts[i].x() - omegaE_*t;
        const scalar cosh_ = Foam::cosh(waveNumber_*(depth_ + pts[i].z()))/chDen;
        const scalar sinh_ = Foam::sinh(waveNumber_*(depth_ + pts[i].z()))/chDen;

        uI[i] = vector
        (
            Ua*cosh_*Foam::cos(theta),
            0,
            Ua*sinh_*Foam::sin(theta)
        );
    }

    return tuI;
}


void Foam::linBodyMotionFvPatchScalarField::bodyLoads
(
    vector& force,
    vector& moment
) const
{
    // Sf points out of the fluid, i.e. into the body, so the force on the body
    // is +sum(p*Sf) -- the same convention that gives buoyancy the right sign.
    const vectorField& Sf = patch().Sf();
    const vectorField& Cf = patch().Cf();

    const auto& p = db().lookupObject<volScalarField>(pName_);
    const scalarField& pb = p.boundaryField()[patch().index()];

    const vector rRef(xG_, 0, 0);

    force = Zero;
    moment = Zero;

    forAll(Sf, i)
    {
        const vector f(rho_*pb[i]*Sf[i]);
        force += f;
        moment += (Cf[i] - rRef) ^ f;
    }

    reduce(force, sumOp<vector>());
    reduce(moment, sumOp<vector>());
}


void Foam::linBodyMotionFvPatchScalarField::steadyRestoringMatrix
(
    tensor6& Ks
) const
{
    Ks = zero66();

    if (!db().foundObject<volScalarField>(pSName_))
    {
        FatalErrorInFunction
            << "steadyRestoring is on but field " << pSName_ << " is not"
               " registered.  Either run a solver that provides the steady"
               " dynamic pressure, or set steadyRestoring false in"
               " constant/bodyMotionProperties."
            << exit(FatalError);
    }

    // Sf points out of the fluid, into the body, exactly as in bodyLoads()
    const vectorField& Sf = patch().Sf();
    const vectorField& Cf = patch().Cf();

    const auto& pS = db().lookupObject<volScalarField>(pSName_);
    const scalarField& pSb = pS.boundaryField()[patch().index()];

    const auto& Us = db().lookupObject<volVectorField>(UsName_);
    const vectorField& W = Us.boundaryField()[patch().index()];

    const auto& gradUs = db().lookupObject<volTensorField>(gradUsName_);
    const tensorField& gradW = gradUs.boundaryField()[patch().index()];

    const vector rRef(xG_, 0, 0);

    // Column k is the load from a unit motion in mode k.  The integrals must
    // be formed in global axes, so the mode goes out through rotate() and the
    // resulting load comes back through it.
    for (int k = 0; k < 6; ++k)
    {
        vector sK(Zero), thK(Zero);
        if (k < 3) { sK[k] = 1; } else { thK[k - 3] = 1; }
        rotate(sK, thK, -heading_);

        vector force(Zero);
        vector moment(Zero);

        forAll(Sf, i)
        {
            const vector r(Cf[i] - rRef);
            const vector S(sK + (thK ^ r));

            // grad(p_S) = -(gradUs & Us), since p_S = -(|W|^2 - |Uinf|^2)/2
            const vector gradPs(-(gradW[i] & W[i]));

            const vector f
            (
                rho_
               *(
                    (S & gradPs)*Sf[i]      // transfer of p_S onto the hull
                  + pSb[i]*(thK ^ Sf[i])    // rotation of the normal in p_S
                )
            );

            force += f;
            moment += (r ^ f) + rho_*pSb[i]*(S ^ Sf[i]);   // lever arm moves
        }

        // Unconditional on every rank, before any master-only work
        reduce(force, sumOp<vector>());
        reduce(moment, sumOp<vector>());

        rotate(force, moment, heading_);

        // The load is F_j = +G_jk eta_k, so on the left-hand side K gains -G
        Ks[0][k] = -force.x();   Ks[1][k] = -force.y();   Ks[2][k] = -force.z();
        Ks[3][k] = -moment.x();  Ks[4][k] = -moment.y();  Ks[5][k] = -moment.z();
    }
}


void Foam::linBodyMotionFvPatchScalarField::assembleSystem
(
    tensor6& M,
    tensor6& C,
    tensor6& K
) const
{
    M = zero66();
    C = zero66();
    K = zero66();

    // Rigid-body inertia about the reference point, including the coupling
    // introduced by a centre of gravity off the waterline
    M[0][0] = mass_;   M[0][4] =  mass_*zG_;
    M[1][1] = mass_;   M[1][3] = -mass_*zG_;
    M[2][2] = mass_;
    M[3][1] = -mass_*zG_;  M[3][3] = I44_;  M[3][5] = -I46_;
    M[4][0] =  mass_*zG_;  M[4][4] = I55_;
    M[5][3] = -I46_;       M[5][5] = I66_;

    // Hydrostatic restoring
    K[2][2] = C33_;
    K[3][3] = C44_;
    K[4][4] = C55_;
    K[2][4] = C35_;
    K[4][2] = C35_;

    // Roll has no wave damping at zero speed in a strictly linear potential
    // model, so a small fraction of critical is added to keep it bounded
    if (C44_ > SMALL)
    {
        C[3][3] = rollDampingFraction_*2.0*Foam::sqrt(I44_*C44_);
    }

    // Soft mooring on surge, sway and yaw, which have no hydrostatic
    // restoring.  Chosen so the planar natural period sits well above the
    // highest wave period, holding the mean position without touching the
    // wave-frequency response.  Zero leaves the DOF free.
    K[0][0] += mooringK11_;
    K[1][1] += mooringK22_;
    K[5][5] += mooringK66_;

    // The other half of the first-order pressure transfer.  C33..C55 above
    // are its hydrostatic part; this is the part carried by the steady
    // DYNAMIC pressure, and the two add.  It contains the Munk moment, so it
    // has negative diagonal entries in pitch and yaw -- which is exactly why
    // it belongs here on the left rather than as a lagged force on the right.
    if (steadyRestoring_)
    {
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                K[i][j] += Ks_[i][j];
            }
        }
    }

    // Critical damping formed from the DRY inertia: the fluid supplies its own
    // added mass through the loads, so this is a numerical bleed on the slow
    // drift oscillation rather than a physical damping model.
    if (mooringDampingFraction_ > SMALL)
    {
        const FixedList<label, 3> moored({0, 1, 5});
        const FixedList<scalar, 3> inertia({mass_, mass_, I66_});
        const FixedList<scalar, 3> stiff
            ({mooringK11_, mooringK22_, mooringK66_});

        forAll(moored, i)
        {
            if (stiff[i] > SMALL)
            {
                C[moored[i]][moored[i]] +=
                    mooringDampingFraction_*2.0
                   *Foam::sqrt(inertia[i]*stiff[i]);
            }
        }
    }
}


void Foam::linBodyMotionFvPatchScalarField::rotate
(
    vector& lin,
    vector& ang,
    const scalar angle
) const
{
    const scalar c = Foam::cos(angle);
    const scalar s = Foam::sin(angle);

    const vector l(lin), a(ang);
    lin = vector(l.x()*c + l.y()*s, -l.x()*s + l.y()*c, l.z());
    ang = vector(a.x()*c + a.y()*s, -a.x()*s + a.y()*c, a.z());
}


void Foam::linBodyMotionFvPatchScalarField::solveMotion()
{
    const scalar dt = db().time().deltaTValue();

    // W, p_S and the mean hull are all steady, so the steady-flow restoring
    // is a constant matrix.  Build it on the first step rather than in
    // readProperties(), which runs at construction, before the solver has
    // solved the basis flow and while Us and pS are still zero.
    if (steadyRestoring_ && !KsValid_)
    {
        steadyRestoringMatrix(Ks_);
        KsValid_ = true;

        if (Pstream::master())
        {
            Info<< "    steady-flow restoring, diag =";
            for (int i = 0; i < 6; ++i) Info<< ' ' << Ks_[i][i];
            Info<< endl;
        }
    }

    tensor6 M, C, K;
    assembleSystem(M, C, K);

    // Artificial added mass, per DOF as a fraction of the diagonal inertia.
    // The body condition sets dphi_D/dn = V_b.n, so the pressure -dphi/dt
    // responds to the body acceleration: the coupling is an added-mass one,
    // evaluated lagged, and goes unstable once the fluid added mass exceeds
    // the body mass.  Carrying A on the left against A*a_old on the right
    // suppresses that.
    tensor6 A = zero66();
    {
        const FixedList<scalar, 6> inertia
            ({mass_, mass_, mass_, I44_, I55_, I66_});
        forAll(A, i) A[i][i] = addedMass_[i]*inertia[i];
    }

    // Loads, rotated from global into body axes
    vector force, moment;
    bodyLoads(force, moment);
    rotate(force, moment, heading_);

    vector6 f = zero6();
    f[0] = force.x();  f[1] = force.y();  f[2] = force.z();
    f[3] = moment.x(); f[4] = moment.y(); f[5] = moment.z();

    // Newmark predictors
    const vector6 qStar =
        q_ + dt*qDot_ + (dt*dt*(0.5 - beta_))*qDdot_;
    const vector6 qDotStar =
        qDot_ + (dt*(1.0 - gamma_))*qDdot_;

    tensor6 lhs = zero66();
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            lhs[i][j] = M[i][j] + A[i][j] + gamma_*dt*C[i][j] + beta_*dt*dt*K[i][j];

    // See addedMassFromTimeStep_: the previous time step stabilises, the
    // previous corrector is consistent.
    const vector6& aDrive =
        addedMassFromTimeStep_ ? qDdot_ : qDdotPrevIter_;

    const vector6 rhs =
        (f + mult6(A, aDrive)) - mult6(C, qDotStar) - mult6(K, qStar);

    qDdotNew_ = solve6(lhs, rhs);

    // Optional under-relaxation / decay on the acceleration
    if (accelRelaxation_ < 1 - SMALL || accelDamping_ < 1 - SMALL)
    {
        forAll(qDdotNew_, i)
        {
            qDdotNew_[i] = accelDamping_
               *(
                    accelRelaxation_*qDdotNew_[i]
                  + (1 - accelRelaxation_)*qDdotPrevIter_[i]
                );
        }
    }

    if (fixBody_)
    {
        // Restrained: zero acceleration keeps the predictors at zero for all
        // time, so the condition reduces to the pure diffraction problem
        qDdotNew_ = zero6();
    }

    qDotNew_ = qDotStar + (gamma_*dt)*qDdotNew_;
    qNew_    = qStar + (beta_*dt*dt)*qDdotNew_;

    qDdotPrevIter_ = qDdotNew_;
}


void Foam::linBodyMotionFvPatchScalarField::writeMotion()
{
    if (!Pstream::master()) return;

    if (!motionFilePtr_)
    {
        const fileName dir(db().time().globalPath()/"postProcessing"/"bodyMotion");
        mkDir(dir);
        motionFilePtr_.reset(new OFstream(dir/"motion.dat"));
        *motionFilePtr_
            << "# Linearised rigid-body motion, patch " << patch().name() << nl
            << "# Time\tsurge\tsway\theave\troll\tpitch\tyaw"
               "\tsurgeDot\tswayDot\theaveDot\trollDot\tpitchDot\tyawDot"
            << endl;
    }

    *motionFilePtr_ << db().time().timeOutputValue();
    for (int i = 0; i < 6; ++i) *motionFilePtr_ << '\t' << q_[i];
    for (int i = 0; i < 6; ++i) *motionFilePtr_ << '\t' << qDot_[i];
    *motionFilePtr_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linBodyMotionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().value();

    // Shift the state once per time step; the non-orthogonal correctors then
    // re-solve the motion against an updated pressure without advancing time.
    if (db().time().timeIndex() != lastTimeIndex_)
    {
        q_ = qNew_;
        qDot_ = qDotNew_;
        qDdot_ = qDdotNew_;
        qDdotPrevIter_ = qDdot_;

        writeMotion();
        lastTimeIndex_ = db().time().timeIndex();
    }

    solveMotion();

    // --- Assemble the normal gradient ------------------------------------
    const vectorField n(patch().nf());
    const vectorField& Cf = patch().Cf();
    const vector rRef(xG_, 0, 0);

    // Motion back into global axes
    vector Ulin(qDotNew_[0], qDotNew_[1], qDotNew_[2]);
    vector omega(qDotNew_[3], qDotNew_[4], qDotNew_[5]);
    rotate(Ulin, omega, -heading_);

    vector X(qNew_[0], qNew_[1], qNew_[2]);
    vector theta(qNew_[3], qNew_[4], qNew_[5]);
    rotate(X, theta, -heading_);

    const tmp<vectorField> tuI = incidentVelocity(Cf, t);
    const vectorField& uI = tuI();

    const auto& Us = db().lookupObject<volVectorField>(UsName_);
    const vectorField& W = Us.boundaryField()[patch().index()];

    const auto& gradUs = db().lookupObject<volTensorField>(gradUsName_);
    const tensorField& gradW = gradUs.boundaryField()[patch().index()];

    scalarField& g = gradient();
    g.setSize(patch().size());

    forAll(g, i)
    {
        const vector r(Cf[i] - rRef);

        // Body velocity and displacement at this face
        const vector Vb(Ulin + (omega ^ r));
        const vector S(X + (theta ^ r));

        // m-terms: (n.grad)W is grad(W)^T & n with OpenFOAM's grad(U)_ij = d_i U_j
        const vector nGradW(gradW[i].T() & n[i]);

        g[i] =
            (Vb & n[i])                 // body motion
          - (uI[i] & n[i])              // incident flux to be cancelled
          - (S & nGradW)                // m1..m3 and part of m4..m6
          - (W[i] & (theta ^ n[i]));    // rotation of the normal in W
    }

    if (Pstream::master())
    {
        Info<< "    body: heave " << qNew_[2]
            << "  pitch " << qNew_[4] << endl;
    }

    fixedGradientFvPatchField<scalar>::updateCoeffs();
}


void Foam::linBodyMotionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("PhiI", "PhiI", PhiIName_);
    os.writeEntryIfDifferent<word>("Us", "Us", UsName_);
    os.writeEntryIfDifferent<word>("gradUs", "gradUs", gradUsName_);
    os.writeEntryIfDifferent<word>("pS", "pS", pSName_);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        linBodyMotionFvPatchScalarField
    );
}

// ************************************************************************* //
