// Wave and current conditions, shared with the boundary conditions
IOdictionary waveConditions
(
    IOobject
    (
        "waveConditions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const scalar gMag = mag(g.value());

const scalar steepness   = waveConditions.get<scalar>("steepness");    // wave steepness, H/L
const scalar waveLength  = waveConditions.get<scalar>("waveLength");   // Wave length, m
const scalar waterDepth  = waveConditions.get<scalar>("waterDepth");   // Water depth, m
const scalar currentSpeed = waveConditions.get<scalar>("currentSpeed"); // Forward current speed, m/s
const scalar headingAngle = waveConditions.get<scalar>("headingAngle"); // Heading angle of the current
const scalar rampPeriods  = waveConditions.get<scalar>("rampPeriods");  // Ramp time in wave periods
const scalar swaySpeed    = waveConditions.getOrDefault<scalar>("swaySpeed", 0); // Sway current speed, m/s

// Basis-flow linearisation.
//   DBL  double body: solve the steady BVP, W is the ship-shaped stream.
//        W.n = 0 on the hull and Wz = 0 on z = 0, and the difference

//   NKL  Neumann-Kelvin: no steady solve at all, W = Uinf everywhere.  Then
//        grad(W) = 0 and pS = 0, so m1..m3, the steady restoring and every
//        (W - Uinf) forcing vanish identically; only the -n x W half of
//        m4..m6 survives, which is the classical uniform-stream result.
const word basisLinearisation
(
    waveConditions.getOrDefault<word>("linearisation", "DBL")
);

if (basisLinearisation != "DBL" && basisLinearisation != "NKL")
{
    FatalIOErrorInFunction(waveConditions)
        << "linearisation must be DBL or NKL, got " << basisLinearisation
        << exit(FatalIOError);
}

// Derived wave properties.  omega is the absolute frequency from the
// dispersion relation; omegaE is the encounter frequency
const scalar waveAmp = 0.5*steepness*waveLength;
const scalar waveNumber = constant::mathematical::twoPi/waveLength;
const scalar omega = Foam::sqrt(gMag*waveNumber*Foam::tanh(waveNumber*waterDepth)); // absolute wave frequency

// Freestream velocity vector
const vector Uinf
(
    currentSpeed*vector(Foam::cos(headingAngle), Foam::sin(headingAngle), 0) + 
    swaySpeed*vector(Foam::cos(headingAngle + constant::mathematical::pi/2), Foam::sin(headingAngle + constant::mathematical::pi/2), 0)
);

const scalar omegaE = omega + waveNumber*Uinf.x(); // Encounter frequency
const scalar wavePeriod = constant::mathematical::twoPi/omega;
const scalar rampTime = rampPeriods*wavePeriod;

Info<< "Incident wave" << nl
    << "    length      " << waveLength << " m" << nl
    << "    amplitude   " << waveAmp << " m" << nl
    << "    wavenumber  " << waveNumber << " 1/m" << nl
    << "    omega       " << omega << " rad/s   (period " << wavePeriod << " s)" << nl
    << "    omega_e     " << omegaE << " rad/s" << nl
    << "    basis flow  " << basisLinearisation
    << (basisLinearisation == "NKL" ? "  (uniform stream, no steady solve)" : "  (double body)") << nl
    << "    Uinf        " << Uinf << " m/s   (tau = "
    << omegaE*currentSpeed/gMag << ")" << nl
    << "    depth       " << waterDepth << " m" << nl
    << "    ramp        " << rampTime << " s" << nl << nl
    << "    sway speed   " << swaySpeed << " m/s" << endl;


// --- Potentials ------------------------------------------------------------
// Convention throughout: u = +grad(Phi)

Info<< "Reading steady basis potential PhiS\n" << endl;
volScalarField PhiS
(
    IOobject("PhiS", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Creating incident potential PhiI\n" << endl;
volScalarField PhiI
(
    IOobject("PhiI", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero)
);

Info<< "Reading disturbance potential PhiD\n" << endl;
volScalarField PhiD
(
    IOobject("PhiD", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Creating total unsteady potential Phi\n" << endl;
volScalarField Phi
(
    IOobject("Phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero)
);


// --- Velocities ------------------------------------------------------------

Info<< "Creating steady basis velocity Us\n" << endl;
volVectorField Us
(
    IOobject("Us", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

volVectorField UI
(
    IOobject("UI", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

volVectorField UD
(
    IOobject("UD", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);


// --- Free-surface elevation ------------------------------------------------
// Only the free-surface patch values are meaningful; the internal field is
// carried so the boundary condition can use oldTime().

volScalarField zetaI
(
    IOobject("zetaI", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimLength, Zero)
);

volScalarField zetaD
(
    IOobject("zetaD", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimLength, Zero)
);

volScalarField zeta
(
    IOobject("zeta", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimLength, Zero)
);


// --- Pressures -------------------------------------------------------------
// Kinematic, p/rho.  Unsteady:  p = -(dPhi/dt + W.grad(Phi))
//                    Steady  :  pS = -1/2 (|W|^2 - |Uinf|^2)

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero)
);

volScalarField pS
(
    IOobject("pS", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero)
);

// Analytic dPhiI/dt: the incident potential is prescribed, so its time
// derivative is taken in closed form rather than differenced.  This keeps the
// Froude-Krylov pressure exact.
volScalarField dPhiIdt
(
    IOobject("dPhiIdt", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero)
);


// --- Steady basis flow derivatives -----------------------------------------

Info<< "Creating steady basis velocity gradient gradUs\n" << endl;
volTensorField gradUs
(
    IOobject("gradUs", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::grad(Us)
);

// dWz/dz on the free surface.  Wz is odd in z for a double-body flow, so this
// does not vanish there, unlike the tangential derivatives.
volScalarField dUsdz
(
    IOobject("dUsdz", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

// The standard forces function object asks for a viscous stress model.  There
// is none in a potential flow, so a zero-viscosity transportProperties is
// registered to satisfy it; the viscous contribution is then identically zero.
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::REGISTER
    )
);

if (!transportProperties.found("nu"))
{
    transportProperties.add("nu", dimensionedScalar(dimViscosity, Zero));
}


label PhiRefCell = 0;
scalar PhiRefValue = 0;
setRefCell(PhiD, potentialFlow.dict(), PhiRefCell, PhiRefValue);
