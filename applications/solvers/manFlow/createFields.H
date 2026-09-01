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

const scalar steepness   = waveConditions.get<scalar>("steepness");
const scalar waveLength  = waveConditions.get<scalar>("waveLength");
const scalar waterDepth  = waveConditions.get<scalar>("waterDepth");
const scalar currentSpeed = waveConditions.get<scalar>("currentSpeed");
const scalar headingAngle = waveConditions.get<scalar>("headingAngle");
const scalar rampPeriods  = waveConditions.get<scalar>("rampPeriods");

// Derived wave properties.  omega is the absolute frequency from the
// dispersion relation; omegaE is what the body sees once the uniform stream
// Doppler-shifts it.
const scalar waveAmp = 0.5*steepness*waveLength;
const scalar waveNumber = constant::mathematical::twoPi/waveLength;
const scalar omega =
    Foam::sqrt(gMag*waveNumber*Foam::tanh(waveNumber*waterDepth));

const vector Uinf
(
    currentSpeed*vector(Foam::cos(headingAngle), Foam::sin(headingAngle), 0)
);

const scalar omegaE = omega + waveNumber*Uinf.x();
const scalar wavePeriod = constant::mathematical::twoPi/omega;
const scalar rampTime = rampPeriods*wavePeriod;

Info<< "Incident wave" << nl
    << "    length      " << waveLength << " m" << nl
    << "    amplitude   " << waveAmp << " m" << nl
    << "    wavenumber  " << waveNumber << " 1/m" << nl
    << "    omega       " << omega << " rad/s   (period " << wavePeriod << " s)" << nl
    << "    omega_e     " << omegaE << " rad/s" << nl
    << "    Uinf        " << Uinf << " m/s   (tau = "
    << omegaE*currentSpeed/gMag << ")" << nl
    << "    depth       " << waterDepth << " m" << nl
    << "    ramp        " << rampTime << " s" << nl << endl;


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
// does not vanish there, unlike the tangential derivatives.  Taken from
// continuity rather than from a one-sided normal difference, because the
// tangential derivatives are the well-resolved ones on a flat patch.
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
