// Prescribed motion and incident wave (constant/prescribedMotion,
// constant/waveConditions), shared with the boundary conditions
const prescribedShipMotion motion(mesh);

const vector xc(motion.rotationCentre());

IOdictionary motionDict
(
    IOobject
    (
        "prescribedMotion",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const word hullPatchName
(
    motionDict.getOrDefault<word>("hullPatch", "sphere")
);

// Reference density, only for the added-mass printout
const scalar rhoRef
(
    IOdictionary
    (
        IOobject
        (
            "bodyMotionProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ).getOrDefault<scalar>("rhoInf", 1000)
);

Info<< "Prescribed motion (MMG axes at the rotation centre " << xc << ")" << nl
    << "    u = " << motion.u() << " m/s   v = " << motion.v()
    << " m/s   r = " << motion.r() << " rad/s" << nl
    << "    frame velocity V0 = " << motion.V0()
    << "   Omega = " << motion.Omega() << "  (mesh axes)" << nl
    << "Incident wave" << nl
    << "    amplitude   " << motion.amplitude() << " m" << nl
    << "    wavenumber  " << motion.waveNumber() << " 1/m" << nl
    << "    omega       " << motion.omega() << " rad/s" << nl
    << "    direction at t = 0 (mesh axes) " << motion.waveDir(0)
    << "   encounter angle " << motion.encounterAngle(0) << " deg" << nl
    << endl;


// --- Potentials ------------------------------------------------------------
// Convention throughout: u = +grad(Phi)

// Unit basis potentials: fixedGradient on the hull (set by the solver),
// fixedValue 0 on the outer boundaries, zeroGradient on upper and lower
Info<< "Reading unit basis potentials Phi1, Phi2, Phi6\n" << endl;
volScalarField Phi1
(
    IOobject("Phi1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
    mesh
);
volScalarField Phi2
(
    IOobject("Phi2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
    mesh
);
volScalarField Phi6
(
    IOobject("Phi6", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
    mesh
);

// Absolute double-body potential of the prescribed motion
volScalarField PhiS
(
    IOobject("PhiS", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero)
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

// Frame velocity V_S(x): the earth velocity of the ship-fixed point x
volVectorField VS
(
    IOobject("VS", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

// Steady flow relative to the ship, W = grad(PhiS) - V_S
Info<< "Creating steady relative basis velocity Us\n" << endl;
volVectorField Us
(
    IOobject("Us", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
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
// Kinematic, p/rho.  Unsteady:  p = -(dPhi/dt|ship + W.grad(Phi))
//                    Steady  :  pS = -1/2 (|W|^2 - |V_S|^2)

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero)
);

volScalarField pS
(
    IOobject("pS", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero)
);

// Analytic dPhiI/dt at a fixed mesh point (local encounter frequency)
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
    IOobject("gradUs", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedTensor(dimless/dimTime, Zero)
);

// dWz/dz on the free surface.  Wz is odd in z for a double-body flow, so this
// does not vanish there, unlike the tangential derivatives.
volScalarField dUsdz
(
    IOobject("dUsdz", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
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
