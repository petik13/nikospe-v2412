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

Application
    manFlow

Group
    grpBasicSolvers

Description
    Time-domain linear seakeeping solver for a body advancing in waves, on a
    double-body basis flow.

    Velocity convention is u = +grad(Phi) throughout.

    The total potential is split three ways,

        Phi_total = PhiS + PhiI + PhiD

    with PhiS the steady double-body basis flow (solved once, at startup),
    PhiI the incident wave (prescribed analytically at every step) and PhiD the
    disturbance, which is the only unknown.  PhiD satisfies Laplace's equation

        laplacian(PhiD) = 0

    subject to the linearised free-surface condition on the upper patch
    (potForwardSpeedBC) and the body condition on the hull (linBodyMotion).
    Both are explicit in time, so each step is a single elliptic solve.

    Prescribing PhiI rather than generating it numerically keeps the incident
    wave exact everywhere, and lets the sponge layers damp the disturbance alone.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Time-domain linear seakeeping solver for a body advancing in waves"
    );

    // Accepted for compatibility with the run scripts; a solver executes its
    // function objects from runTime.loop() regardless.
    argList::addBoolOption("withFunctionObjects", "Execute functionObjects");

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    pisoControl potentialFlow(mesh, "potentialFlow"); // unsteady flow
    pisoControl basisFlow(mesh, "basisFlow");         // steady flow

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (basisLinearisation == "NKL")
    {
        Info<< nl << "Basis flow: Neumann-Kelvin, uniform stream (no steady"
               " solve)" << endl;

        const dimensionedVector UinfD("Uinf", dimVelocity, Uinf);

        PhiS == (UinfD & mesh.C());
        Us == UinfD;
    }
    else
    {
        Info<< nl << "Solving the steady basis flow PhiS" << endl;

        while (basisFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiSEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), PhiS) == Zero
            );

            PhiSEqn.solve();
        }

        Us = fvc::grad(PhiS);
    }
    gradUs = fvc::grad(Us);

    // Steady dynamic pressure due to the basis flow, pS = -1/2 (|W|^2 - |Uinf|^2).
    pS = -0.5*(magSqr(Us) - dimensionedScalar(sqr(dimVelocity), magSqr(Uinf)));

    // dUsdz = gradUs.component(tensor::ZZ); // dWz/dz on the free surface, for the linearised kinematic FSBC.
    // Take dUsdz from continuity
    dUsdz = -gradUs.component(tensor::XX) - gradUs.component(tensor::YY);

    Us.write();
    gradUs.write();
    pS.write();
    dUsdz.write();

    Info<< "    max |Us| = " << max(mag(Us)).value()
        << "   max |dUsdz| = " << max(mag(dUsdz)).value() << nl << endl;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Starting time loop" << nl << endl;

    runTime.functionObjects().start();

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName()
            << "  deltaT = " << runTime.deltaTValue() << nl << endl;

        const scalar t = runTime.value();

        // --- Incident wave ------------------------------------------------
        // Prescribed in closed form on cells and on every patch, so that the
        // Froude-Krylov part of the loads carries no discretisation error.
        {
            // Ramping of the incident wave to avoid ill-conditioning at startup. - Until rampTime.
            const scalar rampFactor =
                (rampTime > SMALL && t < rampTime)
              ? 0.5*(1.0 - Foam::cos(constant::mathematical::pi*t/rampTime))
              : 1.0;

            // The ramp multiplies the AMPLITUDE, so d(phi_I)/dt picks up a
            // term in rdot as well as the one in the phase. 
            const scalar rampFactorDot =
                (rampTime > SMALL && t < rampTime)
              ? 0.5*constant::mathematical::pi/rampTime
                   *Foam::sin(constant::mathematical::pi*t/rampTime)
              : 0.0;

            const scalar coshDen = Foam::cosh(waveNumber*waterDepth);
            const scalar Pa = waveAmp*gMag/omega*rampFactor;          // phi_I
            const scalar Ua = Pa*waveNumber; // u_I
            const scalar dPa = -Pa*omegaE;  // dphi_I/dt
            const scalar dPaR = waveAmp*gMag/omega*rampFactorDot;      // rdot part

            // phi_I = A g/w0 cosh k(h+z)/cosh kh sin(kx - we t)
            auto setIncident = [&]
            (
                const vector& C,
                scalar& phiI,
                vector& uI,
                scalar& dphiIdt,
                scalar& elevation
            )
            {
                const scalar theta = waveNumber*C.x() - omegaE*t;
                const scalar cth = Foam::cos(theta);
                const scalar sth = Foam::sin(theta);
                const scalar ch = Foam::cosh(waveNumber*(waterDepth + C.z()))/coshDen;
                const scalar sh = Foam::sinh(waveNumber*(waterDepth + C.z()))/coshDen;

                phiI = Pa*ch*sth;
                uI = vector(Ua*ch*cth, 0, Ua*sh*sth);
                dphiIdt = dPa*ch*cth + dPaR*ch*sth;
                elevation = waveAmp*cth*rampFactor;
            };

            scalarField& PhiIc = PhiI.primitiveFieldRef();
            vectorField& UIc = UI.primitiveFieldRef();
            scalarField& dPhiIdtc = dPhiIdt.primitiveFieldRef();
            scalarField& zetaIc = zetaI.primitiveFieldRef();

            const vectorField& C = mesh.C();
            forAll(C, celli)
            {
                setIncident
                (
                    C[celli], PhiIc[celli], UIc[celli],
                    dPhiIdtc[celli], zetaIc[celli]
                );
            }

            forAll(mesh.boundary(), patchi)
            {
                const vectorField& Cf = mesh.boundary()[patchi].Cf();

                scalarField& PhiIb = PhiI.boundaryFieldRef()[patchi];
                vectorField& UIb = UI.boundaryFieldRef()[patchi];
                scalarField& dPhiIdtb = dPhiIdt.boundaryFieldRef()[patchi];
                scalarField& zetaIb = zetaI.boundaryFieldRef()[patchi];

                forAll(Cf, facei)
                {
                    setIncident
                    (
                        Cf[facei], PhiIb[facei], UIb[facei],
                        dPhiIdtb[facei], zetaIb[facei]
                    );
                }
            }
        }

        // --- Disturbance potential ----------------------------------------
        // The free-surface and body conditions are explicit, so this is a
        // single Laplace solve per step; the correctors only chase mesh
        // non-orthogonality and the body-motion coupling.
        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiDEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), PhiD) == Zero
            );

            PhiDEqn.setReference(PhiRefCell, PhiRefValue);
            PhiDEqn.solve();

            UD = fvc::grad(PhiD);

            Phi = PhiI + PhiD;
            U = UI + UD;
            zeta = zetaI + zetaD;

            // Linearised Bernoulli, kinematic:  p = -(dPhi/dt + W.grad(Phi))
            // dPhiI/dt is analytic; only the disturbance part is differenced.
            p = -(dPhiIdt + fvc::ddt(PhiD) + (Us & U));
        }

        Info<< "    continuity error = "
            << mag(fvc::div(U))().weightedAverage(mesh.V()).value() << endl;

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    runTime.functionObjects().end();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
