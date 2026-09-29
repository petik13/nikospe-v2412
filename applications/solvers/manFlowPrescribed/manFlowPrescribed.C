/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Application
    manFlowPrescribed

Description
    manFlow for a ship on a PRESCRIBED planar motion: a circular motion test
    (constant u, v, r) in waves.  Time-domain linear seakeeping on a
    double-body basis flow, in the ship-fixed (rotating) frame.

    Velocity convention is u = +grad(Phi) throughout.

    The mesh is ship-fixed (bow at -x) and never moves.  The frame velocity
    V_S(x) = (-u, v, 0) + Omega x (x - x_c) and the incident wave are given by
    prescribedShipMotion (constant/prescribedMotion, constant/waveConditions).

    Basis flow
    ----------
    PhiS is now the ABSOLUTE double-body potential (fluid at rest far away),
    split Kirchhoff-wise into unit potentials solved once at startup,

        PhiS = V0x Phi1 + V0y Phi2 + Omega_z Phi6,
        dPhi1/dn = n_x,  dPhi2/dn = n_y,  dPhi6/dn = (x - x_c) n_y - (y - y_c) n_x

    on the hull, Phi = 0 on the outer boundaries and zero normal gradient on
    z = 0 (double body) and the bottom.  The flow relative to the ship is

        W = Us = grad(PhiS) - V_S(x)

    which is what every term of manFlow calls Us.  For a pure translation it
    is identical to manFlow's grad(PhiS), with PhiS -> Uinf.x far away.  With
    rotation W is not a potential flow (curl W = -2 Omega); the free-surface
    and body conditions are written with it directly, following Yang & Kim
    (2026, eqs 4-6).

    Steady pressure (Bernoulli for the absolute potential in the ship frame):

        pS = -1/2 (|W|^2 - |V_S|^2)

    and the unsteady kinematic pressure is unchanged,

        p = -(dPhi/dt|ship + W.grad(Phi))

    with dPhiI/dt|ship analytic (local encounter frequency).

    Yaw-rate onset
    --------------
    The ship runs straight until yawOnsetTime (so the motions settle) and the
    yaw rate is then ramped in over yawRampTime (constant/prescribedMotion).
    During the ramp the basis flow is rebuilt every step (updateBasisFlow.H,
    quasi-steady: dPhiS/dt is not in the pressure) and basisFlowIndex is
    incremented, so linBodyMotionRot rebuilds its m-terms and steady-flow
    restoring and potRotatingFrameBC its upwind schemes.  frameOmega carries
    the current Omega for middleFieldFormRot.

    Output
    ------
    postProcessing/shipTrajectory/trajectory.dat: time, heading, earth
    position, encounter angle (manModel convention) and u, v, r(t).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "mathematicalConstants.H"
#include "fixedGradientFvPatchFields.H"
#include "OFstream.H"
#include "uniformDimensionedFields.H"
#include "prescribedShipMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Time-domain linear seakeeping for a ship on a prescribed planar motion"
        " (circular motion test in waves)"
    );

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

    // --- Unit basis potentials --------------------------------------------
    {
        const label hullID = mesh.boundaryMesh().findPatchID(hullPatchName);
        if (hullID < 0)
        {
            FatalErrorInFunction
                << "Hull patch " << hullPatchName << " not found"
                << exit(FatalError);
        }

        const vectorField nf(mesh.boundary()[hullID].nf());
        const vectorField d(mesh.boundary()[hullID].Cf() - xc);

        // Patch normal points out of the fluid, into the body.  The body
        // condition grad(Phi).n = V_S.n is the same for either orientation.
        const scalarField g1(nf.component(vector::X));
        const scalarField g2(nf.component(vector::Y));
        const scalarField g6
        (
            d.component(vector::X)*nf.component(vector::Y)
          - d.component(vector::Y)*nf.component(vector::X)
        );

        const FixedList<const scalarField*, 3> gradients({&g1, &g2, &g6});
        const FixedList<volScalarField*, 3> units({&Phi1, &Phi2, &Phi6});

        forAll(units, k)
        {
            volScalarField& Phik = *units[k];

            if (!isA<fixedGradientFvPatchScalarField>(Phik.boundaryField()[hullID]))
            {
                FatalErrorInFunction
                    << Phik.name() << " must be fixedGradient on "
                    << hullPatchName << exit(FatalError);
            }

            refCast<fixedGradientFvPatchScalarField>
            (
                Phik.boundaryFieldRef()[hullID]
            ).gradient() = *gradients[k];

            Info<< "Solving unit basis potential " << Phik.name() << endl;

            while (basisFlow.correctNonOrthogonal())
            {
                fvScalarMatrix PhikEqn
                (
                    fvm::laplacian(dimensionedScalar("1", dimless, 1), Phik)
                 == Zero
                );
                PhikEqn.solve();
            }

            Phik.write();
        }
    }

    // Rigid-lid (double-body) added masses of the underwater hull from the
    // unit potentials, a check against addedMassDict (Motora) and the zero-
    // frequency limit of WAMIT:
    //
    //     m_jk = -rho int_H Phi_j dPhi_k/dn_b dS,   n_b into the fluid
    //
    // The patch normal is n_p = -n_b and snGrad is along n_p, so
    // m_jk = +rho int_H Phi_j snGrad(Phi_k) dS.  (Sphere check: Phi_1 =
    // -a cos(th)/2 on r = a gives m_11 = 2/3 pi rho a^3.)  Mesh axes: x aft,
    // y starboard; m66 is about the rotation centre.
    {
        const label hullID = mesh.boundaryMesh().findPatchID(hullPatchName);
        const scalarField& mA = mesh.boundary()[hullID].magSf();
        const scalarField p1(Phi1.boundaryField()[hullID]);
        const scalarField p2(Phi2.boundaryField()[hullID]);
        const scalarField p6(Phi6.boundaryField()[hullID]);
        const scalarField g1(Phi1.boundaryField()[hullID].snGrad());
        const scalarField g2(Phi2.boundaryField()[hullID].snGrad());
        const scalarField g6(Phi6.boundaryField()[hullID].snGrad());

        const scalar m11 = rhoRef*gSum(p1*g1*mA);
        const scalar m22 = rhoRef*gSum(p2*g2*mA);
        const scalar m66 = rhoRef*gSum(p6*g6*mA);
        const scalar m26 = rhoRef*gSum(p2*g6*mA);

        Info<< nl << "Rigid-lid added masses of the underwater hull (rho = "
            << rhoRef << "):" << nl
            << "    m11 = " << m11 << "   m22 = " << m22
            << "   m66 = " << m66 << "   m26 = " << m26 << nl << endl;
    }

    // --- Basis flow for the prescribed motion at the start time -------------
    {
        const scalar tBasis = runTime.value();
        #include "updateBasisFlow.H"
    }

    PhiS.write();
    Us.write();
    gradUs.write();
    pS.write();
    dUsdz.write();
    VS.write();

    Info<< "    max |Us| = " << max(mag(Us)).value()
        << "   max |V_S| = " << max(mag(VS)).value()
        << "   max |dUsdz| = " << max(mag(dUsdz)).value() << nl << endl;


    // --- Trajectory output --------------------------------------------------
    autoPtr<OFstream> trajFile;
    if (Pstream::master())
    {
        const fileName dir(runTime.globalPath()/"postProcessing"/"shipTrajectory");
        mkDir(dir);
        trajFile.reset(new OFstream(dir/"trajectory.dat"));
        trajFile()
            << "# Prescribed motion, MMG axes at the rotation centre" << nl
            << "# u = " << motion.u() << "  v = " << motion.v()
            << "  r = " << motion.rTarget()
            << "  (yaw onset at t = " << motion.yawOnsetTime()
            << " s, ramp " << motion.yawRampTime() << " s)" << nl
            << "# Time\tpsi[deg]\tX[m]\tY[m]\tmu[deg]\tu\tv\tr" << endl;
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Starting time loop" << nl << endl;

    runTime.functionObjects().start();

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName()
            << "  deltaT = " << runTime.deltaTValue() << nl << endl;

        const scalar t = runTime.value();

        // --- Basis flow: rebuild while the yaw rate changes (onset ramp) ---
        if (mag(motion.Omega(t).z() - OzBasis) > 1e-12*max(mag(motion.rTarget()), 1e-12))
        {
            const scalar tBasis = t;
            #include "updateBasisFlow.H"

            Info<< "    basis flow updated: r = " << motion.r(t)
                << " rad/s (index " << basisFlowIndex.value() << ")" << endl;
        }

        // --- Incident wave ------------------------------------------------
        {
            scalarField& PhiIc = PhiI.primitiveFieldRef();
            vectorField& UIc = UI.primitiveFieldRef();
            scalarField& dPhiIdtc = dPhiIdt.primitiveFieldRef();
            scalarField& zetaIc = zetaI.primitiveFieldRef();

            vector gz;
            const vectorField& C = mesh.C();
            forAll(C, celli)
            {
                motion.incident
                (
                    C[celli], t, PhiIc[celli], UIc[celli],
                    dPhiIdtc[celli], zetaIc[celli], gz
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
                    motion.incident
                    (
                        Cf[facei], t, PhiIb[facei], UIb[facei],
                        dPhiIdtb[facei], zetaIb[facei], gz
                    );
                }
            }
        }

        // --- Disturbance potential ----------------------------------------
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

            // Linearised Bernoulli in the ship frame:  p = -(dPhi/dt + W.grad(Phi))
            p = -(dPhiIdt + fvc::ddt(PhiD) + (Us & U));
        }

        Info<< "    continuity error = "
            << mag(fvc::div(U))().weightedAverage(mesh.V()).value() << endl;

        if (Pstream::master())
        {
            scalar X, Y;
            motion.position(t, X, Y);
            trajFile()
                << t << tab << radToDeg(motion.psi(t)) << tab << X << tab << Y
                << tab << motion.encounterAngle(t)
                << tab << motion.u() << tab << motion.v() << tab << motion.r(t)
                << endl;
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    runTime.functionObjects().end();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
