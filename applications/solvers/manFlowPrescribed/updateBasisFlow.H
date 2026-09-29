// Basis flow for the frame motion at time tBasis (quasi-steady).
//
// PhiS = V0x Phi1 + V0y Phi2 + Omega_z(t) Phi6, the frame velocity V_S, the
// relative flow Us = W = grad(PhiS) - V_S, gradUs, pS = -1/2(|W|^2 - |V_S|^2)
// and dUsdz.  Called once at startup and then on every step on which the yaw
// rate changes (the yaw-rate onset ramp).  basisFlowIndex is incremented so
// the boundary conditions rebuild what they cache from the basis flow (the
// m-terms, the steady-flow restoring and the free-surface upwind schemes).
//
// With a time-dependent basis flow, dPhiS/dt|ship = dOmega/dt Phi6 is left
// out of the pressure: the basis flow is treated as quasi-steady during the
// ramp, which lasts several encounter periods.  The loads during the ramp
// are therefore transitional and are not used for analysis.
{
    const scalar Oznow = motion.Omega(tBasis).z();

    // Phi1, Phi2 have dimensions of length (dPhi/dn is a direction cosine),
    // Phi6 of length^2 (dPhi/dn is a lever arm)
    const dimensionedScalar V0x("V0x", dimVelocity, motion.V0().x());
    const dimensionedScalar V0y("V0y", dimVelocity, motion.V0().y());
    const dimensionedScalar Oz("Oz", dimless/dimTime, Oznow);

    PhiS == V0x*Phi1 + V0y*Phi2 + Oz*Phi6;

    // Frame velocity on cells and faces
    {
        vectorField& VSc = VS.primitiveFieldRef();
        const vectorField& C = mesh.C();
        forAll(C, i) VSc[i] = motion.VS(C[i], tBasis);

        forAll(mesh.boundary(), patchi)
        {
            vectorField& VSb = VS.boundaryFieldRef()[patchi];
            const vectorField& Cf = mesh.boundary()[patchi].Cf();
            forAll(Cf, i) VSb[i] = motion.VS(Cf[i], tBasis);
        }
    }

    // W = grad(PhiS) - V_S, cells and faces.  On the hull the normal component
    // of grad(PhiS) is V_S.n by construction, so W.n = 0 there.
    Us == fvc::grad(PhiS) - VS;

    gradUs = fvc::grad(Us);

    // Steady dynamic pressure, pS = -1/2 (|W|^2 - |V_S|^2)
    pS == -0.5*(magSqr(Us) - magSqr(VS));

    // dWz/dz on the free surface from continuity; the rigid rotation has no
    // divergence and no diagonal gradient, so this is still exact.
    dUsdz = -gradUs.component(tensor::XX) - gradUs.component(tensor::YY);

    OzBasis = Oznow;
    basisFlowIndex.value() += 1;
    frameOmega.value() = vector(0, 0, Oznow);
}
