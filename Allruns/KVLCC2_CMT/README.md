# KVLCC2 circular motion test (CMT) in regular waves

Prescribed planar motion (constant u, v, r) of the KVLCC2 in regular waves,
solved with the time-domain linear seakeeping solver in the **ship-fixed
rotating frame**. First step of the prescribed-motion / coupled manoeuvring
programme (see `longWavelengthDiscrepancies.ipynb`, research plan).

Nothing in `manFlow`, `potForwardSpeedBC`, `linBodyMotionMj` or the existing
cases is modified. New code:

| What | Where | Built as |
|---|---|---|
| Motion + incident wave (header only) | `src/prescribedShipMotion/prescribedShipMotion.H` | – |
| Solver | `applications/solvers/manFlowPrescribed` | `manFlowPrescribed` |
| Free-surface BC | `src/finiteVolume/.../potRotatingFrameBC` | `libpotRotatingFrameBC.so` |
| Body BC with m-terms | `src/finiteVolume/.../linBodyMotionRot` | `liblinBodyMotionRot.so` |
| Near-field mean loads (existing source, now built) | `src/functionObjects/forces/nearFieldForm` | `libnearFieldForm.so` |

## Formulation

Mesh axes: x aft, y starboard, z up, bow at -x (the heading-0 geometry of the
KVLCC2 sweep); the hull never rotates in the mesh. The ship motion enters
through the frame velocity

    V_S(x) = (-u, v, 0) + Omega x (x - x_c),   Omega = (0, 0, -r)

(u, v, r MMG at midship: x forward, y starboard, r > 0 to starboard).

* **Basis flow**: absolute double-body potential, Kirchhoff split into unit
  potentials solved once: `PhiS = V0x Phi1 + V0y Phi2 + Omega_z Phi6`. The
  relative flow `Us = W = grad(PhiS) - V_S` (not irrotational with yaw), and
  `pS = -1/2(|W|^2 - |V_S|^2)`. For r = 0 this is exactly manFlow's basis flow.
* **Incident wave**: analytic, earth-fixed propagation direction, rotating in
  ship axes with psi(t); `dPhiI/dt` at a fixed mesh point uses the local
  encounter frequency `-omega + k e.V_S(x)`.
* **Free surface** (`potRotatingFrameBC`): forcing by the absolute basis
  velocity `grad(PhiS) = W + V_S` instead of `W - Uinf`; incident terms in
  both horizontal directions. Upwinding follows the local sign of W.
* **Body** (`linBodyMotionRot`, m-terms kept): tangential m-term
  `[(n.grad)W]_t = -K W + 2 n x Omega` (vorticity of the relative flow),
  `n.grad(pS)` gains `V_S.(Omega x n)`, Coriolis/centrifugal on surge and
  sway (`rotatingFrameInertia`).
* **Mean loads**: `meanLoadsNear` (near-field, instantaneous, frame
  independent) is the one to use when r != 0. `meanLoads` (midfield) is only
  valid for r = 0: with yaw rate the control-volume balance misses the storage
  and rotation terms.

The solver prints the rigid-lid added masses m11, m22, m66 from the unit
potentials at startup: compare with `addedMassDict` (Motora):
m_x' = 0.022, m_y' = 0.223, J_z' = 0.011, i.e. for L = 3.2 m, d = 0.208 m,
rho = 1000: m11 ~ 23 kg, m22 ~ 237 kg, m66 ~ 120 kg m^2.

## Running

    cp -r baseCase lam3.2_cmt_r0.08 && cd lam3.2_cmt_r0.08
    python3 runCMT.py --lam 3.2 --u 0.318 --v -0.097 --r 0.08
    python3 cmtPost.py --table <manModel case>/constant/waveData.dat --U 0.33

`runCMT.py` meshes exactly like `runsim.py` (heading 0), sets
`constant/prescribedMotion`, and picks the timestep with the largest frame
velocity near and far from the hull. `cmtPost.py` averages the loads over one
local encounter period and writes them in the waveData table axes against
the equivalent table heading, optionally with the table itself.

## Verification sequence

1. **r = 0, v = 0** (`--r 0 --v 0 --u 0.33`): must reproduce the existing
   manFlow head-sea run (`lam3.2_U0.33_head0`), both midfield and near-field.
2. **r = 0, v != 0**: must reproduce the beta runs (manFlow with
   `swaySpeed = -v`).
3. **Added masses** from the startup printout vs addedMassDict / WAMIT.
4. **r != 0**: CMT; compare with the table along the heading history
   (`cmtPost.py --table`), near-field vs midfield (the gap is the missing
   rotating-frame midfield terms).
