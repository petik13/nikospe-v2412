# KVLCC2 circular motion test (CMT) in regular waves

Prescribed planar motion (u, v, r) of the KVLCC2 in regular waves, solved with
the time-domain linear seakeeping solver in the **ship-fixed rotating frame**.
First step of the prescribed-motion / coupled manoeuvring programme (see
`longWavelengthDiscrepancies.ipynb`, research plan and the section "Midfield
mean wave loads in a rotating ship frame").

Nothing in `manFlow`, `potForwardSpeedBC`, `linBodyMotionMj`,
`middleFieldForm` or the existing cases is modified. New code:

| What | Where | Built as |
|---|---|---|
| Motion + incident wave (header only) | `src/prescribedShipMotion/prescribedShipMotion.H` | – |
| Solver | `applications/solvers/manFlowPrescribed` | `manFlowPrescribed` |
| Free-surface BC | `src/finiteVolume/.../potRotatingFrameBC` | `libpotRotatingFrameBC.so` |
| Body BC with m-terms | `src/finiteVolume/.../linBodyMotionRot` | `liblinBodyMotionRot.so` |
| Midfield mean loads, rotating control volume | `src/functionObjects/forces/middleFieldFormRot` | `libmiddleFieldFormRot.so` |
| Near-field mean loads (existing source, now built) | `src/functionObjects/forces/nearFieldForm` | `libnearFieldForm.so` |

All are in `Allmake.sh`.

## Formulation

Mesh axes: x aft, y starboard, z up, bow at -x (the heading-0 geometry of the
KVLCC2 sweep); the hull never rotates in the mesh. The ship motion enters
through the frame velocity

    V_S(x, t) = (-u, v, 0) + Omega(t) x (x - x_c),   Omega = (0, 0, -r(t))

(u, v, r MMG at midship: x forward, y starboard, r > 0 to starboard).

* **Yaw-rate onset**: the ship runs straight for `--yawOnsetPeriods`
  (default 10) encounter periods so the motions settle, then r is ramped in
  (half-cosine) over `--yawRampPeriods` (default 2). During the ramp the basis
  flow is rebuilt every step (quasi-steady) and the boundary conditions
  rebuild their m-terms, steady-flow restoring and upwind schemes. The mean
  loads are output from the onset; during the ramp they are transitional
  (shaded in the `cmtPost.py` plot).
* **Basis flow**: absolute double-body potential, Kirchhoff split into unit
  potentials solved once: `PhiS = V0x Phi1 + V0y Phi2 + Omega_z Phi6`. The
  relative flow `Us = W = grad(PhiS) - V_S` (not irrotational with yaw), and
  `pS = -1/2(|W|^2 - |V_S|^2)`. For r = 0 this is exactly manFlow's basis flow.
* **Incident wave**: analytic, earth-fixed propagation direction, rotating in
  ship axes with psi(t); `dPhiI/dt` at a fixed mesh point uses the local
  encounter frequency `-omega + k e.V_S(x, t)`.
* **Free surface** (`potRotatingFrameBC`): forcing by the absolute basis
  velocity `grad(PhiS) = W + V_S` instead of `W - Uinf`; incident terms in
  both horizontal directions. Upwinding follows the local sign of W.
* **Body** (`linBodyMotionRot`, m-terms kept): tangential m-term
  `[(n.grad)W]_t = -K W + 2 n x Omega` (vorticity of the relative flow),
  `n.grad(pS)` gains `V_S.(Omega x n)`, Coriolis/centrifugal on surge and
  sway (`rotatingFrameInertia`).
* **Mean loads**: `meanLoadsRot` (`middleFieldFormRot`) = Chen's midfield
  (`middleFieldForm`, unchanged) **plus the storage and Coriolis terms** of a
  control volume rotating with the ship,

      F = F_Chen[W] - dP/dt - 2 Omega x P,   Mz = Mz_Chen[W] - dH/dt - 2 Omega Q

  with P the second-order wave momentum in the control volume (free-surface
  strip + hull displacement). Without them the midfield is wrong by about
  2.2 rho g A^2 B^2/L in the KVLCC2 CMT at lambda/L = 1 (the incident wave's
  momentum in the rotating control volume). `meanLoadsRot2` is the same on
  the larger control volume (+-2 L/2); the two must agree. `meanLoads`
  (uncorrected) and `meanLoadsNear` (near-field, not reliable: 24% off for a
  restrained hull at zero speed, ~2x with motions or speed) are kept for
  comparison.

The solver prints the rigid-lid added masses m11, m22, m66 from the unit
potentials at startup: compare with `addedMassDict` (Motora):
m_x' = 0.022, m_y' = 0.223, J_z' = 0.011, i.e. for L = 3.2 m, d = 0.208 m,
rho = 1000: m11 ~ 23 kg, m22 ~ 237 kg, m66 ~ 120 kg m^2 (measured on the
KVLCC2 mesh: 19, 260, 139).

## Running

    cp -r baseCase lam3.2_cmt_r0.08 && cd lam3.2_cmt_r0.08
    python3 runCMT.py --lam 3.2 --u 0.318 --v -0.097 --r 0.08
    python3 cmtPost.py --table <manModel case>/constant/waveData.dat --U 0.33

`runCMT.py` meshes exactly like `runsim.py` (heading 0), builds both control
volumes, sets `constant/prescribedMotion` (including the yaw onset), and
picks the timestep with the largest frame velocity near and far from the
hull. Default end time for r != 0: onset + ramp/2 + one full turn (2 pi/|r|)
+ the half-width of the mean-load filter (2 Te_max), so that `cmtPost.py`
gives mean loads for a full turn, starting at the yaw onset (heading change 0).
The yaw-rate ramp is shaded in the plot (dOmega/dt != 0 there, not in the
midfield formula).
`cmtPost.py` low-pass filters the loads (Gaussian, from the yaw onset) and
reports everything in MMG axes (x forward, y starboard, z down, midship):
F1 = X > 0 forward (added resistance is X < 0), F2 = Y > 0 starboard,
Mz = N > 0 bow to starboard, against the encounter angle
mu = (waveDirection + 180 - psi) mod 360 (0 head sea, 90 waves from
starboard), as manModel. The function objects write in mesh axes (x aft,
y starboard, z up): X = -F_x, Y = F_y, N = -M_z. `results.csv` (manModel
waveData format, with the `# convention: MMG` marker line) has one row per
10 deg of mu along the turn. A table given with `--table` is compared as it
is if it has the marker, and converted from the old convention (heading 90 =
waves from port, F1 > 0 added resistance) if not. It also prints the
control-volume dependence (max |rot - rot2|). `plotSeries.py` and
`meanLoads.py` report in the same axes.

## Verification sequence

1. **r = 0, v = 0** (`--r 0 --v 0 --u 0.33`): reproduces the manFlow head-sea
   run (done: added resistance 1.90 against the table's 1.95, i.e.
   X = -1.90 in MMG axes). `meanLoadsRot` total must
   equal `meanLoads` here.
2. **r = 0, v != 0**: must reproduce the beta runs (manFlow with
   `swaySpeed = -v`).
3. **Added masses** from the startup printout vs addedMassDict / WAMIT (done).
4. **r != 0**: CMT. `meanLoadsRot` vs `meanLoadsRot2` (control-volume
   independence), then against the straight-course table along the heading
   history (`cmtPost.py --table`).
