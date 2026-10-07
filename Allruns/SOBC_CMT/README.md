# SOBC circular motion test (CMT) in regular waves

Same setup as `KVLCC2_CMT` (see `KVLCC2_CMT/README.md` for the method, the
yaw onset and ramp, the rotating-frame midfield and the post-processing), with
the SOBC hull. Built from `KVLCC2_CMT/baseCase` and `SOBC/baseCase` by
`Allruns/makeSOBC_CMT.py`.

Hull-specific (from `SOBC/baseCase`): `constant/triSurface/SOBC.stl`,
`extendedFeatureEdgeMesh`, `snappyHexMeshDict`, `surfaceFeatureExtractDict`,
`bodyMotionProperties` (SOBC mass, inertia, restoring, mooring; plus
`rotatingFrameInertia true`), `linMotions.py`, `meanLoads.py`, `meanVal.py`.

`runCMT.py`: scale 32, Lpp 5.9375 m, B 1.00625 m, draft 0.34375 m; near-hull
refinement boxes as `SOBC/baseCase/runsim.py` (l3 = max(2.2 L/2, lambda);
boxes 4-7 at 2.5, 1.4, 1.3, 1.2 L/2); control volumes +-1.2 L/2 and +-2 L/2;
default 52 processors.

`cmtPost.py` is the same as in `KVLCC2_CMT`: the mean loads are in the same
axes for both hulls (F1 > 0 added resistance, F2 > 0 to starboard, Mz > 0 bow
to port), the convention of the existing SOBC and KVLCC2 waveData tables.
Note: `SOBC/baseCase/meanLoads.py` (and its copy here) currently divides the
forces by -rho g A^2 B^2/L, which flips F1 and F2 relative to those tables.

usage, from a copy of `baseCase`:

    python3 runCMT.py --lam 5.9375 --u 0.3086 --v 0 --r 0.05
    python3 cmtPost.py --table ../waveData.dat      # SOBC straight-course table

`cmtPost.py` writes `results.csv` (lam/L, U, V, heading, F1mean, F2mean,
Mzmean, eta1..eta6) for every 10 deg of heading along the turn.
