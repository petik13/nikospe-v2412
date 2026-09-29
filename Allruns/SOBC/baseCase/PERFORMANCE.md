# SOBC performance settings

The default is GAMG with DIC smoothing, PhiD tolerance 1e-6, relTol 0 and
nNonOrthogonalCorrectors 2. This is the user's accepted beta0_perf setup.
Two non-orthogonal corrections mean THREE solves per step. PhiS retains its
1e-11 tolerance and ten non-orthogonal corrections; startup tuning is outside
this review. No solver, gradient, timestep, mesh or output-frequency change
was applied by this performance update.

## Measured comparison

Both supplied runs finished, used 56 ranks on ggl, and have byte-identical
points, faces, owner, neighbour and boundary files (1,819,935 cells).

| Metric | beta0 | beta0_perf |
|---|---:|---:|
| Whole-run wall time, including startup | 1089 s | 477 s |
| Number of completed steps | 6276 | 4036 |
| dt during the 20–30 s comparison | 0.00816 s | 0.00831 s |
| Wall time / step over 20–30 s | 0.17796 s | 0.11388 s |
| Wall time / simulated second over 20–30 s | 21.81 s | 13.70 s |
| Mean iterations by solve over 20–30 s | 8.22 / 5.00 / 3.00 / 1.09 | 5.35 / 2.33 / 1.00 |

The mature-run step cost is about 36% lower (1.56x speedup). Total wall time
is 56% lower, but that includes fewer timesteps. beta0 starts with dt=0.00416
and switches to 0.00816 at t=18.387; beta0_perf uses 0.00831 throughout.
The final controlDict alone does not capture this baseline history. Timing
is from one run per configuration, not repeated tests under controlled load.

Over the common last-six-period window 22.404404837 .. 33.53916 s, using the
same harmonic fitting and force rotation as meanLoads.py:

| Output | beta0 | beta0_perf | Change |
|---|---:|---:|---:|
| Mean F1, N | 3.25014 | 3.20634 | −1.35% |
| Mean F2, N | 3.55991 | 3.51152 | −1.36% |
| Mean F3, N | 8.05436 | 8.49251 | +5.44% |
| Mean M6, N m | −5.56831 | −5.36358 | magnitude −3.68% |

Motion amplitudes differ by at most 0.461%; phase differences are below
0.20 degrees. Thus motions are very close, but mean loads are not identical.
These are differences between runs, not errors against an exact solution;
they cannot be assigned solely to tolerance or corrector count because dt
also differs. The user accepted the new baseline; future tuning should
still track the more sensitive mean loads.

## Is GAMG + DIC suitable?

Yes: the equation is a scalar Laplace problem and the recorded convergence
is healthy. A GAMG iteration includes work on several grid levels, smoothing
and a coarse solve. Iteration count is not a cross-smoother cost measure.
DICGaussSeidel literally runs DIC followed by Gauss-Seidel in the installed
v2412 source, so fewer cycles need not mean less elapsed time.

Relevant primary references:
[GAMG v2412 API](https://api.openfoam.com/2412/classFoam_1_1GAMGSolver.html),
[OpenCFD solver controls](https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control).
The installed source was checked directly for defaults and implementation.

## Controlled tuning sequence

Use separate prepared copies of the SAME meshed case. Keep dt fixed at the
chosen baseline value for the entire run, with the same horizon, initial
condition, output settings and rank count. The existing runsim.py cleans and
rebuilds cases, so do not use it to launch a fixed-mesh solver-only comparison.
Prepare and launch those copies using the existing decomposition workflow.
Retain a copy of the baseline fvSolution and start each candidate from it;
do not accumulate the changes below in one case.

1. Test `smoother symGaussSeidel;` against DIC. Optionally test
   `DICGaussSeidel` separately, but do not assume fewer cycles wins.
2. With DIC retained, test `nFinestSweeps 1;` (installed default is 2).
   This reduces fine-grid smoothing per cycle, but may increase cycle count.
3. With other defaults retained, test `nCellsInCoarsestLevel 20;` versus 100.
   This changes multigrid work distribution; smaller is not always faster.
4. Test `processorAgglomerator masterCoarsest;` separately. It consolidates
   the coarse-grid work and is worth measuring on 56 ranks. It is supported
   in v2412, but has not been benchmarked for this case. Keep DIC, coarse-cell
   target and sweep counts at baseline for the first comparison.
5. Test 28 versus 56 physical cores on the same reconstructed mesh, by
   re-decomposing it rather than remeshing. The current 56-rank average is
   about 32,499 cells/rank. There is no evidence yet that 28 is faster; test
   elapsed time and core-hours. Recheck motion/loads since the custom
   free-surface stencil also depends on decomposition.

Examples, applied inside an isolated candidate case:

```bash
foamDictionary system/fvSolution -entry solvers.PhiD.smoother -set symGaussSeidel
# Or, in a fresh candidate:
foamDictionary system/fvSolution -entry solvers.PhiD.nFinestSweeps -add 1
# Or, in a fresh candidate:
foamDictionary system/fvSolution -entry solvers.PhiD.nCellsInCoarsestLevel -set 20
# Or, in a fresh candidate:
foamDictionary system/fvSolution -entry solvers.PhiD.processorAgglomerator -add masterCoarsest
```

PhiS inherits PhiD entries through `$PhiD` in this template. To isolate a
transient-solver experiment completely, expand/freeze the baseline PhiS
subdictionary in the candidate before changing PhiD. Its tolerance remains
explicitly 1e-11 in all cases. No candidate above has replaced the default.

The solver does not set a final-iteration flag in this loop, so adding a
PhiDFinal dictionary alone is not an effective way to impose a tighter last
solve. Keep relTol 0 unless relative stopping is separately validated.
The custom body update also runs each corrector, so further reducing the
count changes body/fluid coupling as well as non-orthogonal correction.

Repeat promising candidates at least twice and compare a mature interval.
For a quick screening run, use the same horizon long enough to measure after
the wave ramp; validate winners over the original full duration. Compare
loads and motions before keeping a candidate.

## Read-only timing helper

From this template directory, or using the script's absolute path:

```bash
python3 comparePerformance.py /path/to/baseline /path/to/candidate --from-time 20 --to-time 30
```

It reads existing logs only. It reports timing with startup excluded, solve
iterations, all timesteps found in the log and those within the chosen
window. It rejects concatenated/restarted timing histories instead of
silently subtracting unrelated clocks. ClockTime is integer-rounded in these
logs, so compare sufficiently long windows. The script cannot determine
numerical accuracy or guarantee identical hardware load.

Full comparison data and the analysis script are retained in
`doc/manFlow-sobc-performance-assets` at the repository root.
