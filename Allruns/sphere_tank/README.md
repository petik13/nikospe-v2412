# sphere_tank: the restrained hemisphere in the towing tank of Zhao et al. (1988)

Zhao, Faltinsen, Krokstad & Aanesland (1988), "Wave-current interaction effects on
large-volume structures", BOSS '88, pp. 623-638.  Their measured mean drift force on
the restrained hemisphere (Fig. 9) is dominated by the side walls of the tank, so the
open-water `sphere` case can only be compared with their open-water theory (Fig. 4),
not with the experiment.  This case adds the tank.

## What is the same as the experiment

| | Zhao et al. | here |
|---|---|---|
| body | hemisphere D = 1 m, centre at the still water level, restrained | `sphere` geometry, `fixBody true` |
| tank breadth | 10.5 m, model on the centreline | walls at y = +-5.25 m (`front`, `back`) |
| water depth | 10 m | 10 m (`waterDepth`, mesh floor) |
| waves | regular, H = 0.1 m, along the tank | `--H 0.1` (default), heading 0 |
| current | towing against the waves, U = -0.1, 0, 0.1, 0.2 m/s | `--Ucur`, stream along the waves (U > 0: omega_e > omega_0) |
| wave given as | encounter period, 1.14 ... 1.57 s | `--Te` (lambda from omega_e = omega_0 + k U, finite depth) or `--lam` |

The walls: zero normal gradient for PhiD (the waves reflect), PhiS (the double-body
stream follows the walls) and psibar (`psiBarFarFieldPatches (inlet outlet)`: no
return flow through the walls); no side damping zone (`ySide` 1000).

## What is not modelled

The tank length (85 m), wave maker, beach and the towing transient.  The run is in the
body frame with damping zones at the inlet and outlet, i.e. an infinitely long tank.
Zhao et al.'s tank-wall theory makes the same assumption (13 image singularities on
each side).  Near the tank's natural periods (theirs, without the model: 1.16, 1.30,
1.50 s) the response converges slowly, in their theory and here; runsim.py adds
3 wall-to-wall crossings at the group velocity to endTime for the reflections.

## Running Zhao et al.'s conditions

    python3 runsim.py --Te 1.43 --Ucur 0.1        # one case (folder named as you like)

Encounter periods in the paper: 1.14, 1.17, 1.22, 1.30, 1.32, 1.43, 1.51, 1.57 s,
at U = -0.1, 0, 0.1, 0.2 m/s (U/sqrt(gD) = -0.032, 0, 0.032, 0.064).  Theory and
experiment agree best away from the tank resonances, at 1.14, 1.22, 1.43 and 1.57 s.
The authors estimate the measured drift force to +-14 %; at U = 0.2 m/s the flow may
separate (up to ~15 % viscous).

`meanLoads.py` prints and writes the Zhao et al. form (F1 / (1/2 rho g A^2 D) against
omega_e^2 R / g); the drift force along the waves is -X at heading 0.  It is the same
file as in `../sphere/baseCase`.
