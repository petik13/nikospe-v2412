#!/usr/bin/env python3
"""
Diagnostic plot of the motion RAOs: the six body motions (surge, sway,
heave, roll, pitch, yaw) in a 3x2 figure (translational dofs in the left
column, rotational in the right).

Not part of prepPost.sh -- run it by hand when a case looks wrong.

Everything is in MMG axes (x forward, y starboard, z down, midship):
    motions  surge > 0 forward, sway > 0 starboard, heave > 0 down, roll,
             pitch, yaw about x forward, y starboard, z down.  The solver
             writes them in its body frame (x forward, y port, z up), so
             sway, heave, pitch and yaw change sign.
    angle    heading = psi - waveDirection + 180 = 360 - mu: 0 head seas,
             increasing as the ship turns to starboard, 90 waves from port,
             270 waves from starboard (the runsim --heading convention; with
             the runCMT.py defaults psi0 = 0, waveDirection = 180 it is the
             ship heading psi).  The tables (results.csv) use the encounter
             angle mu = 360 - heading instead.

Non-dimensional by default: motions eta1..3 / A and eta4..6 (rad) / (k A),
as the eta columns of results.csv.  --raw gives m and deg.

Motions are drawn with their envelope (red: running max / min over one
local encounter period) and the half peak-to-peak value (black dashed),
which is what the eta columns of results.csv measure.

Circular motion test (constant/prescribedMotion present): the local
encounter period follows the heading from the trajectory, the top axes give
the heading (unwrapped along the turn, shown mod 360), and dotted
lines mark the yaw-rate onset and the end of its ramp.

usage:
    python3 plotMotions.py                     # opens a window
    python3 plotMotions.py --save              # also write timeSeries.png
    python3 plotMotions.py --save foo.pdf
    python3 plotMotions.py --tmin 4 --tmax 10  # zoom, e.g. onto a blow-up
    python3 plotMotions.py --raw               # motions in m/deg
"""

import argparse
import os
import sys

import numpy as np

import matplotlib

# Opens a window by default; falls back to writing a file where there is no
# display (over ssh, on a cluster node).  The backend must be chosen before
# pyplot is imported.
HEADLESS = not (os.environ.get("DISPLAY") or sys.platform == "darwin")
if HEADLESS:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

from meanLoads import G, read_dict, load


# --------------------------------------------------------------- helpers ----
def panel(ax, t, y, label, colour="C0", lw=0.9):
    ax.plot(t, y, colour, lw=lw)
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)
    ax.axhline(0, color="k", lw=0.5, alpha=0.4)


def mark(ax, t_ramp, lo=None, hi=None, onset=()):
    """Shade the wave ramp and (optionally) the averaging window; dotted
    lines at the yaw-rate onset and the end of its ramp."""
    if t_ramp > 0:
        ax.axvspan(0, t_ramp, color="0.85", zorder=0)
    if lo is not None:
        ax.axvspan(lo, hi, color="C2", alpha=0.10, zorder=0)
    for tt in onset:
        ax.axvline(tt, color="k", lw=0.8, ls=":", alpha=0.7)


def envelope(t, y, npts=800):
    """Upper and lower envelope of y: every local maximum and every local
    minimum of the raw samples, connected by straight lines."""
    imax = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    imin = np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1
    if len(imax) < 2 or len(imin) < 2:
        return t, np.full(t.shape, np.nan), np.full(t.shape, np.nan)
    lo_, hi_ = max(t[imax[0]], t[imin[0]]), min(t[imax[-1]], t[imin[-1]])
    tc = np.linspace(lo_, hi_, npts)
    return tc, np.interp(tc, t[imax], y[imax]), np.interp(tc, t[imin], y[imin])


def heading_axis(ax, cmt, step=30.0):
    """Top axis with the heading at the times it is reached."""
    t, h = cmt["t"], cmt["h"]

    sel = t >= cmt["tOn"]
    if sel.sum() < 2:
        return
    t, h = t[sel], h[sel]
    sgn = 1.0 if h[-1] >= h[0] else -1.0
    hh = np.maximum.accumulate(sgn * h)                  # increasing along the turn
    targets = np.arange(np.ceil(hh[0] / step - 1e-6) * step, hh[-1] + 1e-9, step)
    ticks = [t[np.searchsorted(hh, x - 1e-9)] for x in targets]
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    lo, hi = ax.get_xlim()
    keep = [(tk, x) for tk, x in zip(ticks, targets) if lo <= tk <= hi]
    top.set_xticks([tk for tk, _ in keep])
    top.set_xticklabels([f"{(sgn * x) % 360:.0f}" for _, x in keep], fontsize=8)
    top.set_xlabel("heading [deg] (0 head seas, 90 waves from port)", fontsize=9)


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", nargs="?", const="timeSeries.png", default=None,
                    metavar="FILE", help="also write the figure to a file")
    ap.add_argument("--raw", action="store_true",
                    help="plot motions in m/deg instead of non-dimensional")
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)

    a = ap.parse_args()

    wc = read_dict("constant/waveConditions")
    bd = read_dict("constant/bodyMotionProperties")

    waveDat = np.loadtxt("../waveData10_u033.dat", delimiter=",", skiprows=1)

    lam, steep = wc["waveLength"], wc["steepness"]
    U0, head, h = wc["currentSpeed"], wc["headingAngle"], wc["waterDepth"]
    L = bd.get("Lpp", 1.0)

    k = 2.0 * np.pi / lam
    w0 = np.sqrt(G * k * np.tanh(k * h))
    we = w0 + k * U0 * np.cos(head)
    Te = 2.0 * np.pi / we
    A = 0.5 * steep * lam
    t_ramp = wc.get("rampPeriods", 0.0) * 2.0 * np.pi / w0

    mot = load("postProcessing/bodyMotion/motion.dat")

    if mot is None:
        sys.exit("ERROR: no postProcessing output found -- run from the case "
                 "directory")

    tmax = mot[-1, 0]

    # Local encounter period: constant for a fixed-heading run.  Circular
    # motion test: follows the heading psi(t) of the trajectory, as
    # cmtPost.py: omega_e = |-omega + k e.V0|, e = (-cos a, sin a), a = chi - psi;
    # the load means are fitted at its value at the END of the run.
    Te_at = lambda tt: np.full_like(np.asarray(tt, dtype=float), Te)
    cmt, onset = None, ()
    tr = load("postProcessing/shipTrajectory/trajectory.dat")
    if os.path.isfile("constant/prescribedMotion") and tr is not None:
        pm = read_dict("constant/prescribedMotion")
        chi = pm.get("waveDirection", 180.0)
        uu, vv = pm.get("u", 0.0), pm.get("v", 0.0)

        def Te_at(tt):
            al = np.radians(chi - np.interp(tt, tr[:, 0], tr[:, 1]))
            return 2.0 * np.pi / np.abs(-w0 + k * (np.cos(al) * uu + np.sin(al) * vv))

        tOn, tYr = pm.get("yawOnsetTime", 0.0), pm.get("yawRampTime", 0.0)
        cmt = {"t": tr[:, 0], "h": tr[:, 1] - chi + 180.0, "tOn": tOn}   # heading, unwrapped
        onset = (tOn, tOn + tYr) if tYr > 0 else (tOn,)
        Te = float(Te_at(tmax))
        we = 2.0 * np.pi / Te
        hd_end = (np.interp(tmax, tr[:, 0], tr[:, 1]) - chi + 180.0) % 360.0
        head_lbl = f"CMT, at the end: heading {hd_end:.1f} deg"
    else:
        hd0 = (-np.degrees(head)) % 360.0 + 0.0      # straight run: heading = -headingAngle
        head_lbl = f"heading {hd0:.4g} deg"

    fig, ax = plt.subplots(3, 2, figsize=(15, 12.5), sharex=True,
                           constrained_layout=True)

    case = os.path.basename(os.getcwd())
    fig.suptitle(
        f"{case}      lambda/L {lam/L:.3g}   A {A:.4g} m   "
        f"U {U0:+.4g} m/s (Fn {U0/np.sqrt(G*L):.3g})   "
        f"{head_lbl}   "
        f"Te {Te:.3g} s   ramp {t_ramp:.3g} s\n"
        f"grey = wave ramp"
        + (",  dotted = yaw onset / end of ramp" if onset else "") + "\n"
        f"motions: red = envelope,  black dashed = half peak-to-peak,  "
        f"green = steady-state reference",
        fontsize=10)

    # --- motion RAOs (MMG) with their envelope ------------------------------
    # non-dimensional (eta1..3 / A, eta4..6 / kA, as results.csv) unless --raw.
    # The solver's body frame is x forward, y port, z up: MMG sign per dof.
    # Translational dofs (surge, sway, heave) in the left column, rotational
    # (roll, pitch, yaw) in the right.
    dofs = [("surge", 1, 1.0, "m", A, "A", 1.0), ("sway", 2, 1.0, "m", A, "A", -1.0),
            ("heave", 3, 1.0, "m", A, "A", -1.0),
            ("roll", 4, np.degrees(1), "deg", k * A, "kA", 1.0),
            ("pitch", 5, np.degrees(1), "deg", k * A, "kA", -1.0),
            ("yaw", 6, np.degrees(1), "deg", k * A, "kA", -1.0)]

    if mot.shape[1] >= 7:
        t = mot[:, 0]
        still = np.allclose(mot[:, 1:7], 0)
        for i, (name, col, deg, u, norm, nlbl, sgn) in enumerate(dofs):
            row, cidx = i % 3, i // 3
            axi = ax[row][cidx]
            y = sgn * mot[:, col] * (deg if a.raw else 1.0 / norm)
            ylab = f"{name} [{u}]" if a.raw else f"eta{col} = {name} / {nlbl}"
            panel(axi, t, y, ylab, colour="C1", lw=0.6)
            mark(axi, t_ramp, onset=onset)
            if not still:
                tc, up, dn = envelope(t, y)
                axi.plot(tc, up, "C3", lw=1.1, label="envelope")
                axi.plot(tc, dn, "C3", lw=1.1)
                t_h = []
                for tab_h in waveDat[:, 4]:
                    t_h.append(np.interp(tab_h, cmt["h"], cmt["t"]))
                axi.plot(t_h, waveDat[:, 8 + i], lw=1.5, color="C2", label="steady-state")
        if still:
            ax[0][0].set_title("body restrained (fixBody true)", fontsize=9,
                               loc="left", color="C3")
        else:
            ax[0][0].legend(fontsize=7, loc="upper left")
    else:
        for i in range(6):
            ax[i % 3][i // 3].text(0.5, 0.5, "no motion output", ha="center",
                          transform=ax[i % 3][i // 3].transAxes)

    for j in range(2):
        ax[2][j].set_xlabel("t [s]")
    if a.tmin is not None or a.tmax is not None:
        ax[0][0].set_xlim(a.tmin, a.tmax)
    if cmt is not None:
        for j in range(2):
            heading_axis(ax[0][j], cmt)

    out = a.save or ("timeSeries.png" if HEADLESS else None)
    if out:
        fig.savefig(out, dpi=110)
        print(f"wrote {out}")
    if HEADLESS:
        print("no display -- wrote a file instead of opening a window")
    else:
        plt.show()


if __name__ == "__main__":
    main()
