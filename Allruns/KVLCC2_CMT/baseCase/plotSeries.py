#!/usr/bin/env python3
"""
Diagnostic plot of the whole time series: drift loads on the left, body
motions on the right.

Not part of prepPost.sh -- run it by hand when a case looks wrong.

The load means drawn here come from the same harmonic fits as meanLoads.py,
so for a fixed-heading run the levels marked on these axes are exactly the
numbers in the summary and in the sweep file.

Everything is in MMG axes (x forward, y starboard, z down, midship):
    loads    X > 0 forward (added resistance is X < 0), Y > 0 starboard,
             Z > 0 down, N about z down (> 0 bow to starboard).  The function
             objects write in mesh axes; with h = -headingAngle (0 for the
             circular motion test, where the hull is aligned with the mesh,
             bow at -x):
                 X = -(F_x cos h - F_y sin h),  Y = F_x sin h + F_y cos h,
                 Z = -F_z,  N = -M_z
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

Non-dimensional by default: loads / (rho g A^2 B^2 / L), moments
/ (rho g A^2 B^2); motions eta1..3 / A and eta4..6 (rad) / (k A), as the
eta columns of results.csv.  --raw gives N, N m, m and deg.

Motions are drawn with their envelope (red: running max / min over one local
encounter period) and the half peak-to-peak value (black dashed), which is
what the eta columns of results.csv measure.

The contributions panel splits the F1 of the chosen load object into the
terms it writes, read from the force.dat header: surface, elevation and strip
(Chen's midfield, middleFieldForm), plus coriolis, storage and centripetal for
the rotating midfield (middleFieldFormRot); the legend gives the mean of each
over the averaging window.

Circular motion test (constant/prescribedMotion present): the local
encounter period follows the heading from the trajectory, the top axes give
the heading (unwrapped along the turn, shown mod 360), and dotted
lines mark the yaw-rate onset and the end of its ramp.

usage:
    python3 plotSeries.py                     # opens a window
    python3 plotSeries.py --save              # also write timeSeries.png
    python3 plotSeries.py --save foo.pdf
    python3 plotSeries.py --fo meanLoadsRot   # load object (default meanLoads)
    python3 plotSeries.py --tmin 4 --tmax 10  # zoom, e.g. onto a blow-up
    python3 plotSeries.py --raw               # loads in N, motions in m/deg
"""

import argparse
import glob
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

from meanLoads import RHO, G, read_dict, load, fit


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


def level(ax, value, fmt="{:+.4f}"):
    ax.axhline(value, color="C3", lw=1.0, ls="--")
    ax.annotate(fmt.format(value), xy=(0.995, value), xycoords=("axes fraction",
                "data"), ha="right", va="bottom", fontsize=8, color="C3")


def header_columns(pattern):
    """{name: column} of the vectors in a force.dat / moment.dat, from its
    '# Time ...' header line (column of the x component; Time is column 0)."""
    files = sorted(glob.glob(pattern))
    if not files:
        return {}
    cols = {}
    with open(files[0]) as f:
        for line in f:
            if not line.startswith("#"):
                break
            tok = line[1:].split()
            if tok and tok[0] == "Time":
                for i, name in enumerate(tok):
                    if name.endswith("_x"):
                        cols[name[:-2]] = i
    return cols


def running_mean(t, y, w, width, npts=250):
    """Harmonic-fit mean over a sliding window of the given width."""
    lo = t[0] + width
    if t[-1] <= lo:
        return None, None
    centres = np.linspace(lo, t[-1], npts)
    out = []
    for c in centres:
        m = (t >= c - width) & (t <= c)
        out.append(fit(t[m], y[m], w, 2)[0] if m.sum() >= 20 else np.nan)
    return centres, np.array(out)


def envelope(t, y, T, npts=800):
    """Upper and lower envelope of y through its peaks: the maximum / minimum
    in a centred window of one local period T(t) (array on the t grid) picks
    the peaks, and the envelope interpolates linearly between them."""
    tc = np.linspace(t[0], t[-1], npts)
    Tc = np.interp(tc, t, T)
    imax, imin = set(), set()
    for c, w in zip(tc, Tc):
        i0 = np.searchsorted(t, c - 0.5 * w)
        i1 = np.searchsorted(t, c + 0.5 * w, side="right")
        if i1 - i0 >= 3 and c - 0.5 * w >= t[0] and c + 0.5 * w <= t[-1]:
            imax.add(i0 + int(np.argmax(y[i0:i1])))
            imin.add(i0 + int(np.argmin(y[i0:i1])))
    if len(imax) < 2 or len(imin) < 2:
        return tc, np.full(npts, np.nan), np.full(npts, np.nan)
    imax, imin = np.array(sorted(imax)), np.array(sorted(imin))
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
    ap.add_argument("--nper", type=int, default=6,
                    help="encounter periods in the averaging window")
    ap.add_argument("--nharm", type=int, default=3)
    ap.add_argument("--fo", default="meanLoads")
    ap.add_argument("--save", nargs="?", const="timeSeries.png", default=None,
                    metavar="FILE", help="also write the figure to a file")
    ap.add_argument("--raw", action="store_true",
                    help="plot loads in N and motions in m/deg")
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=None)
    a = ap.parse_args()

    wc = read_dict("constant/waveConditions")
    bd = read_dict("constant/bodyMotionProperties")

    lam, steep = wc["waveLength"], wc["steepness"]
    U0, head, h = wc["currentSpeed"], wc["headingAngle"], wc["waterDepth"]
    L = bd.get("Lpp", 1.0)
    B = bd.get("beam", 1.0)

    k = 2.0 * np.pi / lam
    w0 = np.sqrt(G * k * np.tanh(k * h))
    we = w0 + k * U0 * np.cos(head)
    Te = 2.0 * np.pi / we
    A = 0.5 * steep * lam
    t_ramp = wc.get("rampPeriods", 0.0) * 2.0 * np.pi / w0

    den_F = 1.0 if a.raw else RHO * G * A**2 * B**2 / L
    den_M = 1.0 if a.raw else RHO * G * A**2 * B**2

    F = load(f"postProcessing/{a.fo}/*/force.dat")
    M = load(f"postProcessing/{a.fo}/*/moment.dat")
    mot = load("postProcessing/bodyMotion/motion.dat")

    if F is None and mot is None:
        sys.exit("ERROR: no postProcessing output found -- run from the case "
                 "directory")

    tmax = max(d[-1, 0] for d in (F, M, mot) if d is not None)

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

    nper = min(a.nper, max(int(np.floor((tmax - t_ramp) / Te)), 1))
    lo, hi = tmax - nper * Te, tmax

    fig, ax = plt.subplots(6, 2, figsize=(12, 12.5), sharex=True,
                           constrained_layout=True)

    case = os.path.basename(os.getcwd())
    fig.suptitle(
        f"{case}      lambda/L {lam/L:.3g}   A {A:.4g} m   "
        f"U {U0:+.4g} m/s (Fn {U0/np.sqrt(G*L):.3g})   "
        f"{head_lbl}   "
        f"Te {Te:.3g} s   ramp {t_ramp:.3g} s\n"
        f"loads: grey = wave ramp,  green = averaging window (last {nper} Te),  "
        f"red dashed = fitted mean"
        + (",  dotted = yaw onset / end of ramp" if onset else "") + "\n"
        f"motions: red = envelope,  black dashed = half peak-to-peak",
        fontsize=10)

    # --- left column: drift loads, MMG axes --------------------------------
    # mesh -> MMG: h = -headingAngle (0 for the CMT, hull bow at -x)
    hcos, hsin = np.cos(-head), np.sin(-head)
    mmgX = lambda fx, fy: -(fx * hcos - fy * hsin)
    mmgY = lambda fx, fy: fx * hsin + fy * hcos

    unit = " [N]" if a.raw else ""
    if F is not None:
        t = F[:, 0]
        m = (t >= lo) & (t <= hi)
        FX, FY, FZ = mmgX(F[:, 1], F[:, 2]), mmgY(F[:, 1], F[:, 2]), -F[:, 3]
        rows = [("X  surge (> 0 fwd)", FX, den_F),
                ("Y  sway (> 0 stbd)", FY, den_F),
                ("Z  heave (> 0 down)", FZ, den_F)]
        for i, (lbl, y, den) in enumerate(rows):
            panel(ax[i][0], t, y / den, lbl + unit)
            mark(ax[i][0], t_ramp, lo, hi, onset)
            level(ax[i][0], fit(t[m], y[m], 2 * we, a.nharm)[0] / den)

        NZ = None
        if M is not None:
            tm = M[:, 0]
            mm = (tm >= lo) & (tm <= hi)
            NZ = -M[:, 3]
            panel(ax[3][0], tm, NZ / den_M,
                  "N  yaw (> 0 bow stbd)" + (" [N m]" if a.raw else ""))
            mark(ax[3][0], t_ramp, lo, hi, onset)
            level(ax[3][0], fit(tm[mm], NZ[mm], 2 * we, a.nharm)[0] / den_M)

        # the total is a difference of larger terms -- show them, with their
        # means over the averaging window.  Columns from the header:
        # middleFieldForm  total surface elevation strip
        # middleFieldFormRot  total chen surface elevation strip coriolis
        #                     storage [centripetal]
        hcols = header_columns(f"postProcessing/{a.fo}/*/force.dat")
        terms = (("surface", "C0"), ("elevation", "C1"), ("strip", "C4"),
                 ("coriolis", "C5"), ("storage", "C8"), ("centripetal", "C6"))
        for name, c in terms + (("total", "k"),):
            col = hcols.get(name, 1 if name == "total" else None)
            if col is None or col + 1 >= F.shape[1]:
                continue
            yx = mmgX(F[:, col], F[:, col + 1])
            mean = fit(t[m], yx[m], 2 * we, a.nharm)[0] / den_F if m.sum() >= 20 else np.nan
            ax[4][0].plot(t, yx / den_F, c, lw=1.2 if name == "total" else 0.8,
                          label=f"{name} {mean:+.3g}")
        ax[4][0].set_ylabel("X contributions" + unit)
        ax[4][0].grid(alpha=0.3)
        ax[4][0].legend(fontsize=7, ncol=4, loc="best",
                        title="mean over the averaging window", title_fontsize=7)
        mark(ax[4][0], t_ramp, lo, hi, onset)

        # Convergence of the mean: a sliding 2-Te fit.  In head seas F2 and F6
        # are three orders down on F1, so they get their own axis or they plot
        # flat.
        tw = ax[5][0].twinx()
        for axis, series in ((ax[5][0], ((t, FX, den_F, "X", "C0"),)),
                             (tw, ((t, FY, den_F, "Y", "C1"),
                                   (M[:, 0] if M is not None else None, NZ, den_M, "N", "C4")))):
            for tt, y, den, lbl, c in series:
                if y is None:
                    continue
                tc, yc = running_mean(tt, y / den, 2 * we, 2 * Te)
                if tc is not None:
                    axis.plot(tc, yc, c, lw=1.2, label=lbl)
        ax[5][0].set_ylabel("running mean (2 Te): X" + unit, color="C0")
        tw.set_ylabel("Y, N", color="C1")
        ax[5][0].grid(alpha=0.3)
        ax[5][0].axhline(0, color="k", lw=0.5, alpha=0.4)
        h1, l1 = ax[5][0].get_legend_handles_labels()
        h2, l2 = tw.get_legend_handles_labels()
        ax[5][0].legend(h1 + h2, l1 + l2, fontsize=8, ncol=3,
                        loc="lower right", framealpha=0.9)
        mark(ax[5][0], t_ramp, lo, hi, onset)
    else:
        for i in range(6):
            ax[i][0].text(0.5, 0.5, "no drift-load output", ha="center",
                          transform=ax[i][0].transAxes)

    # --- right column: body motions (MMG) with their envelope --------------
    # non-dimensional (eta1..3 / A, eta4..6 / kA, as results.csv) unless --raw.
    # The solver's body frame is x forward, y port, z up: MMG sign per dof.
    dofs = [("surge", 1, 1.0, "m", A, "A", 1.0), ("sway", 2, 1.0, "m", A, "A", -1.0),
            ("heave", 3, 1.0, "m", A, "A", -1.0),
            ("roll", 4, np.degrees(1), "deg", k * A, "kA", 1.0),
            ("pitch", 5, np.degrees(1), "deg", k * A, "kA", -1.0),
            ("yaw", 6, np.degrees(1), "deg", k * A, "kA", -1.0)]

    if mot is not None and mot.shape[1] >= 7:
        t = mot[:, 0]
        still = np.allclose(mot[:, 1:7], 0)
        T_loc = Te_at(t)
        for i, (name, col, deg, u, norm, nlbl, sgn) in enumerate(dofs):
            y = sgn * mot[:, col] * (deg if a.raw else 1.0 / norm)
            ylab = f"{name} [{u}]" if a.raw else f"eta{col} = {name} / {nlbl}"
            panel(ax[i][1], t, y, ylab, colour="C1", lw=0.6)
            mark(ax[i][1], t_ramp, onset=onset)
            if not still:
                tc, up, dn = envelope(t, y, T_loc)
                ax[i][1].plot(tc, up, "C3", lw=1.1, label="envelope")
                ax[i][1].plot(tc, dn, "C3", lw=1.1)
                ax[i][1].plot(tc, 0.5 * (up - dn), "k--", lw=1.0, label="half peak-to-peak")
        if still:
            ax[0][1].set_title("body restrained (fixBody true)", fontsize=9,
                               loc="left", color="C3")
        else:
            ax[0][1].legend(fontsize=7, loc="upper left")
    else:
        for i in range(6):
            ax[i][1].text(0.5, 0.5, "no motion output", ha="center",
                          transform=ax[i][1].transAxes)

    for j in range(2):
        ax[5][j].set_xlabel("t [s]")
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
