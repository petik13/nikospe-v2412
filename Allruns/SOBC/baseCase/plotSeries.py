#!/usr/bin/env python3
"""
Diagnostic plot of the whole time series: drift loads on the left, body
motions on the right.

Not part of prepPost.sh -- run it by hand when a case looks wrong.

The means and RAOs drawn here come from the same harmonic fits as
meanLoads.py, so the levels marked on these axes are exactly the numbers in
the summary and in the sweep file.

usage:
    python3 plotSeries.py                     # opens a window
    python3 plotSeries.py --save              # also write timeSeries.png
    python3 plotSeries.py --save foo.pdf
    python3 plotSeries.py --tmin 4 --tmax 10  # zoom, e.g. onto a blow-up
    python3 plotSeries.py --raw               # loads in N, motions in m/deg
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

from meanLoads import RHO, G, read_dict, load, fit


# --------------------------------------------------------------- helpers ----
def panel(ax, t, y, label, colour="C0"):
    ax.plot(t, y, colour, lw=0.9)
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)
    ax.axhline(0, color="k", lw=0.5, alpha=0.4)


def mark(ax, t_ramp, lo, hi):
    """Shade the ramp and the averaging window."""
    if t_ramp > 0:
        ax.axvspan(0, t_ramp, color="0.85", zorder=0)
    ax.axvspan(lo, hi, color="C2", alpha=0.10, zorder=0)


def level(ax, value, fmt="{:+.4f}"):
    ax.axhline(value, color="C3", lw=1.0, ls="--")
    ax.annotate(fmt.format(value), xy=(0.995, value), xycoords=("axes fraction",
                "data"), ha="right", va="bottom", fontsize=8, color="C3")


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
    nper = min(a.nper, max(int(np.floor((tmax - t_ramp) / Te)), 1))
    lo, hi = tmax - nper * Te, tmax

    fig, ax = plt.subplots(6, 2, figsize=(15, 15), sharex=True,
                           constrained_layout=True)

    case = os.path.basename(os.getcwd())
    fig.suptitle(
        f"{case}      lambda/L {lam/L:.3g}   A {A:.4g} m   "
        f"U {U0:+.4g} m/s (Fn {U0/np.sqrt(G*L):.3g})   "
        f"heading {np.degrees(head)+0.0:+.4g} deg   "
        f"Te {Te:.3g} s   ramp {t_ramp:.3g} s\n"
        f"grey = ramp,  green = averaging window (last {nper} Te),  "
        f"red dashed = fitted mean / amplitude",
        fontsize=11)

    # --- left column: drift loads ----------------------------------------
    unit = " [N]" if a.raw else ""
    if F is not None:
        t = F[:, 0]
        m = (t >= lo) & (t <= hi)
        rows = [("F1  surge (added R)", 1, den_F),
                ("F2  sway", 2, den_F),
                ("F3  heave", 3, den_F)]
        for i, (lbl, col, den) in enumerate(rows):
            panel(ax[i][0], t, F[:, col] / den, lbl + unit)
            mark(ax[i][0], t_ramp, lo, hi)
            level(ax[i][0], fit(t[m], F[m, col], 2 * we, a.nharm)[0] / den)

        if M is not None:
            tm = M[:, 0]
            mm = (tm >= lo) & (tm <= hi)
            panel(ax[3][0], tm, M[:, 3] / den_M,
                  "F6  yaw" + (" [N m]" if a.raw else ""))
            mark(ax[3][0], t_ramp, lo, hi)
            level(ax[3][0], fit(tm[mm], M[mm, 3], 2 * we, a.nharm)[0] / den_M)

        # the total is a difference of larger terms -- show them
        for off, lbl, c in ((3, "surface", "C0"), (6, "elevation", "C1"),
                            (9, "strip", "C4")):
            ax[4][0].plot(t, F[:, 1 + off] / den_F, c, lw=0.9, label=lbl)
        ax[4][0].plot(t, F[:, 1] / den_F, "k", lw=1.2, label="total")
        ax[4][0].set_ylabel("F1 contributions" + unit)
        ax[4][0].grid(alpha=0.3)
        ax[4][0].legend(fontsize=8, ncol=4, loc="best")
        mark(ax[4][0], t_ramp, lo, hi)

        # Convergence of the mean: a sliding 2-Te fit.  In head seas F2 and F6
        # are three orders down on F1, so they get their own axis or they plot
        # flat.
        tw = ax[5][0].twinx()
        for axis, series in ((ax[5][0], ((1, F, den_F, "F1", "C0"),)),
                             (tw, ((2, F, den_F, "F2", "C1"),
                                   (3, M, den_M, "F6", "C4")))):
            for col, d, den, lbl, c in series:
                if d is None:
                    continue
                tc, yc = running_mean(d[:, 0], d[:, col] / den, 2 * we, 2 * Te)
                if tc is not None:
                    axis.plot(tc, yc, c, lw=1.2, label=lbl)
        ax[5][0].set_ylabel("running mean (2 Te): F1" + unit, color="C0")
        tw.set_ylabel("F2, F6", color="C1")
        ax[5][0].grid(alpha=0.3)
        ax[5][0].axhline(0, color="k", lw=0.5, alpha=0.4)
        h1, l1 = ax[5][0].get_legend_handles_labels()
        h2, l2 = tw.get_legend_handles_labels()
        ax[5][0].legend(h1 + h2, l1 + l2, fontsize=8, ncol=3,
                        loc="lower right", framealpha=0.9)
        mark(ax[5][0], t_ramp, lo, hi)
    else:
        for i in range(6):
            ax[i][0].text(0.5, 0.5, "no drift-load output", ha="center",
                          transform=ax[i][0].transAxes)

    # --- right column: body motions --------------------------------------
    dofs = [("surge", 1, 1.0, "m", A, "A"), ("sway", 2, 1.0, "m", A, "A"),
            ("heave", 3, 1.0, "m", A, "A"),
            ("roll", 4, np.degrees(1), "deg", k * A, "kA"),
            ("pitch", 5, np.degrees(1), "deg", k * A, "kA"),
            ("yaw", 6, np.degrees(1), "deg", k * A, "kA")]

    if mot is not None and mot.shape[1] >= 7:
        t = mot[:, 0]
        m = (t >= lo) & (t <= hi)
        still = np.allclose(mot[:, 1:7], 0)
        for i, (name, col, scale, u, norm, nlbl) in enumerate(dofs):
            y = mot[:, col]
            panel(ax[i][1], t, y * scale, f"{name} [{u}]", colour="C1")
            mark(ax[i][1], t_ramp, lo, hi)
            if not still and m.sum() >= 20:
                mean, amp, ph = fit(t[m], y[m], we, a.nharm)
                level(ax[i][1], (mean + amp) * scale, "{:+.4g}")
                level(ax[i][1], (mean - amp) * scale, "{:+.4g}")
                ax[i][1].set_title(
                    f"{name}:  amp {amp*scale:.4g} {u}   "
                    f"RAO {amp/norm:.4f} ({nlbl})   phase "
                    f"{np.degrees(ph):+.1f} deg",
                    fontsize=9, loc="left")
        if still:
            ax[0][1].set_title("body restrained (fixBody true)", fontsize=9,
                               loc="left", color="C3")
    else:
        for i in range(6):
            ax[i][1].text(0.5, 0.5, "no motion output", ha="center",
                          transform=ax[i][1].transAxes)

    for j in range(2):
        ax[5][j].set_xlabel("t [s]")
    if a.tmin is not None or a.tmax is not None:
        ax[0][0].set_xlim(a.tmin, a.tmax)

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
