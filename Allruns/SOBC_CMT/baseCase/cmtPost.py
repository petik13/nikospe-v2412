#!/usr/bin/env python3
"""
Post-processing of a circular motion test (manFlowPrescribed).

Mean (second-order) wave loads during the turn.  In the rotating frame the
local encounter frequency varies over the control volume,

    omega_e(x) = |-omega + k e.V_S(x)|,   V_S = V0 + Omega x r,

so the quadratic integrands oscillate at a BAND of frequencies 2 omega_e(x),
+-2 k Omega r_max wide (+-0.1 Hz for the small, +-0.24 Hz for the large
KVLCC2 control volume at r = 0.08).  Their amplitude is 5-20 times the mean
(storage and strip terms in particular).  A boxcar over one midship encounter
period only removes the centre frequency and leaves a ripple at 2 omega_e.

Default smoother: a centred Gaussian (zero phase), sigma = --sigma x the
largest encounter period (default 0.5).  Its response exp(-2 pi^2 sigma^2 f^2)
is ~1e-6 at the lowest 2 omega_e and >= 0.95 at the harmonics of the turn rate
(<= 0.07 Hz) that make up the physical mean.  --smoother boxcar gives the
one-period boxcar (the encounter frequency at midship is
|dtheta/dt| = |-omega + k e(t).V0|).

Loads (all non-dimensional, rho g A^2 B^2 / L and rho g A^2 B^2):
    rot      meanLoadsRot 'total': midfield for the rotating control volume
             (Chen + Coriolis + storage + centripetal).  THE result.
    rot2     the same on the larger control volume (must agree with rot)
    chen     Chen's midfield alone (= middleFieldForm); wrong while turning
    cor, sto the Coriolis and storage parts of rot (they include the
             waterline corner through P)
    cen      the centripetal part of rot (waterline corner only; written by
             middleFieldFormRot with cornerTerm on, absent in older runs)
    near     near-field (not reliable, for comparison only)
    table    the straight-course table at the same encounter heading

Axes: the hull is always aligned with the mesh (bow at -x, y starboard), so
the function-object loads are already in the axes of the waveData tables
(meanLoads.py rotates the fixed-heading runs into exactly these axes):
F1 > 0 is added resistance, and F2, Mz are as in the table columns
F1mean, F2mean, Mzmean.

Equivalent table heading h (0 = head sea, table convention) at time t:
    h = 180 - (chi - psi(t))      wrapped to (-180, 180]
The table is tabulated for h in [0, 180]; for h < 0 the mirror image is used,
F1(h) = F1(|h|), F2(h) = -F2(|h|), Mz(h) = -Mz(|h|).
In the plot, the second half of a turn (h > 180, i.e. h < 0 in (-180, 180])
is folded onto the table range: plotted at 360 - h = |h| (running 180 -> 0),
with F2 and Mz sign-flipped, in a different colour.  The CSV keeps the
unfolded values and adds h_fold = |h| and mirror = sign.

Per-heading results (results.csv, same format as the meanLoads.py --csv
collection, plus V): one row per heading 0, dh, 2 dh, ... (default dh = 10 deg)
along the turn, from the yaw onset (heading 0 is always there) to 360 deg if
the run gets that far:

    lam/L,U,V,heading,F1mean,F2mean,Mzmean,eta1,eta2,eta3,eta4,eta5,eta6

    lam/L      wavelength / Lpp
    U, V       prescribed surge and sway speed at midship [m/s] (MMG axes)
    heading    table heading h [deg] (0 = head sea), unwrapped along the turn
               and reported in [0, 360]; NOT folded, so F2 and Mz for
               h > 180 have the signs of the actual heading
    F1..Mz     filtered mean loads of the rotating midfield (meanLoadsRot,
               --resultsCV rot2 for the large control volume), at the time
               the heading is reached; F / (rho g A^2 B^2 / L), Mz / (rho g A^2 B^2)
    eta1..6    first-order motion amplitudes at that time (harmonic fit at the
               local midship encounter frequency over +-1 encounter period),
               eta1..3 / A, eta4..6 / (k A), as meanLoads.py

FORCE_SIGN (below) is the sign convention of this hull's meanLoads.py for F1
and F2 (KVLCC2 +1, SOBC -1), so the CMT loads and tables agree.

usage:
    python3 cmtPost.py
    python3 cmtPost.py --table /path/to/waveData.dat --U 0.33
    python3 cmtPost.py --dh 5 --csv ../results_CMT.csv   # also merge into a collection
Writes cmt_meanLoads.csv and cmt_meanLoads.png.  By default the output starts
at the yaw onset (heading change 0); the function objects run from t = 0, so
the centred filter has the straight run before the onset to work on.  The
yaw-rate ramp is shaded in the plot: there dOmega/dt != 0 (Euler force,
unsteady basis flow), which the midfield formula, derived for constant Omega,
does not contain.
"""

import argparse
import glob
import os
import re
import numpy as np

# Sign of F1, F2 in this hull's tables (the den_F of its meanLoads.py):
# KVLCC2 +1; SOBC -1 (SOBC/baseCase/meanLoads.py divides by -rho g A^2 B^2/L).
# Mz is not affected.
FORCE_SIGN = -1.0

RESULTS_COLS = (["lam/L", "U", "V", "heading", "F1mean", "F2mean", "Mzmean"]
                + [f"eta{i}" for i in range(1, 7)])


def read_dict(path):
    d = {}
    with open(path) as f:
        for line in f:
            s = line.split("//")[0].strip().rstrip(";")
            m = re.match(r"^(\w+)\s+(.+)$", s)
            if m:
                d[m.group(1)] = m.group(2).strip()
    return d


def num(d, k, default=None):
    return float(d[k]) if k in d else default


def load(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    rows = []
    for fn in files:
        with open(fn) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                vals = line.replace("(", " ").replace(")", " ").split()
                try:
                    rows.append([float(x) for x in vals])
                except ValueError:
                    pass
    if not rows:
        return None
    n = min(len(r) for r in rows)
    a = np.array([r[:n] for r in rows])
    _, idx = np.unique(a[:, 0], return_index=True)     # restarts
    return a[idx]


def boxcar(t, y, T):
    """Centred average of y over one period T(t) (T may vary slowly)."""
    c = np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*np.diff(t))])
    out = np.full_like(y, np.nan)
    lo = np.interp(t - 0.5*T, t, c, left=np.nan, right=np.nan)
    hi = np.interp(t + 0.5*T, t, c, left=np.nan, right=np.nan)
    ok = np.isfinite(lo) & np.isfinite(hi)
    out[ok] = (hi[ok] - lo[ok])/T[ok]
    return out


def gauss(t, y, sigma):
    """Centred Gaussian smoothing with standard deviation sigma [s], on a
    uniform resampling of t.  Points closer than 3 sigma to either end are
    returned as NaN (the kernel would be truncated there)."""
    dt = np.median(np.diff(t))
    tu = np.arange(t[0], t[-1] + 0.5*dt, dt)
    yu = np.interp(tu, t, y)
    n = int(np.ceil(4*sigma/dt))
    s = np.arange(-n, n + 1)*dt
    w = np.exp(-0.5*(s/sigma)**2)
    w /= w.sum()
    ys = np.convolve(yu, w, mode="same")
    ys[:n] = np.nan
    ys[-n:] = np.nan
    return np.interp(t, tu, ys, left=np.nan, right=np.nan)


def motion_amplitude(tm, ym, tc, wc, Tc):
    """First-harmonic amplitude of y around time tc: least-squares fit of
    a + b (t - tc) + c cos(wc t) + d sin(wc t) over tc +- Tc.  The linear term
    absorbs the slow drift of the soft-moored modes."""
    sel = (tm >= tc - Tc) & (tm <= tc + Tc)
    if sel.sum() < 8:
        return np.nan
    tt, yy = tm[sel], ym[sel]
    X = np.column_stack([np.ones_like(tt), tt - tc, np.cos(wc*tt), np.sin(wc*tt)])
    coef, *_ = np.linalg.lstsq(X, yy, rcond=None)
    return float(np.hypot(coef[2], coef[3]))


def write_results(path, rows):
    """Write rows (dicts with RESULTS_COLS) to path.  If path exists, rows with
    the same (lam/L, U, V, heading) are replaced and the rest kept."""
    def key(vals):
        return tuple(round(float(x), 6) for x in vals[:4])

    table = {}
    if os.path.isfile(path):
        with open(path) as f:
            head = f.readline().strip().split(",")
            if head == RESULTS_COLS:
                for line in f:
                    vals = line.strip().split(",")
                    if len(vals) == len(RESULTS_COLS):
                        table[key(vals)] = vals
    for row in rows:
        vals = [f"{row[c]:.6g}" for c in RESULTS_COLS]
        table[key(vals)] = vals
    with open(path, "w") as f:
        f.write(",".join(RESULTS_COLS) + "\n")
        for kk in sorted(table):
            f.write(",".join(table[kk]) + "\n")


def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0


def table_loads(path, U, h):
    """F1, F2, Mz from a waveData table at speed U (V = 0), heading |h|,
    mirrored for h < 0."""
    with open(path) as f:
        head = [s.strip() for s in f.readline().split(",")]
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    col = {name: i for i, name in enumerate(head)}
    sel = np.ones(len(data), bool)
    if "V" in col:
        sel &= np.isclose(data[:, col["V"]], 0.0)
    speeds = np.unique(data[sel, col["U"]])
    Us = speeds[np.argmin(np.abs(speeds - U))]
    sel &= np.isclose(data[:, col["U"]], Us)
    d = data[sel]
    d = d[np.argsort(d[:, col["heading"]])]
    hh = d[:, col["heading"]]
    ah = np.abs(h)
    sgn = np.where(h < 0, -1.0, 1.0)
    F1 = np.interp(ah, hh, d[:, col["F1mean"]])
    F2 = sgn*np.interp(ah, hh, d[:, col["F2mean"]])
    Mz = sgn*np.interp(ah, hh, d[:, col["Mzmean"]])
    return Us, F1, F2, Mz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=None, help="waveData table for comparison")
    ap.add_argument("--U", type=float, default=None, help="table speed (default: prescribed u)")
    ap.add_argument("--tskip", type=float, default=None,
                    help="discard t < tskip (default: the yaw onset,"
                         " or wave ramp + 3 periods for r = 0)")
    ap.add_argument("--smoother", choices=("gauss", "boxcar"), default="gauss")
    ap.add_argument("--sigma", type=float, default=0.5,
                    help="Gaussian standard deviation in largest encounter periods")
    ap.add_argument("--tag", default="",
                    help="suffix for the output files, cmt_meanLoads<tag>.csv/.png")
    ap.add_argument("--dh", type=float, default=10.0,
                    help="heading step for results.csv [deg]")
    ap.add_argument("--resultsCV", choices=("rot", "rot2"), default="rot",
                    help="control volume for results.csv")
    ap.add_argument("--csv", default=None, metavar="PATH",
                    help="also merge the results.csv rows into this collection file")
    a = ap.parse_args()

    wc = read_dict("constant/waveConditions")
    pm = read_dict("constant/prescribedMotion")
    bd = read_dict("constant/bodyMotionProperties")

    RHO, G = 1000.0, 9.81
    lam = num(wc, "waveLength")
    A = 0.5*num(wc, "steepness")*lam
    h = num(wc, "waterDepth")
    k = 2*np.pi/lam
    w0 = np.sqrt(G*k*np.tanh(k*h))
    L = num(bd, "Lpp", 1.0)
    B = num(bd, "beam", 1.0)
    u, v, r = num(pm, "u"), num(pm, "v"), num(pm, "r")
    chi = num(pm, "waveDirection")
    tOn = num(pm, "yawOnsetTime", 0.0)
    tRamp = num(pm, "yawRampTime", 0.0)
    denF = RHO*G*A**2*B**2/L
    denM = RHO*G*A**2*B**2

    tr = load("postProcessing/shipTrajectory/trajectory.dat")
    t, psi = tr[:, 0], tr[:, 1]

    # Local encounter frequency at midship, e = (-cos a, sin a), V0 = (-u, v)
    al = np.radians(chi - psi)
    eV0 = np.cos(al)*u + np.sin(al)*v
    we = np.abs(-w0 + k*eV0)
    Te = 2*np.pi/we

    hdg = wrap180(180.0 - (chi - psi))

    # Fold onto the table range [0, 180]: the hull is symmetric, so heading
    # h > 180 (hdg < 0) is the mirror image of 360 - h = |hdg|, with F2 and Mz
    # changing sign.  mirror = -1 marks those samples.
    hfold = np.abs(hdg)
    mirror = np.where(hdg < 0, -1.0, 1.0)

    tskip = a.tskip
    if tskip is None:
        if r:
            tskip = tOn
        else:
            tskip = (num(wc, "rampPeriods", 3.0) + 3.0)*2*np.pi/w0

    out = {"t": t, "psi": psi, "h_table": hdg, "h_fold": hfold, "mirror": mirror, "Te": Te}

    sigma = a.sigma*np.max(Te[t >= min(tskip, t[-1])])
    if a.smoother == "gauss":
        smooth = lambda y: gauss(t, y, sigma)
        print(f"  smoother: Gaussian, sigma = {sigma:.3f} s")
    else:
        smooth = lambda y: boxcar(t, y, Te)
        print("  smoother: boxcar over one midship encounter period")

    def add(name, F, M, cF, cM):
        """F, M loaded arrays; cF: first column of the force vector, cM: column
        of the moment z component"""
        out[f"F1_{name}"] = FORCE_SIGN*smooth(np.interp(t, F[:, 0], F[:, cF]))/denF
        out[f"F2_{name}"] = FORCE_SIGN*smooth(np.interp(t, F[:, 0], F[:, cF + 1]))/denF
        out[f"Mz_{name}"] = smooth(np.interp(t, M[:, 0], M[:, cM]))/denM

    # Rotating-frame midfield: Time total chen surface elevation strip coriolis
    # storage [centripetal]
    for fo, tag in (("meanLoadsRot", "rot"), ("meanLoadsRot2", "rot2")):
        F = load(f"postProcessing/{fo}/*/force.dat")
        M = load(f"postProcessing/{fo}/*/moment.dat")
        if F is None or M is None:
            continue
        add(tag, F, M, 1, 3)
        if tag == "rot":
            add("chen", F, M, 4, 6)
            add("cor", F, M, 16, 18)
            add("sto", F, M, 19, 21)
            if F.shape[1] >= 25 and M.shape[1] >= 25:
                add("cen", F, M, 22, 24)

    # Old objects, if present: first vector is the total
    if "F1_chen" not in out:
        F = load("postProcessing/meanLoads/*/force.dat")
        M = load("postProcessing/meanLoads/*/moment.dat")
        if F is not None and M is not None:
            add("chen", F, M, 1, 3)
    F = load("postProcessing/meanLoadsNear/*/force.dat")
    M = load("postProcessing/meanLoadsNear/*/moment.dat")
    if F is not None and M is not None:
        add("near", F, M, 1, 3)

    if a.table:
        Uq = u if a.U is None else a.U
        Us, T1, T2, T6 = table_loads(a.table, Uq, hdg)
        out["F1_table"], out["F2_table"], out["Mz_table"] = T1, T2, T6
        print(f"  table {a.table} at U = {Us}")

    keep = t >= tskip
    cols = list(out)
    np.savetxt(f"cmt_meanLoads{a.tag}.csv", np.column_stack([out[c][keep] for c in cols]),
               delimiter=",", header=",".join(cols), comments="")
    print(f"  u {u}  v {v}  r {r} (onset {tOn:.2f} s, ramp {tRamp:.2f} s)"
          f"   lam {lam}  A {A:.4g}   denF {denF:.4g} N  denM {denM:.4g} Nm")
    if keep.any():
        print(f"  analysed t >= {tskip:.2f} s: heading {psi[keep][0]:.1f} -> {psi[keep][-1]:.1f} deg,"
              f" table heading {hdg[keep][0]:.1f} -> {hdg[keep][-1]:.1f} deg")
    if "F1_rot" in out and "F1_rot2" in out and keep.any():
        for q in ("F1", "F2", "Mz"):
            d = out[f"{q}_rot"][keep] - out[f"{q}_rot2"][keep]
            print(f"  control-volume dependence {q}: max |rot - rot2| = {np.nanmax(np.abs(d)):.3f}")
    print(f"  wrote cmt_meanLoads{a.tag}.csv")

    # -- results.csv: one row per heading step along the turn ----------------
    cv = a.resultsCV
    if f"F1_{cv}" in out:
        hu = np.degrees(np.unwrap(np.radians(180.0 - (chi - psi))))
        mot = load("postProcessing/bodyMotion/motion.dat")
        ok = keep & np.isfinite(out[f"F1_{cv}"])
        rows = []
        if r and ok.any():
            sgn = 1.0 if r > 0 else -1.0
            tt, hh = t[ok], sgn*hu[ok]                      # hh increases along the turn
            hh = np.maximum.accumulate(hh)
            h0, h1 = hh[0], hh[-1]
            first = np.ceil(h0/a.dh - 1e-6)*a.dh
            targets = np.arange(first, h1 + 1e-9, a.dh)
            # t at which each heading is reached (first sample at or beyond it)
            tk = np.array([tt[np.searchsorted(hh, ht - 1e-9)] for ht in targets])
            for ht, tc in zip(targets, tk):
                hrep = (sgn*ht) % 360.0
                if hrep < 1e-6 and abs(ht - h0) > 1.0:
                    hrep = 360.0
                rows.append((hrep, tc))
        elif ok.any():
            tc = t[ok][-1]
            rows.append((hu[ok][-1] % 360.0, tc))

        res = []
        for hrep, tc in rows:
            i = int(np.argmin(np.abs(t - tc)))
            row = {"lam/L": lam/L, "U": u, "V": v, "heading": hrep,
                   "F1mean": out[f"F1_{cv}"][i], "F2mean": out[f"F2_{cv}"][i],
                   "Mzmean": out[f"Mz_{cv}"][i]}
            for j in range(6):
                if mot is not None:
                    amp = motion_amplitude(mot[:, 0], mot[:, 1 + j], tc, we[i], Te[i])
                else:
                    amp = np.nan
                row[f"eta{j + 1}"] = amp/(A if j < 3 else k*A)
            res.append(row)

        if res:
            write_results("results.csv", res)
            print(f"  wrote results.csv: {len(res)} headings"
                  f" ({res[0]['heading']:.0f} .. {res[-1]['heading']:.0f} deg, step {a.dh:g},"
                  f" control volume {cv})")
            if a.csv:
                write_results(a.csv, res)
                print(f"  merged into {a.csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    # name, line style, colour for h <= 180, colour for h > 180 (mirrored), label
    styles = (("rot", "-", "tab:blue", "tab:orange", "midfield, rotating CV"),
              ("rot2", "--", "tab:cyan", "tab:brown", "same, large CV"),
              ("chen", ":", "tab:green", "tab:olive", "Chen only (uncorrected)"),
              ("table", "-", "k", "k", "straight-course table"),
              ("near", ":", "tab:red", "tab:pink", "near-field"),
              ("cen", ":", "tab:purple", "tab:gray", "centripetal part (waterline corner)"))

    # contiguous pieces with the same mirror sign, so the folded curve runs
    # 0 -> 180 and then back 180 -> 0 without a jump
    idx = np.flatnonzero(keep)
    pieces = []
    if idx.size:
        breaks = np.flatnonzero(np.diff(mirror[idx]) != 0) + 1
        pieces = np.split(idx, breaks)

    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for ax, q, lab in zip(axs, ("F1", "F2", "Mz"),
                          (r"$F_1/(\rho g A^2 B^2/L)$", r"$F_2/(\rho g A^2 B^2/L)$",
                           r"$M_z/(\rho g A^2 B^2)$")):
        for name, ls, c1, c2, leg in styles:
            key = f"{q}_{name}"
            if key not in out:
                continue
            labelled = set()
            for pc in pieces:
                m = mirror[pc[0]]
                y = out[key][pc] if q == "F1" else m*out[key][pc]
                if name == "table":
                    col, lbl = c1, (leg if "table" not in labelled else None)
                    labelled.add("table")
                elif m > 0:
                    col, lbl = c1, (leg if "a" not in labelled else None)
                    labelled.add("a")
                else:
                    col, lbl = c2, (f"{leg}, h > 180 (mirrored)" if "b" not in labelled else None)
                    labelled.add("b")
                ax.plot(hfold[pc], y, ls, color=col, label=lbl, lw=1.2)
        if r and tRamp > 0:
            hr = np.interp([tOn, tOn + tRamp], t, hfold)
            ax.axvspan(min(hr), max(hr), color="0.85", zorder=0,
                       label="yaw-rate ramp" if q == "F1" else None)
        ax.set_ylabel(lab)
        ax.set_xlim(0, 180)
        ax.grid(True)
    axs[0].legend(fontsize=7)
    axs[-1].set_xlabel("table heading h [deg] (0 = head sea); h > 180 folded to 360 - h,"
                       " with F2 and Mz sign-flipped")
    axs[0].set_title(f"CMT  u={u} v={v} r={r}   lam={lam} m")
    fig.tight_layout()
    fig.savefig(f"cmt_meanLoads{a.tag}.png", dpi=150)
    print(f"  wrote cmt_meanLoads{a.tag}.png")


if __name__ == "__main__":
    main()
