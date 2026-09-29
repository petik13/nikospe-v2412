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
             (Chen + Coriolis + storage).  THE result.
    rot2     the same on the larger control volume (must agree with rot)
    chen     Chen's midfield alone (= middleFieldForm); wrong while turning
    cor, sto the Coriolis and storage parts of rot
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

usage:
    python3 cmtPost.py
    python3 cmtPost.py --table /path/to/waveData.dat --U 0.33
Writes cmt_meanLoads.csv and cmt_meanLoads.png.  By default the analysis
starts two encounter periods after the end of the yaw-rate ramp.
"""

import argparse
import glob
import re
import numpy as np


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
                    help="discard t < tskip (default: end of the yaw ramp + 2 periods,"
                         " or wave ramp + 3 periods for r = 0)")
    ap.add_argument("--smoother", choices=("gauss", "boxcar"), default="gauss")
    ap.add_argument("--sigma", type=float, default=0.5,
                    help="Gaussian standard deviation in largest encounter periods")
    ap.add_argument("--tag", default="",
                    help="suffix for the output files, cmt_meanLoads<tag>.csv/.png")
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

    tskip = a.tskip
    if tskip is None:
        if r:
            tskip = tOn + tRamp + 2*Te[0]
        else:
            tskip = (num(wc, "rampPeriods", 3.0) + 3.0)*2*np.pi/w0

    out = {"t": t, "psi": psi, "h_table": hdg, "Te": Te}

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
        out[f"F1_{name}"] = smooth(np.interp(t, F[:, 0], F[:, cF]))/denF
        out[f"F2_{name}"] = smooth(np.interp(t, F[:, 0], F[:, cF + 1]))/denF
        out[f"Mz_{name}"] = smooth(np.interp(t, M[:, 0], M[:, cM]))/denM

    # Rotating-frame midfield: Time total chen surface elevation strip coriolis storage
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

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    styles = (("rot", "b-", "midfield, rotating CV"),
              ("rot2", "c--", "same, large CV"),
              ("chen", "g:", "Chen only (uncorrected)"),
              ("table", "k-", "straight-course table"),
              ("near", "r:", "near-field"))
    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for ax, q, lab in zip(axs, ("F1", "F2", "Mz"),
                          (r"$F_1/(\rho g A^2 B^2/L)$", r"$F_2/(\rho g A^2 B^2/L)$",
                           r"$M_z/(\rho g A^2 B^2)$")):
        for name, sty, leg in styles:
            key = f"{q}_{name}"
            if key in out:
                ax.plot(hdg[keep], out[key][keep], sty, label=leg, lw=1.2)
        ax.set_ylabel(lab)
        ax.grid(True)
    axs[0].legend(fontsize=8)
    axs[-1].set_xlabel("equivalent table heading h [deg] (0 = head sea)")
    axs[0].set_title(f"CMT  u={u} v={v} r={r}   lam={lam} m")
    fig.tight_layout()
    fig.savefig(f"cmt_meanLoads{a.tag}.png", dpi=150)
    print(f"  wrote cmt_meanLoads{a.tag}.png")


if __name__ == "__main__":
    main()
