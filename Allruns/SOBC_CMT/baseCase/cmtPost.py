#!/usr/bin/env python3
"""
Post-processing of a circular motion test (manFlowPrescribed).

EVERYTHING REPORTED IS IN MMG AXES (x forward, y starboard, z down, origin at
midship = CofR), the convention of manModel:

    heading  encounter angle mu [deg] = (waveDirection + 180 - psi) mod 360,
             0 = head sea, 90 = waves from starboard, 180 = following,
             270 = waves from port (manModel's waveLoads encounter angle;
             waveDirection = propagation direction and psi = heading, both
             clockwise from north).  A starboard turn (r > 0) from head seas
             runs mu = 0, 350, 340, ...: the waves come from port first.
    U, V     surge and sway speed at midship, V > 0 to starboard (= u, v of
             constant/prescribedMotion)
    F1 = X   > 0 forward          (added resistance is X < 0)
    F2 = Y   > 0 to starboard
    Mz = N   about z down through midship, > 0 turns the bow to starboard
    F / (rho g A^2 B^2 / L),  N / (rho g A^2 B^2),  A the wave amplitude.

The function objects write in mesh axes (x aft, y starboard, z up, bow at -x;
the hull never rotates in the mesh), so

    X = -F_x,mesh,   Y = F_y,mesh,   N = -M_z,mesh.

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

Loads (all non-dimensional, MMG axes):
    rot      meanLoadsRot 'total': midfield for the rotating control volume
             (Chen + Coriolis + storage + centripetal).  THE result.
    rot2     the same on the larger control volume (must agree with rot)
    chen     Chen's midfield alone (= middleFieldForm); wrong while turning
    cor, sto the Coriolis and storage parts of rot (they include the
             waterline corner through P)
    cen      the centripetal part of rot (waterline corner only; written by
             middleFieldFormRot with cornerTerm on, absent in older runs)
    near     near-field (not reliable, for comparison only)
    table    the straight-course table at the same encounter angle

Plot: folded onto mu in [0, 180] (waves from starboard).  Waves from port
(mu > 180) are plotted at 360 - mu with Y and N sign-flipped (the mirror
image, exact for the port/starboard-symmetric hull), in a different colour.
cmt_meanLoads.csv keeps the unfolded values and adds mu_fold and mirror.

Per-heading results (results.csv, the manModel waveData format):

    lam/L,U,V,r,heading,F1mean,F2mean,Mzmean,eta1,eta2,eta3,eta4,eta5,eta6
    # convention: MMG ...                       (marker line, see MMG_MARKER)

U, V, r are the prescribed surge and sway speed at midship [m/s] and the
steady yaw rate [rad/s] (> 0 turning to starboard) of the test, model scale;
the rows within the yaw-rate ramp right after the onset are transitional.
One row per dh (default 10 deg) of encounter angle along the turn, from the
yaw onset (head seas, always there) to a full turn if the run gets that far.
Head seas at the onset and after a full turn are written as 360 and 0 for a
starboard turn (mu decreasing) and as 0 and 360 for a port turn, so that each
keeps its neighbours in the table.  F1..Mz are the filtered
loads of the rotating midfield (meanLoadsRot; --resultsCV rot2 for the large
control volume) at the time the heading is reached; eta1..6 the first-order
motion amplitudes there (harmonic fit at the local midship encounter
frequency over +-1 encounter period), eta1..3 / A, eta4..6 / (k A).

--table: a straight-course table for comparison.  An MMG table (with the
marker line) is used as it is; an older table without it is taken to be in
the old OpenFOAM convention (heading h = 180 - (chi - psi), 90 = waves from
port; F1 > 0 added resistance; F2 > 0 starboard; Mz about z up) and converted.

usage:
    python3 cmtPost.py
    python3 cmtPost.py --table /path/to/waveData.dat --U 0.33
    python3 cmtPost.py --dh 5 --csv ../results_CMT.csv   # also merge into a collection
Writes cmt_meanLoads.csv, cmt_meanLoads.png and results.csv.  By default the
output starts at the yaw onset; the function objects run from t = 0, so the
centred filter has the straight run before the onset to work on.  The
yaw-rate ramp is shaded in the plot: there dOmega/dt != 0 (Euler force,
unsteady basis flow), which the midfield formula, derived for constant Omega,
does not contain.
"""

import argparse
import glob
import os
import re
import shutil
import numpy as np

RESULTS_COLS = (["lam/L", "U", "V", "r", "heading", "F1mean", "F2mean", "Mzmean"]
                + [f"eta{i}" for i in range(1, 7)])
N_KEYS = 5                                   # lam/L, U, V, r, heading identify a row

MMG_MARKER = ("# convention: MMG body axes (x forward, y starboard, z down, midship);"
              " heading = encounter angle mu = (waveDirection + 180 - psi) mod 360 [deg],"
              " 0 head sea, 90 waves from starboard;"
              " F1mean = X > 0 forward, F2mean = Y > 0 starboard,"
              " Mzmean = N > 0 bow to starboard; U, V surge and sway at midship (V > 0 starboard)"
              " [m/s]; r yaw rate [rad/s] (> 0 turning to starboard)")


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


# ------------------------------------------------- tables and conventions ----
def is_mmg(path):
    """True if the table carries the MMG marker line."""
    with open(path) as f:
        for line in f:
            if line.startswith("#") and "convention: MMG" in line:
                return True
    return False


def old_to_mmg(head, data):
    """Rows of a table in the old OpenFOAM convention (heading h, 90 = waves
    from port; F1 > 0 aft; F2 > 0 starboard; Mz about z up) -> MMG, physical
    (no mirroring): mu = 360 - h, X = -F1, Y = F2, N = -Mz.  No modulo: in a
    circular-motion table h = 0 (yaw onset) and h = 360 (end of the turn)
    become 360 and 0 and keep their neighbours."""
    ih = head.index("heading")                # loads in the 3 columns after it
    d = data.copy()
    h = d[:, ih]
    d[:, ih] = (360.0 - h) if np.all((h >= 0.0) & (h <= 360.0)) else (360.0 - h) % 360.0
    d[:, ih + 1] *= -1.0                      # X = -F1
    d[:, ih + 3] *= -1.0                      # N = -Mz; Y = F2
    return d


def canonical_half(head, data):
    """A table covering only waves from port (no mu in (0, 180)) is mirrored
    onto waves from starboard: mu -> 360 - mu, V -> -V, r -> -r, Y -> -Y,
    N -> -N."""
    ih = head.index("heading")
    mu = data[:, ih] % 360.0
    if np.any((mu > 1e-6) & (mu < 180.0 - 1e-6)) or not np.any(mu > 180.0 + 1e-6):
        return data
    d = data.copy()
    d[:, ih] = (360.0 - mu) % 360.0
    for j in [ih + 2, ih + 3] + [head.index(n) for n in ("V", "r") if n in head[:ih]]:
        d[:, j] *= -1.0
    return d


def write_results(path, rows, merge=False):
    """Write rows (dicts with RESULTS_COLS) to path in the MMG convention.
    merge: keep the other rows of an existing file (an existing file without
    the MMG marker is converted from the old convention first, after a backup
    to <path>.oldConvention; one without the yaw-rate column r, whose rows'
    yaw rate is unknown, is not merged but moved to <path>.noYawRate)."""
    def key(vals):
        return tuple(round(float(x), 6) for x in vals[:N_KEYS])

    table = {}
    if merge and os.path.isfile(path):
        with open(path) as f:
            head = f.readline().strip().split(",")
        if head != RESULTS_COLS and "r" not in head:
            shutil.move(path, path + ".noYawRate")
            print(f"  {path}: no yaw-rate column r, not merged (moved to {path}.noYawRate)")
        elif head == RESULTS_COLS:
            data = np.atleast_2d(np.genfromtxt(path, delimiter=",", skip_header=1))
            if data.size:
                if not is_mmg(path):
                    shutil.copy2(path, path + ".oldConvention")
                    data = old_to_mmg(head, data)
                    print(f"  {path}: old convention, converted to MMG"
                          f" (backup {path}.oldConvention)")
                for r_ in data:
                    vals = [f"{x:.6g}" for x in r_]
                    table[key(vals)] = vals
    for row in rows:
        vals = [f"{row[c]:.6g}" for c in RESULTS_COLS]
        table[key(vals)] = vals
    with open(path, "w") as f:
        f.write(",".join(RESULTS_COLS) + "\n")
        f.write(MMG_MARKER + "\n")
        for kk in sorted(table):
            f.write(",".join(table[kk]) + "\n")


def table_loads(path, U, mu):
    """X, Y, N (MMG) from a straight-course table at speed U (the V and r
    closest to 0) and encounter angle mu [deg, 0..360].  A half table (waves
    from starboard, mu in [0, 180]) is mirrored for waves from port."""
    with open(path) as f:
        head = [s.strip() for s in f.readline().split(",")]
    data = np.atleast_2d(np.genfromtxt(path, delimiter=",", skip_header=1))[:, :len(head)]
    if not is_mmg(path):
        print(f"  table {path}: no MMG marker, read as the old OpenFOAM convention"
              " and converted")
        data = old_to_mmg(head, data)
    data = canonical_half(head, data)
    ih = head.index("heading")                # lam/L, U, [V], [r], heading, X, Y, N
    sel = np.ones(len(data), bool)
    for name in ("V", "r"):
        if name in head[:ih]:
            j = head.index(name)
            vals = np.unique(data[sel, j])
            sel &= np.isclose(data[:, j], vals[np.argmin(np.abs(vals))])
    speeds = np.unique(data[sel, 1])
    Us = speeds[np.argmin(np.abs(speeds - U))]
    sel &= np.isclose(data[:, 1], Us)
    d = data[sel]
    hh = d[:, ih]
    o = np.argsort(hh)
    hh, d = hh[o], d[o]
    mu = np.asarray(mu) % 360.0
    full = np.any(hh > 180.0 + 1e-6)
    if full:                                  # periodic: close the circle
        if hh[0] < 1e-6 and hh[-1] < 360.0 - 1e-6:
            hh, d = np.append(hh, 360.0), np.vstack([d, d[:1]])
        elif hh[0] > 1e-6 and hh[-1] > 360.0 - 1e-6:
            hh, d = np.insert(hh, 0, 0.0), np.vstack([d[-1:], d])
        q, sgn = mu, np.ones_like(mu)
    else:                                     # half table: mirror waves from port
        q = np.where(mu > 180.0, 360.0 - mu, mu)
        sgn = np.where(mu > 180.0, -1.0, 1.0)
    X = np.interp(q, hh, d[:, ih + 1])
    Y = sgn*np.interp(q, hh, d[:, ih + 2])
    N = sgn*np.interp(q, hh, d[:, ih + 3])
    return Us, X, Y, N


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
                    help="encounter-angle step for results.csv [deg]")
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
    # in mesh axes, a = chi - psi
    al = np.radians(chi - psi)
    eV0 = np.cos(al)*u + np.sin(al)*v
    we = np.abs(-w0 + k*eV0)
    Te = 2*np.pi/we

    # MMG encounter angle, and its fold onto waves from starboard [0, 180]
    mu = (180.0 + chi - psi) % 360.0
    mu_fold = np.where(mu > 180.0, 360.0 - mu, mu)
    mirror = np.where(mu > 180.0, -1.0, 1.0)            # -1: waves from port

    tskip = a.tskip
    if tskip is None:
        if r:
            tskip = tOn
        else:
            tskip = (num(wc, "rampPeriods", 3.0) + 3.0)*2*np.pi/w0

    out = {"t": t, "psi": psi, "mu": mu, "mu_fold": mu_fold, "mirror": mirror, "Te": Te}

    sigma = a.sigma*np.max(Te[t >= min(tskip, t[-1])])
    if a.smoother == "gauss":
        smooth = lambda y: gauss(t, y, sigma)
        print(f"  smoother: Gaussian, sigma = {sigma:.3f} s")
    else:
        smooth = lambda y: boxcar(t, y, Te)
        print("  smoother: boxcar over one midship encounter period")

    def add(name, F, M, cF, cM):
        """F, M loaded arrays (mesh axes); cF: first column of the force
        vector, cM: column of the moment z component.  Stored in MMG axes:
        X = -F_x, Y = F_y, N = -M_z."""
        out[f"F1_{name}"] = -smooth(np.interp(t, F[:, 0], F[:, cF]))/denF
        out[f"F2_{name}"] = smooth(np.interp(t, F[:, 0], F[:, cF + 1]))/denF
        out[f"Mz_{name}"] = -smooth(np.interp(t, M[:, 0], M[:, cM]))/denM

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
        Us, T1, T2, T6 = table_loads(a.table, Uq, mu)
        out["F1_table"], out["F2_table"], out["Mz_table"] = T1, T2, T6
        print(f"  table {a.table} at U = {Us}")

    keep = t >= tskip
    cols = list(out)
    np.savetxt(f"cmt_meanLoads{a.tag}.csv", np.column_stack([out[c][keep] for c in cols]),
               delimiter=",", header=",".join(cols), comments="")
    print(f"  u {u}  v {v}  r {r} (onset {tOn:.2f} s, ramp {tRamp:.2f} s)"
          f"   lam {lam}  A {A:.4g}   denF {denF:.4g} N  denM {denM:.4g} Nm")
    if keep.any():
        print(f"  analysed t >= {tskip:.2f} s: heading psi {psi[keep][0]:.1f} -> {psi[keep][-1]:.1f} deg,"
              f" encounter angle mu {mu[keep][0]:.1f} -> {mu[keep][-1]:.1f} deg")
    if "F1_rot" in out and "F1_rot2" in out and keep.any():
        for q in ("F1", "F2", "Mz"):
            d = out[f"{q}_rot"][keep] - out[f"{q}_rot2"][keep]
            print(f"  control-volume dependence {q}: max |rot - rot2| = {np.nanmax(np.abs(d)):.3f}")
    print(f"  wrote cmt_meanLoads{a.tag}.csv  (MMG axes)")

    # -- results.csv: one row per dh of encounter angle along the turn -------
    cv = a.resultsCV
    if f"F1_{cv}" in out:
        mu_u = np.degrees(np.unwrap(np.radians(180.0 + chi - psi)))
        mot = load("postProcessing/bodyMotion/motion.dat")
        ok = keep & np.isfinite(out[f"F1_{cv}"])
        rows = []
        if r and ok.any():
            tt, mm = t[ok], mu_u[ok]
            sgn = 1.0 if mm[-1] >= mm[0] else -1.0           # mu decreases for r > 0
            hh = np.maximum.accumulate(sgn*mm)               # increasing along the turn
            h0, h1 = hh[0], hh[-1]
            first = np.ceil(h0/a.dh - 1e-6)*a.dh
            targets = np.arange(first, h1 + 1e-9, a.dh)
            # t at which each angle is reached (first sample at or beyond it)
            tk = np.array([tt[np.searchsorted(hh, ht - 1e-9)] for ht in targets])
            for ht, tc in zip(targets, tk):
                mrep = (sgn*ht) % 360.0
                if mrep < 1e-6 or mrep > 360.0 - 1e-6:
                    # head seas at the onset and after a full turn: label them
                    # so that each keeps its neighbours in the table.  mu
                    # decreasing (starboard turn): onset 360, end 0; mu
                    # increasing (port turn): onset 0, end 360
                    end = abs(ht - h0) > 1.0
                    mrep = 360.0 if (end == (sgn > 0)) else 0.0
                rows.append((mrep, tc))
        elif ok.any():
            tc = t[ok][-1]
            rows.append((mu[ok][-1], tc))

        res = []
        for mrep, tc in rows:
            i = int(np.argmin(np.abs(t - tc)))
            row = {"lam/L": lam/L, "U": u, "V": v, "r": r, "heading": mrep,
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
            print(f"  wrote results.csv (MMG): {len(res)} encounter angles"
                  f" ({res[0]['heading']:.0f} .. {res[-1]['heading']:.0f} deg, step {a.dh:g},"
                  f" control volume {cv})")
            if a.csv:
                write_results(a.csv, res, merge=True)
                print(f"  merged into {a.csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    # name, line style, colour for waves from starboard, colour for waves from
    # port (mirrored), label
    styles = (("rot", "-", "tab:blue", "tab:orange", "midfield, rotating CV"),
              ("rot2", "--", "tab:cyan", "tab:brown", "same, large CV"),
              ("chen", ":", "tab:green", "tab:olive", "Chen only (uncorrected)"),
              ("table", "-", "k", "k", "straight-course table"),
              ("near", ":", "tab:red", "tab:pink", "near-field"),
              ("cen", ":", "tab:purple", "tab:gray", "centripetal part (waterline corner)"))

    # contiguous pieces with the same side, so the folded curve runs without
    # a jump where the waves change side
    idx = np.flatnonzero(keep)
    pieces = []
    if idx.size:
        breaks = np.flatnonzero(np.diff(mirror[idx]) != 0) + 1
        pieces = np.split(idx, breaks)

    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for ax, q, lab in zip(axs, ("F1", "F2", "Mz"),
                          (r"$X/(\rho g A^2 B^2/L)$  (F1, > 0 forward)",
                           r"$Y/(\rho g A^2 B^2/L)$  (F2, > 0 starboard)",
                           r"$N/(\rho g A^2 B^2)$  (Mz, > 0 bow to stbd)")):
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
                    col, lbl = c1, (f"{leg}, waves from stbd" if "a" not in labelled else None)
                    labelled.add("a")
                else:
                    col, lbl = c2, (f"{leg}, waves from port (mirrored)" if "b" not in labelled else None)
                    labelled.add("b")
                ax.plot(mu_fold[pc], y, ls, color=col, label=lbl, lw=1.2)
        if r and tRamp > 0:
            hr = np.interp([tOn, tOn + tRamp], t, mu_fold)
            ax.axvspan(min(hr), max(hr), color="0.85", zorder=0,
                       label="yaw-rate ramp" if q == "F1" else None)
        ax.set_ylabel(lab, fontsize=9)
        ax.set_xlim(0, 180)
        ax.grid(True)
    axs[0].legend(fontsize=7)
    axs[-1].set_xlabel("encounter angle mu [deg] (MMG: 0 head sea, 90 waves from starboard);"
                       "\nwaves from port folded to 360 - mu with Y and N sign-flipped",
                       fontsize=9)
    axs[0].set_title(f"CMT  u={u} v={v} r={r}   lam={lam} m   (MMG axes)")
    fig.tight_layout()
    fig.savefig(f"cmt_meanLoads{a.tag}.png", dpi=150)
    print(f"  wrote cmt_meanLoads{a.tag}.png")


if __name__ == "__main__":
    main()
