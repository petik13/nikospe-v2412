#!/usr/bin/env python3
"""
First-order wave excitation force on the hemisphere, non-dimensionalised as in

    Zhao, Faltinsen, Krokstad & Aanesland (BOSS),
    "Wave-Current Interaction Effects on Large-Volume Structures"

    Fig. 7 : horizontal excitation   F1 / (rho g zeta_a pi R^2)
    Fig. 8 : vertical   excitation   F3 / (rho g zeta_a pi R^2)
    Fig. 9 : horizontal drift        F1_mean / (0.5 rho g zeta_a^2 D)
    abscissa (all figures)           omega_e^2 R / g

Note on the abscissa: the paper's omega is the ENCOUNTER frequency
(their Eq.5, omega = omega_0 + k U cos(beta-alpha); Eq.6, k = omega_0^2/g),
so omega_e is used here, matching meanVal.py.

The amplitude is taken from a least-squares Fourier fit at omega_e over a whole
number of encounter periods, which removes the mean and the 2*omega_e drift
content cleanly (a peak-picker does not).

usage:
    python3 excForce.py                 # default: forces1, last 6 encounter periods
    python3 excForce.py --nper 10
    python3 excForce.py --obj forces1 --nharm 3
"""

import argparse
import io
import glob
import os
import re
import sys

import numpy as np

RHO = 1000.0
G = 9.81


# ----------------------------------------------------------------- input ----
def read_dict(path):
    """Crude OpenFOAM dictionary scalar reader: 'key  value;'."""
    vals = {}
    if not os.path.isfile(path):
        sys.exit(f"ERROR: {path} not found (run from the case directory)")
    for line in open(path, errors="ignore"):
        line = line.split("//")[0].strip()
        m = re.match(r"^([A-Za-z_]\w*)\s+(-?[\d.eE+-]+)\s*;", line)
        if m:
            try:
                vals[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
    return vals


def force_object_rhoref(obj):
    """rhoRef the forces FO actually applied to force.dat.

    This differs by FO type, and NOT the way the dictionary suggests:

      myFunctionObject : rhoRef is HARD-CODED to 1 in myFunctionObject.C --
          the line reads  rhoRef_ = 1;// dict.getCheck<scalar>("rhoInf",...)
          so an 'rhoInf' entry in the dict is silently ignored and the solver
          log prints 'rho Ref :1'.  With a kinematic p (m^2/s^2) the output is
          therefore force-per-unit-density and must be multiplied by rho.

      meanWaveLoads(Claude) : reads rhoInf properly, so its force.dat is
          already in Newtons.

    Returns (rhoRef, explanation).
    """
    path = "system/controlDict"
    if not os.path.isfile(path):
        return 1.0, "controlDict not found, assuming rhoRef = 1"
    txt = open(path, errors="ignore").read()
    m = re.search(r"^\s*" + re.escape(obj) + r"\s*$", txt, re.M)
    if not m:
        return 1.0, f"'{obj}' block not found in controlDict, assuming rhoRef = 1"
    blk = txt[m.end():]
    end = blk.find("\n\t}")
    if end < 0:
        end = blk.find("\n    }")
    blk = blk[:end if end > 0 else 2000]

    ftype = re.search(r"\btype\s+(\w+)\s*;", blk)
    ftype = ftype.group(1) if ftype else "?"

    if ftype in ("myFunctionObject", "myMeanForce"):   # legacy: rhoRef hard-coded to 1
        return 1.0, (f"type {ftype}: rhoRef_ is hard-coded to 1 in the source, "
                     "any rhoInf entry is ignored")

    r = re.search(r"\brhoInf\s+([-\d.eE+]+)\s*;", blk)
    if r:
        return float(r.group(1)), f"type {ftype}: rhoInf read from controlDict"
    return 1.0, f"type {ftype}: no rhoInf entry -> rhoRef = 1"


def load_forces(obj):
    """Concatenate every postProcessing/<obj>/*/force.dat, sorted, de-duplicated."""
    files = sorted(glob.glob(f"postProcessing/{obj}/*/force.dat"))
    if not files:
        sys.exit(f"ERROR: no postProcessing/{obj}/*/force.dat found")
    rows = []
    for f in files:
        raw = open(f, errors="ignore").read().replace("(", " ").replace(")", " ")
        d = np.loadtxt(io.StringIO(raw), comments="#")
        if d.size:
            rows.append(np.atleast_2d(d))
    d = np.vstack(rows)
    d = d[np.argsort(d[:, 0])]
    _, keep = np.unique(d[:, 0], return_index=True)   # last write wins on restart
    return d[keep], files


# ------------------------------------------------------------- harmonics ----
def harmonic_fit(t, y, w, nharm):
    """Least-squares fit  y = a0 + sum_n [ an cos(n w t) + bn sin(n w t) ].

    Returns (mean, [(amp_n, phase_n), ...]) with y_n = amp_n * cos(n w t + phase_n).
    """
    M = [np.ones_like(t)]
    for n in range(1, nharm + 1):
        M += [np.cos(n * w * t), np.sin(n * w * t)]
    M = np.column_stack(M)
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    out = []
    for n in range(1, nharm + 1):
        a, b = c[2 * n - 1], c[2 * n]
        out.append((np.hypot(a, b), np.arctan2(-b, a)))
    return c[0], out


def window(t, nper, Te, t_ramp):
    """Mask for the last whole number of encounter periods, after the ramp."""
    t_end = t[-1]
    t_beg = t_end - nper * Te
    if t_beg < t_ramp:                       # not enough clean signal: use what there is
        n_avail = int(np.floor((t_end - t_ramp) / Te))
        if n_avail < 1:
            return None, 0
        t_beg = t_end - n_avail * Te
        nper = n_avail
    return (t >= t_beg) & (t <= t_end), nper


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obj", default="forces1", help="postProcessing object name")
    ap.add_argument("--nper", type=int, default=6, help="encounter periods to fit over")
    ap.add_argument("--nharm", type=int, default=3, help="harmonics in the fit")
    ap.add_argument("--fk", default=None, metavar="OBJ",
                    help="FO holding the incident/Froude-Krylov force (e.g. forcesFK). "
                         "The diffraction part is then the COMPLEX difference "
                         "F_D = F_total - F_FK, which is what must be compared "
                         "against theory -- subtracting magnitudes is wrong.")
    a = ap.parse_args()

    wc = read_dict("constant/waveConditions")
    bd = read_dict("constant/bodyMotionProperties")
    lam = wc["waveLength"]
    steep = wc["steepness"]
    U0 = wc["currentSpeed"]
    head = wc["headingAngle"]                     # radians
    h = wc["waterDepth"]
    ramp_per = wc.get("rampPeriods", 0.0)
    R = bd.get("radius", 0.5)
    D = 2.0 * R

    k = 2.0 * np.pi / lam
    w0 = np.sqrt(G * k * np.tanh(k * h))      # absolute frequency
    we = w0 + k * U0 * np.cos(head)           # encounter frequency  (paper Eq.5)
    Te = 2.0 * np.pi / we
    zeta_a = 0.5 * steep * lam
    t_ramp = ramp_per * 2.0 * np.pi / w0
    x_paper = we**2 * R / G

    exc_den = RHO * G * zeta_a * np.pi * R**2      # Figs 7, 8
    drift_den = 0.5 * RHO * G * zeta_a**2 * D      # Fig 9

    d, files = load_forces(a.obj)
    t = d[:, 0]

    # Undo the FO's own density scaling and apply the real one.
    rho_ref, why = force_object_rhoref(a.obj)
    scale = RHO / rho_ref
    F1, F3 = d[:, 1] * scale, d[:, 3] * scale

    print(f"case            : {os.path.basename(os.getcwd())}")
    print(f"files           : {', '.join(files)}")
    print(f"lambda {lam:.4g} m   steepness {steep:.5g}   zeta_a {zeta_a:.5g} m"
          f"   R {R:.4g} m   h {h:.4g} m")
    print(f"U0 {U0:+.4g} m/s   head_ang {head:+.4g} rad"
          f"   tau = we*U/g = {we*U0/G:+.4f}")
    print(f"w0 {w0:.5g}   we {we:.5g} rad/s   Te {Te:.5g} s   ramp ends {t_ramp:.4g} s")
    print(f"t: {t[0]:.4g} -> {t[-1]:.5g} s  ({(t[-1]-t_ramp)/Te:.2f} encounter periods "
          f"after ramp)")
    print(f"rhoRef used by '{a.obj}' = {rho_ref:g}  ({why});"
          f"  forces scaled by rho/rhoRef = {scale:g}")
    print()
    print(f"PAPER ABSCISSA   we^2 R / g = {x_paper:.4f}")
    print()

    m, nper = window(t, a.nper, Te, t_ramp)
    if m is None or m.sum() < 20:
        sys.exit("ERROR: not enough post-ramp signal yet to fit a period. "
                 "Let the run go further.")

    print(f"fit: last {nper} encounter periods, {m.sum()} samples, "
          f"{a.nharm} harmonics")
    print()
    for name, F, in (("F1 (surge)", F1), ("F3 (heave)", F3)):
        mean, hs = harmonic_fit(t[m], F[m], we, a.nharm)
        amp, ph = hs[0]
        amp2 = hs[1][0] if a.nharm > 1 else 0.0
        print(f"  {name}")
        print(f"     amplitude at we      = {amp:.6g} N      phase = "
              f"{np.degrees(ph):+7.2f} deg")
        print(f"     NON-DIM (paper)      = {amp/exc_den:.4f}"
              f"        [ F/(rho g zeta_a pi R^2) ]")
        print(f"     mean                 = {mean:+.6g} N")
        print(f"     2nd harmonic / 1st   = {amp2/amp if amp else 0:.3f}"
              f"   (linearity check, should be small)")
        print()

    mean1, _ = harmonic_fit(t[m], F1[m], we, a.nharm)
    print(f"  drift cross-check (linear pressure only, NOT the CV value)")
    print(f"     mean F1 / (0.5 rho g zeta_a^2 D) = {mean1/drift_den:+.4f}")
    print()

    # --- diffraction part = complex(total) - complex(FK) ---------------------
    if a.fk:
        dfk, _ = load_forces(a.fk)
        rr, why_fk = force_object_rhoref(a.fk)
        sfk = RHO / rr
        tk = dfk[:, 0]
        mk, _ = window(tk, nper, Te, t_ramp)
        if mk is None or mk.sum() < 20:
            print(f"  '{a.fk}': not enough post-ramp signal to fit -- skipping")
            return
        print(f"  DIFFRACTION PART  = complex(total) - complex('{a.fk}')")
        print(f"     ('{a.fk}' rhoRef {rr:g}: {why_fk})")
        for name, col in (("F1 (surge)", 1), ("F3 (heave)", 3)):
            _, ht = harmonic_fit(t[m], (F1 if col == 1 else F3)[m], we, a.nharm)
            _, hk = harmonic_fit(tk[mk], dfk[mk, col] * sfk, we, a.nharm)
            At, pt = ht[0]
            Ak, pk = hk[0]
            zt = At * np.exp(1j * pt)
            zk = Ak * np.exp(1j * pk)
            zd = zt - zk
            print(f"     {name}")
            print(f"        total  |F| = {At/exc_den:.4f}  phase {np.degrees(pt):+7.2f} deg")
            print(f"        FK     |F| = {Ak/exc_den:.4f}  phase {np.degrees(pk):+7.2f} deg")
            print(f"        DIFFR  |F| = {abs(zd)/exc_den:.4f}  phase "
                  f"{np.degrees(np.angle(zd)):+7.2f} deg   <-- the diffraction force")
        print()

    # --- convergence: same fit over successive windows -----------------------
    print("  convergence of |F1|, |F3| over successive 2-Te windows:")
    n_tot = int(np.floor((t[-1] - t_ramp) / Te))
    for j in range(2, n_tot + 1):
        lo, hi = t[-1] - j * Te, t[-1] - (j - 2) * Te
        mm = (t >= lo) & (t <= hi)
        if mm.sum() < 20:
            continue
        _, h1 = harmonic_fit(t[mm], F1[mm], we, a.nharm)
        _, h3 = harmonic_fit(t[mm], F3[mm], we, a.nharm)
        print(f"     t=[{lo:7.3f},{hi:7.3f}]   |F1|nd = {h1[0][0]/exc_den:.4f}"
              f"   |F3|nd = {h3[0][0]/exc_den:.4f}")


if __name__ == "__main__":
    main()
