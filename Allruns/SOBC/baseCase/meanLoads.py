#!/usr/bin/env python3
"""
Post-run summary: mean wave drift loads and motion RAOs, in the
non-dimensional form used by

    Seo, Ha, Nam & Kim (2021), "Experimental and Numerical Analysis of Wave
    Drift Force on SOBC Moving in Oblique Waves", JMSE 9, 136.

    drift force  (Figs 11, 14)   Fx, Fy / (rho g A^2 B^2 / L)
    drift moment (Fig 16)        Mz     / (rho g A^2 B^2)
    motion RAOs  (Figs 9, 10)    xi1..3 / A ,   xi4..6 / (k A)
    abscissa                     lambda / L

Note the two denominators differ by a factor L -- the moment is NOT divided
by L.

Means are taken from a least-squares fit at 2*omega_e (a drift load is
quadratic in a first-order field, so it is DC + 2w only) and RAO amplitudes
from a fit at omega_e.  Both beat averaging between zero crossings, which
biases as soon as the DC approaches the oscillation amplitude.

Axes: the function object integrates in mesh axes, where +x is the wave
propagation direction.  The drift loads are rotated through the heading into
ship axes before being printed, matching meanVal.py:

    F1 =  Fx cos(psi) - Fy sin(psi)        psi = -headingAngle
    F2 =  Fx sin(psi) + Fy cos(psi)

which puts the bow at -x', so F1 > 0 is added resistance -- directly
comparable with Seo et al.  Fz and Mz are about the rotation axis and are
left alone.  In head seas the rotation is 180 deg and only flips signs; at an
oblique heading it mixes x and y, so it matters.

The motion RAOs come from the solver's own body frame, which is this one
turned 180 deg about z.  Amplitudes are unaffected; the surge and sway phases
sit 180 deg from the load axes above.

usage:
    python3 meanLoads.py                  # summary for this case
    python3 meanLoads.py --nper 8
    python3 meanLoads.py --append ../sweep_head.dat
"""

import argparse
import glob
import io
import os
import re
import sys

import numpy as np

RHO = 1000.0
G = 9.81


# ----------------------------------------------------------------- input ----
def read_dict(path):
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


def load(pattern):
    """Load OpenFOAM tabular output, tolerating bracketed vectors."""
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    rows = []
    for f in files:
        raw = open(f, errors="ignore").read().replace("(", " ").replace(")", " ")
        d = np.atleast_2d(np.loadtxt(io.StringIO(raw), comments="#"))
        if d.size:
            rows.append(d)
    if not rows:
        return None
    d = np.vstack(rows)
    d = d[np.argsort(d[:, 0])]
    _, keep = np.unique(d[:, 0], return_index=True)
    return d[keep]


# ------------------------------------------------------------- harmonics ----
def fit(t, y, w, nharm=3):
    """y = a0 + sum_n [an cos(n w t) + bn sin(n w t)].

    Returns (mean, amp1, phase1) with the first harmonic written as
    amp1*cos(w t + phase1).
    """
    M = [np.ones_like(t)]
    for n in range(1, nharm + 1):
        M += [np.cos(n * w * t), np.sin(n * w * t)]
    c, *_ = np.linalg.lstsq(np.column_stack(M), y, rcond=None)
    return c[0], float(np.hypot(c[1], c[2])), float(np.arctan2(-c[2], c[1]))


def rule(title, width=78):
    print()
    print(f"  {title}")
    print("  " + "-" * (width - 2))


# ------------------------------------------------------------ sweep file ----
COLS = ["lam/L", "T", "Te", "F1", "F2", "F6",
        "eta1", "eta2", "eta3", "eta4", "eta5", "eta6"]


def append_row(path, row, case, head, U0, L, B, steep, depth):
    """Append one case to the sweep file, writing a header if it is new.

    A case already in the file is replaced rather than duplicated, so a rerun
    does not leave two rows for the same wavelength, and rows are kept sorted
    by lambda/L so the file plots straight out of the box.
    """
    header = [
        "# SOBC wave-load sweep",
        f"#   heading   {np.degrees(head) + 0.0:+.4g} deg"
        " (0 = head sea, waves along +x)",
        f"#   speed     U {U0:+.5g} m/s   Fn {U0/np.sqrt(G*L):.4g}",
        f"#   ship      L {L:.4g} m   B {B:.4g} m"
        "   (water depth is per case, 0.5 lambda)",
        f"#   wave      H/lambda {steep:.5g}",
        f"#   non-dim   F1,F2 / (rho g A^2 B^2 / L)   F6 / (rho g A^2 B^2)",
        "#             eta1..3 / A                   eta4..6 / (k A)",
        "#   axes      loads in ship axes (bow at -x, so F1 > 0 is added"
        " resistance)",
        "#",
        "#" + " ".join(f"{c:>11s}" for c in COLS)[1:],
    ]

    def key(s):
        try:
            return float(s.split()[0])
        except ValueError:
            return float("inf")

    # Rows are keyed on lambda/L -- one heading and speed per file, so it
    # identifies the case on its own now that the name is not carried.
    lam_L = row["lam/L"]
    old = []
    if os.path.exists(path):
        for line in open(path, errors="ignore"):
            f = line.split()
            if line.startswith("#") or len(f) != len(COLS):
                continue                          # blank, comment or old format
            if abs(key(line) - lam_L) > 1e-6*max(1.0, abs(lam_L)):
                old.append(line.rstrip("\n"))     # else: drop the previous run

    line = " ".join(f"{row.get(c, float('nan')):11.5g}" for c in COLS)

    with open(path, "w") as f:
        f.write("\n".join(header + sorted(old + [line], key=key)) + "\n")

    print(f"  wrote row for {case} (lambda/L {lam_L:.4g}) to {path}"
          f"  ({len(old)+1} case{'s' if old else ''} in the sweep)")


# ------------------------------------------------- collected csv output ----
CSV_COLS = (["lam/L", "U", "heading", "F1mean", "F2mean", "Mzmean"]
            + [f"eta{i}" for i in range(1, 7)])
CSV_KEY = 3                                # lam/L, U, heading identify a row


def write_csv(path, lam_L, U0, heading, row):
    """Append one case to a single csv collecting the whole sweep.

    Keyed on (lambda/L, U, heading) -- all three are needed.  A sweep over
    headings and speeds at one wavelength has lambda/L constant, so keying on
    it alone would make every case overwrite the last.

    heading is in degrees in the same convention as the case name and the
    runsim --heading argument (0 = head seas), i.e. -headingAngle.  Loads are
    non-dimensional and in ship axes, as printed above.  eta1..3 are divided
    by A and eta4..6 by kA; they are amplitudes, so the 180 deg between the
    solver body frame and the load axes does not affect them.  A restrained
    body writes zeros, and a case with no motion output writes nan.
    """
    if path.endswith(os.sep) or os.path.isdir(path):
        path = os.path.join(path, "meanLoads.csv")
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    me = (lam_L, U0, heading)

    def key(fields):
        return tuple(float(v) for v in fields[:CSV_KEY])

    def same(a, b):
        return all(abs(x - y) <= 1e-6*max(1.0, abs(x), abs(y))
                   for x, y in zip(a, b))

    rows = []
    if os.path.exists(path):
        for line in open(path, errors="ignore"):
            f = line.strip().split(",")
            if len(f) != len(CSV_COLS):
                continue
            try:
                k = key(f)
            except ValueError:
                continue                       # header row
            if not same(k, me):                # else: drop the previous run
                rows.append(line.rstrip("\n"))

    nan = float("nan")
    vals = ([lam_L, U0, heading, row.get("F1", nan), row.get("F2", nan),
             row.get("F6", nan)]
            + [row.get(f"eta{i}", nan) for i in range(1, 7)])
    rows.append(",".join(f"{v:.6g}" for v in vals))
    rows.sort(key=lambda line: key(line.split(",")))

    with open(path, "w") as f:
        f.write(",".join(CSV_COLS) + "\n")
        f.write("\n".join(rows) + "\n")

    print(f"  wrote {path}  ({len(rows)} cases collected)")


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nper", type=int, default=6,
                    help="encounter periods to average over")
    ap.add_argument("--nharm", type=int, default=3)
    ap.add_argument("--fo", default="meanLoads", help="drift-load function object")
    ap.add_argument("--append", default=None,
                    help="append one summary row to this file, for a sweep")
    ap.add_argument("--csv", default=None, metavar="PATH",
                    help="also collect lam/L,U,heading,F1,F2,Mz into one csv "
                         "(a directory gets meanLoads.csv inside it)")
    a = ap.parse_args()

    wc = read_dict("constant/waveConditions")
    bd = read_dict("constant/bodyMotionProperties")

    lam = wc["waveLength"]
    steep = wc["steepness"]
    U0 = wc["currentSpeed"]
    head = wc["headingAngle"]
    h = wc["waterDepth"]
    ramp_per = wc.get("rampPeriods", 0.0)

    L = bd.get("Lpp", 1.0)
    B = bd.get("beam", 1.0)

    k = 2.0 * np.pi / lam
    w0 = np.sqrt(G * k * np.tanh(k * h))
    we = w0 + k * U0 * np.cos(head)
    Te = 2.0 * np.pi / we
    A = 0.5 * steep * lam
    t_ramp = ramp_per * 2.0 * np.pi / w0

    # Seo et al. denominators: the moment is NOT divided by L
    den_F = RHO * G * A**2 * B**2 / L
    den_M = RHO * G * A**2 * B**2

    print()
    print("=" * 78)
    print(f"  {os.path.basename(os.getcwd())}")
    print("=" * 78)
    print(f"  wave      lambda {lam:.4g} m   lambda/L {lam/L:.4g}   A {A:.5g} m"
          f"   H/L {2*A/L:.5g}")
    print(f"            k {k:.4g} 1/m   omega {w0:.5g}   omega_e {we:.5g} rad/s"
          f"   Te {Te:.4g} s")
    print(f"  ship      L {L:.4g} m   B {B:.4g} m   U {U0:+.4g} m/s"
          f"   Fn {U0/np.sqrt(G*L):.4g}   tau {we*U0/G:+.4g}")
    print(f"  heading   {np.degrees(head)+0.0:+.4g} deg (0 = head sea)"
          f"   loads rotated into ship axes by {-np.degrees(head)+0.0:+.4g} deg")
    print(f"  denom     force  rho g A^2 B^2 / L = {den_F:.5g} N")
    print(f"            moment rho g A^2 B^2     = {den_M:.5g} N m")

    # --- window ----------------------------------------------------------
    F = load(f"postProcessing/{a.fo}/*/force.dat")
    M = load(f"postProcessing/{a.fo}/*/moment.dat")
    mot = load("postProcessing/bodyMotion/motion.dat")

    avail = [d for d in (F, M, mot) if d is not None]
    if not avail:
        sys.exit("\n  ERROR: no postProcessing output found -- has the case run?\n"
                 f"         looked for postProcessing/{a.fo}/*/force.dat and\n"
                 "         postProcessing/bodyMotion/motion.dat\n")
    tmax = max(d[-1, 0] for d in avail)
    n_avail = int(np.floor((tmax - t_ramp) / Te))
    if n_avail < 1:
        sys.exit(f"\n  ERROR: only {(tmax-t_ramp)/Te:.2f} encounter periods after "
                 f"the ramp (ends {t_ramp:.3g} s, data to {tmax:.3g} s).")
    nper = min(a.nper, n_avail)
    lo, hi = tmax - nper * Te, tmax
    print(f"  window    last {nper} of {n_avail} encounter periods:"
          f"  t = {lo:.4g} .. {hi:.4g} s")

    row = {"lam/L": lam / L, "T": 2.0 * np.pi / w0, "Te": Te}

    # --- drift loads ------------------------------------------------------
    if F is not None:
        rule("MEAN WAVE DRIFT LOAD        (Seo et al. Figs 11, 14, 16)")
        print(f"  {'':<12s} {'value':>13s} {'non-dim':>10s}     contributions"
              f" (surface / elevation / strip)")
        m = (F[:, 0] >= lo) & (F[:, 0] <= hi)

        # The function object integrates in mesh axes, where +x is the wave
        # propagation direction.  The loads are wanted in ship axes, so rotate
        # the horizontal pair through the heading.  In head seas this is a
        # 180 deg turn and only flips signs, which is why it went unnoticed;
        # at an oblique heading it mixes x and y.  Fz and Mz are about the
        # rotation axis and are unaffected.
        #
        # psi = -headingAngle puts the bow at -x', so F1 > 0 is added
        # resistance -- the convention meanVal.py and Seo et al. use.  Note it
        # is the solver's body frame turned 180 deg, so the surge/sway phases
        # in the RAO block below sit 180 deg from these axes.
        psi = -head
        cps, sps = np.cos(psi), np.sin(psi)

        def rot(cx, cy):
            """Mean of the (x, y) force pair in column cx, cy, in ship axes."""
            fx = fit(F[m, 0], F[m, cx], 2 * we, a.nharm)[0]
            fy = fit(F[m, 0], F[m, cy], 2 * we, a.nharm)[0]
            return fx * cps - fy * sps, fx * sps + fy * cps

        tots = {}
        parts = {}
        tots["F1"], tots["F2"] = rot(1, 2)
        parts["F1"], parts["F2"] = zip(*(rot(1 + o, 2 + o) for o in (3, 6, 9)))
        tots["F3"] = fit(F[m, 0], F[m, 3], 2 * we, a.nharm)[0]
        parts["F3"] = [fit(F[m, 0], F[m, 3 + o], 2 * we, a.nharm)[0]
                       for o in (3, 6, 9)]

        for name, key, unit in (("F1  (added R)", "F1", "N"),
                                ("F2  (sway)", "F2", "N"),
                                ("F3", "F3", "N")):
            tot = tots[key]
            print(f"  {name:<12s} {tot:+13.6g} {tot/den_F:+10.4f}     "
                  + " / ".join(f"{p/den_F:+7.4f}" for p in parts[key])
                  + f"   [{unit}]")
            row[key] = tot / den_F
        if M is not None:
            mm = (M[:, 0] >= lo) & (M[:, 0] <= hi)
            tot = fit(M[mm, 0], M[mm, 3], 2 * we, a.nharm)[0]
            parts = [fit(M[mm, 0], M[mm, 3 + off], 2 * we, a.nharm)[0]
                     for off in (3, 6, 9)]
            print(f"  {'F6  (yaw)':<12s} {tot:+13.6g} {tot/den_M:+10.4f}     "
                  + " / ".join(f"{p/den_M:+7.4f}" for p in parts) + "   [N m]")
            row["F6"] = tot / den_M
        print()
        print("  the total is a difference of larger, individually CV-dependent")
        print("  terms, so watch the contributions as well as the sum")

    # --- motion RAOs ------------------------------------------------------
    if mot is not None and mot.shape[1] >= 7:
        rule("MOTION RAO at omega_e       (Seo et al. Figs 9, 10)")
        print(f"  {'':<10s} {'amplitude':>12s} {'phase':>9s} {'RAO':>10s}"
              f"   {'normalised by':<14s}")
        m = (mot[:, 0] >= lo) & (mot[:, 0] <= hi)
        dofs = [("surge", 1, A, "A"), ("sway", 2, A, "A"), ("heave", 3, A, "A"),
                ("roll", 4, k * A, "kA"), ("pitch", 5, k * A, "kA"),
                ("yaw", 6, k * A, "kA")]
        if m.sum() < 20:
            print("  (not enough motion samples in the window)")
        else:
            for name, col, norm, lbl in dofs:
                mean, amp, ph = fit(mot[m, 0], mot[m, col], we, a.nharm)
                print(f"  {name:<10s} {amp:12.6g} {np.degrees(ph):+9.2f}"
                      f" {amp/norm:10.4f}   {lbl:<14s}"
                      + ("" if abs(mean) < 1e-3*max(amp, 1e-30)
                         else f"  mean {mean:+.4g}"))
                row[f"eta{col}"] = amp / norm
            print()
            print("  a restrained body (fixBody true) reports zeros here")

    # --- convergence ------------------------------------------------------
    if F is not None:
        rule("STABILITY  (F1 non-dim over successive 2-Te windows)")
        line = []
        for j in range(2, min(n_avail, 10) + 1):
            a1, b1 = hi - j * Te, hi - (j - 2) * Te
            mm = (F[:, 0] >= a1) & (F[:, 0] <= b1)
            if mm.sum() < 20:
                continue
            line.append(f"{fit(F[mm,0],F[mm,1],2*we,a.nharm)[0]/den_F:+.4f}")
        print("  " + "  ".join(reversed(line)) + "   (oldest -> newest)")

    print()

    if a.append:
        append_row(a.append, row, os.path.basename(os.getcwd()),
                   head, U0, L, B, 2 * A / lam, h)
        print()

    if a.csv:
        if all(kk in row for kk in ("F1", "F2", "F6")):
            write_csv(a.csv, lam / L, U0, -np.degrees(head) + 0.0, row)
            print()
        else:
            print("  no drift loads to write to the csv")
            print()


if __name__ == "__main__":
    main()
