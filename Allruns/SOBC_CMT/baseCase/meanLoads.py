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

Axes: EVERYTHING REPORTED IS IN MMG AXES (x forward, y starboard, z down,
origin at midship), the convention of manModel:

    X > 0 forward (added resistance is X < 0), Y > 0 starboard, Z > 0 down,
    N about z down, > 0 turns the bow to starboard;
    encounter angle mu = (360 - h) mod 360, 0 = head sea, 90 = waves from
    starboard, 270 = waves from port, with h the runsim --heading (0 = head
    sea, waves along +x; h = 90 has the waves from port).

The function object integrates in mesh axes, where +x is the wave
propagation direction and the hull's bow points along (-cos h, sin h).  With
F1 = Fx cos h - Fy sin h (aft) and F2 = Fx sin h + Fy cos h (starboard):

    X = -F1,   Y = F2,   Z = -Fz,   N = -Mz.

The printed summary is for the run itself (encounter angle mu).  The sweep
file (--append) and the collected csv (--csv) are tables for manModel and are
written for waves from starboard, mu in [0, 180]: a run with the waves from
port (mu > 180) is written as its mirror image, at 360 - mu with Y and N
sign-flipped (exact for the port/starboard-symmetric hull).  So a sweep over
h = 0..180 gives the table mu = 0..180.  Both files carry the marker line
"# convention: MMG ..."; an existing file without it (written by an older
meanLoads.py, F1 > 0 added resistance, heading h) is converted once, after a
backup to <file>.oldConvention, before the new row is added.

Seo et al. plot the added resistance, which is -X.

psi-bar.  When the drift-load function object ran with psiBar on
(middleFieldForm), postProcessing/<fo>/*/psiBar.dat holds the contribution of
the mean second-order potential (Grue & Palm 1993), averaged over the same
last 6 encounter periods.  Its last row is added to X, Y and N: the printed
totals, the sweep file and the csv include it (--noPsiBar leaves it out).  Z
is not corrected -- the vertical psi-bar force needs psibar on the hull itself.
The contribution is listed as a fourth part next to surface / elevation /
strip, with its free-surface and hull parts.  Derivation:
researchQuestions/psibarMidfield.ipynb.

The motion RAOs are amplitudes and do not depend on the axes.

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

MMG_MARKER = ("# convention: MMG body axes (x forward, y starboard, z down, midship);"
              " heading = encounter angle mu = (waveDirection + 180 - psi) mod 360 [deg],"
              " 0 head sea, 90 waves from starboard;"
              " F1mean = X > 0 forward, F2mean = Y > 0 starboard,"
              " Mzmean = N > 0 bow to starboard; U, V surge and sway at midship (V > 0 starboard)")


def to_table(X, Y, N, h_deg):
    """MMG loads of a run at runsim heading h -> (mu, X, Y, N) for a table
    written for waves from starboard: a run with the waves from port
    (mu = 360 - h > 180) is written as its mirror image (360 - mu, X, -Y, -N)."""
    mu = (360.0 - h_deg) % 360.0
    if mu > 180.0 + 1e-9:
        return 360.0 - mu, X, -Y, -N
    return mu, X, Y, N


def old_to_table(F1, F2, F6, h_deg):
    """Row of an older meanLoads.py file (F1 > 0 added resistance, F2 > 0
    starboard, F6 about z up, heading h) -> (mu, X, Y, N) in table form."""
    return to_table(-F1, F2, -F6, h_deg)


def has_marker(path):
    for line in open(path, errors="ignore"):
        if line.startswith("#") and "convention: MMG" in line:
            return True
    return False


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
    h_deg = -np.degrees(head) + 0.0
    mu_tab = to_table(0.0, 0.0, 0.0, h_deg)[0]
    header = [
        "# SOBC wave-load sweep",
        f"#   heading   h {h_deg:+.4g} deg (runsim; 0 = head sea, waves along +x)"
        f"  ->  table encounter angle mu {mu_tab:.4g} deg",
        f"#   speed     U {U0:+.5g} m/s   Fn {U0/np.sqrt(G*L):.4g}",
        f"#   ship      L {L:.4g} m   B {B:.4g} m"
        "   (water depth is per case, 0.5 lambda)",
        f"#   wave      H/lambda {steep:.5g}",
        "#   psibar    F1, F2, F6 include the psi-bar terms of runs that computed them"
        " (psiBar.dat; see summary.txt of each case)",
        f"#   non-dim   F1,F2 / (rho g A^2 B^2 / L)   F6 / (rho g A^2 B^2)",
        "#             eta1..3 / A                   eta4..6 / (k A)",
        "#   axes      MMG: F1 = X > 0 forward (added resistance is -X), F2 = Y > 0"
        " starboard, F6 = N > 0 bow to starboard; table form, waves from starboard",
        MMG_MARKER,
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
    if os.path.exists(path) and not has_marker(path):
        # written by an older meanLoads.py: convert its rows to MMG table form
        import shutil
        shutil.copy2(path, path + ".oldConvention")
        h_file = h_deg
        for line in open(path, errors="ignore"):
            m_ = re.match(r"^#\s+heading\s+([-+\d.eE]+)\s+deg", line)
            if m_:
                h_file = -float(m_.group(1))        # old header printed -h
        conv = []
        for line in open(path, errors="ignore"):
            f = line.split()
            if line.startswith("#") or len(f) != len(COLS):
                continue
            v_ = [float(x) for x in f]
            _, v_[3], v_[4], v_[5] = old_to_table(v_[3], v_[4], v_[5], h_file)
            conv.append(" ".join(f"{x:11.5g}" for x in v_))
        with open(path, "w") as fo:
            fo.write("\n".join(conv) + ("\n" if conv else ""))
        print(f"  {path}: old convention, {len(conv)} rows converted to MMG"
              f" (backup {path}.oldConvention)")
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

    heading is the MMG encounter angle mu [deg] of the table (waves from
    starboard, see the module docstring); loads are non-dimensional and in MMG
    axes, F1mean = X, F2mean = Y, Mzmean = N.  eta1..3 are divided
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
    if os.path.exists(path) and not has_marker(path):
        # written by an older meanLoads.py: convert its rows to MMG table form
        import shutil
        shutil.copy2(path, path + ".oldConvention")
        conv = []
        for line in open(path, errors="ignore"):
            f = line.strip().split(",")
            if len(f) != len(CSV_COLS):
                continue
            try:
                v_ = [float(x) for x in f]
            except ValueError:
                continue                       # header row
            v_[2], v_[3], v_[4], v_[5] = old_to_table(v_[3], v_[4], v_[5], v_[2])
            conv.append(",".join(f"{x:.6g}" for x in v_))
        with open(path, "w") as fo:
            fo.write(",".join(CSV_COLS) + "\n" + MMG_MARKER + "\n"
                     + "\n".join(conv) + ("\n" if conv else ""))
        print(f"  {path}: old convention, {len(conv)} rows converted to MMG"
              f" (backup {path}.oldConvention)")
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
        f.write(MMG_MARKER + "\n")
        f.write("\n".join(rows) + "\n")

    print(f"  wrote {path}  ({len(rows)} cases collected)")


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nper", type=int, default=6,
                    help="encounter periods to average over")
    ap.add_argument("--nharm", type=int, default=3)
    ap.add_argument("--fo", default="meanLoads", help="drift-load function object")
    ap.add_argument("--noPsiBar", action="store_true",
                    help="leave the psi-bar terms (psiBar.dat) out of the totals")
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
    Usway = wc.get("swaySpeed", 0)
    h = wc["waterDepth"]
    ramp_per = wc.get("rampPeriods", 0.0)

    L = bd.get("Lpp", 1.0)
    B = bd.get("beam", 1.0)

    k = 2.0 * np.pi / lam
    w0 = np.sqrt(G * k * np.tanh(k * h))
    # Uinf as the solver builds it: surge along the heading, sway 90 deg to it
    Ux = U0*np.cos(head) - Usway*np.sin(head)
    Uy = U0*np.sin(head) + Usway*np.cos(head)
    Umag = np.hypot(Ux, Uy)

    we = w0 + k * Ux
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
    # psi-bar terms: one row per field write and one at the end of the run
    Pb = None if a.noPsiBar else load(f"postProcessing/{a.fo}/*/psiBar.dat")

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

    # --- drift loads (MMG axes) -------------------------------------------
    h_deg = -np.degrees(head) + 0.0           # runsim heading, 0 = head sea
    mu_run = (360.0 - h_deg) % 360.0          # MMG encounter angle of this run
    if F is not None:
        rule(f"MEAN WAVE DRIFT LOAD, MMG axes, encounter angle mu {mu_run:.4g} deg")
        print(f"  {'':<14s} {'value':>13s} {'non-dim':>10s}     contributions"
              f" (surface / elevation / strip"
              + (" / psibar)" if Pb is not None and Pb.shape[1] >= 8 else ")"))
        m = (F[:, 0] >= lo) & (F[:, 0] <= hi)

        # The function object integrates in mesh axes, where +x is the wave
        # propagation direction and the bow points along (-cos h, sin h).
        # F1 = Fx cos h - Fy sin h is aft, F2 = Fx sin h + Fy cos h starboard;
        # MMG: X = -F1, Y = F2, Z = -Fz, N = -Mz.
        psi = -head
        cps, sps = np.cos(psi), np.sin(psi)

        def rot(cx, cy):
            """Mean (X, Y) in MMG axes of the force pair in columns cx, cy."""
            fx = fit(F[m, 0], F[m, cx], 2 * we, a.nharm)[0]
            fy = fit(F[m, 0], F[m, cy], 2 * we, a.nharm)[0]
            return -(fx * cps - fy * sps), fx * sps + fy * cps

        tots = {}
        parts = {}
        tots["X"], tots["Y"] = rot(1, 2)
        parts["X"], parts["Y"] = zip(*(rot(1 + o, 2 + o) for o in (3, 6, 9)))
        parts["X"], parts["Y"] = list(parts["X"]), list(parts["Y"])
        noPsi = {"X": tots["X"], "Y": tots["Y"]}

        # psi-bar terms (last row of psiBar.dat), mesh axes -> MMG as above
        psiB = None
        if Pb is not None and Pb.shape[1] >= 8:
            r_ = Pb[-1]
            psiB = {"t": r_[0], "per": r_[1],
                    "X": -(r_[2] * cps - r_[3] * sps), "Y": r_[2] * sps + r_[3] * cps,
                    "N": -r_[7],
                    "Nfs": -r_[10] if Pb.shape[1] > 13 else np.nan,
                    "Nhull": -r_[13] if Pb.shape[1] > 13 else np.nan,
                    "Q": (r_[14], r_[15]) if Pb.shape[1] > 15 else (np.nan, np.nan)}
            for key in ("X", "Y"):
                tots[key] += psiB[key]
                parts[key].append(psiB[key])
        tots["Z"] = -fit(F[m, 0], F[m, 3], 2 * we, a.nharm)[0]
        parts["Z"] = [-fit(F[m, 0], F[m, 3 + o], 2 * we, a.nharm)[0]
                      for o in (3, 6, 9)]

        for name, key in (("X  (surge)", "X"), ("Y  (sway)", "Y"), ("Z  (heave)", "Z")):
            tot = tots[key]
            print(f"  {name:<14s} {tot:+13.6g} {tot/den_F:+10.4f}     "
                  + " / ".join(f"{p/den_F:+7.4f}" for p in parts[key])
                  + (" /       -" if psiB is not None and key == "Z" else "")
                  + "   [N]")
        print(f"  {'':<14s} added resistance -X/den = {-tots['X']/den_F:+.4f}")
        Nz = None
        if M is not None:
            mm = (M[:, 0] >= lo) & (M[:, 0] <= hi)
            Nz = -fit(M[mm, 0], M[mm, 3], 2 * we, a.nharm)[0]
            parts = [-fit(M[mm, 0], M[mm, 3 + off], 2 * we, a.nharm)[0]
                     for off in (3, 6, 9)]
            noPsi["N"] = Nz
            if psiB is not None:
                Nz += psiB["N"]
                parts.append(psiB["N"])
            print(f"  {'N  (yaw)':<14s} {Nz:+13.6g} {Nz/den_M:+10.4f}     "
                  + " / ".join(f"{p/den_M:+7.4f}" for p in parts) + "   [N m]")
        if psiB is not None:
            print(f"  psibar         psiBar.dat at t = {psiB['t']:.4g} s, averaged over"
                  f" {psiB['per']:.2f} encounter periods; N = free surface"
                  f" {psiB['Nfs']/den_M:+.4f} + hull {psiB['Nhull']/den_M:+.4f}")
            print(f"                 net flux Q: hull {psiB['Q'][0]:+.3g}, free surface"
                  f" {psiB['Q'][1]:+.3g} m^3/s;  without psibar: X {noPsi['X']/den_F:+.4f}"
                  f"  Y {noPsi['Y']/den_F:+.4f}"
                  + (f"  N {noPsi['N']/den_M:+.4f}" if 'N' in noPsi else ""))
            if abs(psiB["t"] - hi) > 0.5 * Te:
                print(f"  WARNING: psiBar.dat ends at t = {psiB['t']:.4g} s, the fit window"
                      f" at {hi:.4g} s")
            if psiB["per"] < nper - 0.5:
                print(f"  WARNING: psibar averaged over {psiB['per']:.2f} periods only"
                      f" (fit window {nper})")
        elif not a.noPsiBar and abs(U0) > 0:
            print("  psibar         not computed (psiBar off in the function object):"
                  " N lacks the O(U) psi-bar term")

        # table form (waves from starboard) for the sweep file and the csv
        mu_tab, Xt, Yt, Nt = to_table(tots["X"] / den_F, tots["Y"] / den_F,
                                      (Nz / den_M) if Nz is not None else np.nan, h_deg)
        row["mu"] = mu_tab
        row["F1"], row["F2"] = Xt, Yt
        if Nz is not None:
            row["F6"] = Nt
        if abs(mu_tab - mu_run) > 1e-6:
            print(f"  table form (waves from starboard): mu {mu_tab:.4g} deg,"
                  f"  X {Xt:+.4f}  Y {Yt:+.4f}  N {Nt:+.4f}  (mirror image of this run)")
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
        rule("STABILITY  (X non-dim, MMG, over successive 2-Te windows)")
        line = []
        for j in range(2, min(n_avail, 10) + 1):
            a1, b1 = hi - j * Te, hi - (j - 2) * Te
            mm = (F[:, 0] >= a1) & (F[:, 0] <= b1)
            if mm.sum() < 20:
                continue
            fx = fit(F[mm, 0], F[mm, 1], 2 * we, a.nharm)[0]
            fy = fit(F[mm, 0], F[mm, 2], 2 * we, a.nharm)[0]
            line.append(f"{-(fx*np.cos(-head) - fy*np.sin(-head))/den_F:+.4f}")
        print("  " + "  ".join(reversed(line)) + "   (oldest -> newest)")

    print()

    if a.append:
        append_row(a.append, row, os.path.basename(os.getcwd()),
                   head, U0, L, B, 2 * A / lam, h)
        print()

    if a.csv:
        if all(kk in row for kk in ("F1", "F2", "F6")):
            write_csv(a.csv, lam / L, U0, row["mu"], row)
            print()
        else:
            print("  no drift loads to write to the csv")
            print()


if __name__ == "__main__":
    main()
