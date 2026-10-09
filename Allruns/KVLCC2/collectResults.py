#!/usr/bin/env python3
"""
Collect the mean wave loads of all cases in a directory into two manModel
tables (MMG axes, the format of meanLoads.py --csv):

    results.csv      X, Y without the psi-bar terms; N with the psi-bar yaw moment
    results_psi.csv  the same, except that X includes the psi-bar surge force

(Why X is kept apart: see the docstring of meanLoads.py -- the psi-bar surge
force comes out as -rho U Q, and the net psibar source Q of the runs so far is
too large to be physical.)

A case is a subdirectory with constant/waveConditions and
postProcessing/<fo>/*/force.dat; baseCase and cases that have not run are
skipped.  Each case is processed once by meanLoads.py --fo <fo> --csv, run in
the case directory: the current baseCase/meanLoads.py if there is one next to
this script, in DIR or in its parent (so not a case's own, possibly older,
copy), else each case's own meanLoads.py.  That gives results.csv.

results_psi.csv is results.csv with the psi-bar surge force added to F1mean.
It is read here from the last row of postProcessing/<fo>/*/psiBar.dat and
converted exactly as meanLoads.py does (mesh axes -> MMG X, divided by
rho g A^2 B^2 / L), so any meanLoads.py version that handles psiBar.dat works.

Both files are written fresh; existing ones are kept as <name>.bak.  A case
without psiBar.dat (psiBar off) gives the same row in both files and no
psi-bar term at all; they are listed in a comment line of the files.  Rows are
keyed on (lam/L, U, heading): if two cases share a key, the later one
(alphabetically) wins and a warning names both.

usage: run it in (or point it at) a directory of cases; the two tables are
written there:
    cd <directory of cases> && python3 collectResults.py      # a copy of it, or ../collectResults.py
    python3 collectResults.py <directory of cases>
    python3 collectResults.py --fo meanLoads2                  # another load function object
"""
import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RHO, G = 1000.0, 9.81


def read_dict(path):
    vals = {}
    for line in open(path, errors="ignore"):
        line = line.split("//")[0].strip()
        m = re.match(r"^([A-Za-z_]\w*)\s+(-?[\d.eE+-]+)\s*;", line)
        if m:
            try:
                vals[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
    return vals


def case_dicts(case):
    return (read_dict(os.path.join(case, "constant", "waveConditions")),
            read_dict(os.path.join(case, "constant", "bodyMotionProperties")))


def table_key(wc, bd):
    """(lam/L, U, mu) of the table row meanLoads.py writes for this case."""
    h = -math.degrees(wc["headingAngle"]) + 0.0
    mu = (360.0 - h) % 360.0
    if mu > 180.0 + 1e-9:
        mu = 360.0 - mu
    return (round(wc["waveLength"] / bd.get("Lpp", 1.0), 6),
            round(wc["currentSpeed"], 6), round(mu, 6))


def psi_surge(case, fo, wc, bd):
    """psi-bar surge force X_psi / (rho g A^2 B^2 / L) from the last row of
    psiBar.dat, as meanLoads.py: X = -(Fx cos psi - Fy sin psi), psi = -heading.
    None without psiBar.dat."""
    rows = []
    for p in glob.glob(os.path.join(case, "postProcessing", fo, "*", "psiBar.dat")):
        for line in open(p, errors="ignore"):
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(v) for v in line.split()])
            except ValueError:
                pass
    rows = [r for r in rows if len(r) >= 8]
    if not rows:
        return None
    r = max(rows, key=lambda r: r[0])
    psi = -wc["headingAngle"]
    X = -(r[2] * math.cos(psi) - r[3] * math.sin(psi))
    A = 0.5 * wc["steepness"] * wc["waveLength"]
    L, B = bd.get("Lpp", 1.0), bd.get("beam", 1.0)
    return X / (RHO * G * A**2 * B**2 / L)


def read_rows(path):
    """{(lam/L, U, heading): row fields} of a meanLoads.py csv."""
    rows = {}
    if not os.path.isfile(path):
        return rows
    for line in open(path):
        f = line.strip().split(",")
        try:
            k = tuple(round(float(v), 6) for v in f[:3])
        except ValueError:
            continue
        rows[k] = f
    return rows


def write_table(lines, out, note):
    """Write lines to out (backup of an existing out), with the note after the
    '# convention: MMG' marker line."""
    lines = list(lines)
    i = next((j for j, l in enumerate(lines) if l.startswith("# convention: MMG")), 0)
    lines[i + 1:i + 1] = note
    if os.path.isfile(out):
        shutil.copy2(out, out + ".bak")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", nargs="?", default=".", help="directory with the cases")
    ap.add_argument("--fo", default="meanLoads", help="drift-load function object")
    ap.add_argument("--meanLoads", default=None,
                    help="meanLoads.py to use (default: baseCase/meanLoads.py next to this"
                         " script or in DIR or its parent, else each case's own copy)")
    ap.add_argument("--out", default="results.csv")
    ap.add_argument("--outPsi", default="results_psi.csv")
    a = ap.parse_args()

    root = os.path.abspath(a.dir)
    if a.meanLoads:
        ml = os.path.abspath(a.meanLoads)
        if not os.path.isfile(ml):
            sys.exit(f"no {ml}")
    else:
        cand = [os.path.join(d, "baseCase", "meanLoads.py")
                for d in (HERE, root, os.path.dirname(root))]
        ml = next((p for p in cand if os.path.isfile(p)), None)   # None: per case
    out, outPsi = os.path.join(root, a.out), os.path.join(root, a.outPsi)
    tmp = out + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)

    cases = []
    for name in sorted(os.listdir(root)):
        c = os.path.join(root, name)
        if name == "baseCase" or not os.path.isdir(c):
            continue
        if not os.path.isfile(os.path.join(c, "constant", "waveConditions")):
            continue
        if not glob.glob(os.path.join(c, "postProcessing", a.fo, "*", "force.dat")):
            print(f"  skip {name}: no postProcessing/{a.fo}/*/force.dat")
            continue
        cases.append((name, c))
    if not cases:
        sys.exit(f"no cases with postProcessing/{a.fo} in {root}")

    print(f"collecting {len(cases)} cases from {root} with "
          + (ml if ml else "each case's own meanLoads.py"))
    print(f"  {'case':38s} {'lam/L':>6s} {'U':>7s} {'mu':>5s} | {'X':>8s} {'X_psi':>8s}"
          f" {'Y':>8s} {'N':>8s}  notes")
    seen, noPsi, xpsi = {}, [], {}
    for name, c in cases:
        wc, bd = case_dicts(c)
        key = table_key(wc, bd)
        if key in seen:
            print(f"  WARNING: {name} and {seen[key]} are the same table row {key};"
                  f" {name} wins")
        seen[key] = name
        notes = []
        mlc = ml or os.path.join(c, "meanLoads.py")
        if not os.path.isfile(mlc):
            print(f"  skip {name}: no {mlc}; give --meanLoads")
            continue
        if "psiBar.dat" not in open(mlc, errors="ignore").read():
            notes.append("meanLoads.py too old: N without psibar")
        r = subprocess.run([sys.executable, mlc, "--fo", a.fo, "--csv", tmp],
                           cwd=c, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  {name:38s} meanLoads.py failed: "
                  + (r.stderr or r.stdout).strip().splitlines()[-1])
            continue
        notes += [l.strip() for l in r.stdout.splitlines() if "WARNING" in l]
        xp = psi_surge(c, a.fo, wc, bd)
        if xp is None:
            noPsi.append(name)
            notes.append("no psiBar.dat")
            xpsi.pop(key, None)
        else:
            xpsi[key] = xp
        row = read_rows(tmp).get(key)
        if row:
            X, Y, N = float(row[3]), float(row[4]), float(row[5])
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |"
                  f" {X:+8.3f} {X + (xp or 0.0):+8.3f} {Y:+8.3f} {N:+8.3f}  " + "; ".join(notes))
        else:
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |  no row  "
                  + "; ".join(notes))

    if not os.path.isfile(tmp):
        sys.exit("no rows written")
    lines = open(tmp).read().splitlines()
    os.remove(tmp)

    # results_psi.csv: F1mean (X) + X_psi on the rows of cases with psiBar.dat
    linesPsi = []
    for line in lines:
        f = line.split(",")
        try:
            k = tuple(round(float(v), 6) for v in f[:3])
        except ValueError:
            linesPsi.append(line)
            continue
        if k in xpsi:
            f[3] = f"{float(f[3]) + xpsi[k]:.6g}"
        linesPsi.append(",".join(f))

    common = ([f"# cases without psiBar.dat (no psi-bar term at all): {', '.join(noPsi)}"]
              if noPsi else [])
    write_table(lines, out, [
        "# psibar: F1mean (X), F2mean (Y) without the psi-bar terms; Mzmean (N) includes"
        " the psi-bar yaw moment (collectResults.py)"] + common)
    write_table(linesPsi, outPsi, [
        "# psibar: F1mean (X) includes the psi-bar surge force; F2mean (Y) does not;"
        " Mzmean (N) includes the psi-bar yaw moment (collectResults.py)"] + common)
    print(f"wrote {out}")
    print(f"wrote {outPsi}")


if __name__ == "__main__":
    main()
