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
skipped.  Each case is processed by meanLoads.py -- by default the current
baseCase/meanLoads.py next to this script, not the case's own (possibly older)
copy -- run in the case directory:

    meanLoads.py --fo <fo> --csv <tmp>                        -> results.csv
    meanLoads.py --fo <fo> --psiBarForces X --csv <tmp_psi>   -> results_psi.csv

Both files are written fresh; existing ones are kept as <name>.bak.  A case
without psiBar.dat (psiBar off) gives the same row in both files and no
psi-bar term at all; they are listed in a comment line of the files.  Rows are
keyed on (lam/L, U, heading): if two cases share a key, the later one
(alphabetically) wins and a warning names both.

usage (from Allruns/KVLCC2):
    python3 collectResults.py                  # cases in the current directory
    python3 collectResults.py runs/sweep90     # cases in another directory
    python3 collectResults.py --fo meanLoads2  # another load function object
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


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


def table_key(case):
    """(lam/L, U, mu) of the table row meanLoads.py writes for this case."""
    import math
    wc = read_dict(os.path.join(case, "constant", "waveConditions"))
    bd = read_dict(os.path.join(case, "constant", "bodyMotionProperties"))
    h = -math.degrees(wc["headingAngle"]) + 0.0
    mu = (360.0 - h) % 360.0
    if mu > 180.0 + 1e-9:
        mu = 360.0 - mu
    return (round(wc["waveLength"] / bd.get("Lpp", 1.0), 6),
            round(wc["currentSpeed"], 6), round(mu, 6))


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


def finish(tmp, out, note):
    """Move tmp to out (backup of an existing out), with a note after the
    '# convention: MMG' marker line."""
    if not os.path.isfile(tmp):
        return False
    lines = open(tmp).read().splitlines()
    i = next((j for j, l in enumerate(lines) if l.startswith("# convention: MMG")), 0)
    lines[i + 1:i + 1] = note
    if os.path.isfile(out):
        shutil.copy2(out, out + ".bak")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.remove(tmp)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", nargs="?", default=".", help="directory with the cases")
    ap.add_argument("--fo", default="meanLoads", help="drift-load function object")
    ap.add_argument("--meanLoads", default=os.path.join(HERE, "baseCase", "meanLoads.py"),
                    help="meanLoads.py to use (default: baseCase/meanLoads.py)")
    ap.add_argument("--out", default="results.csv")
    ap.add_argument("--outPsi", default="results_psi.csv")
    a = ap.parse_args()

    root = os.path.abspath(a.dir)
    ml = os.path.abspath(a.meanLoads)
    if not os.path.isfile(ml):
        sys.exit(f"no {ml}")
    out, outPsi = os.path.join(root, a.out), os.path.join(root, a.outPsi)
    tmp, tmpPsi = out + ".tmp", outPsi + ".tmp"
    for p in (tmp, tmpPsi):
        if os.path.exists(p):
            os.remove(p)

    cases = []
    for name in sorted(os.listdir(root)):
        c = os.path.join(root, name)
        if name == "baseCase" or not os.path.isdir(c):
            continue
        if not os.path.isfile(os.path.join(c, "constant", "waveConditions")):
            continue
        fo_dir = os.path.join(c, "postProcessing", a.fo)
        if not (os.path.isdir(fo_dir) and any(
                os.path.isfile(os.path.join(fo_dir, t, "force.dat")) for t in os.listdir(fo_dir))):
            print(f"  skip {name}: no postProcessing/{a.fo}/*/force.dat")
            continue
        cases.append((name, c))
    if not cases:
        sys.exit(f"no cases with postProcessing/{a.fo} in {root}")

    print(f"collecting {len(cases)} cases from {root} with {ml}")
    print(f"  {'case':38s} {'lam/L':>6s} {'U':>7s} {'mu':>5s} | {'X':>8s} {'X_psi':>8s} {'Y':>8s} {'N':>8s}  notes")
    seen, noPsi = {}, []
    for name, c in cases:
        key = table_key(c)
        if key in seen:
            print(f"  WARNING: {name} and {seen[key]} are the same table row {key};"
                  f" {name} wins")
        seen[key] = name
        notes = []
        for extra, dst in (([], tmp), (["--psiBarForces", "X"], tmpPsi)):
            r = subprocess.run([sys.executable, ml, "--fo", a.fo, "--csv", dst] + extra,
                               cwd=c, capture_output=True, text=True)
            if r.returncode != 0:
                notes.append("meanLoads.py failed: " + (r.stderr or r.stdout).strip().splitlines()[-1])
                break
            if not extra:
                if "psiBar.dat at t" not in r.stdout:
                    noPsi.append(name)
                    notes.append("no psiBar.dat")
                notes += [l.strip() for l in r.stdout.splitlines() if "WARNING" in l]
        row, rowP = read_rows(tmp).get(key), read_rows(tmpPsi).get(key)
        if row and rowP:
            X, Y, N, Xp = float(row[3]), float(row[4]), float(row[5]), float(rowP[3])
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |"
                  f" {X:+8.3f} {Xp:+8.3f} {Y:+8.3f} {N:+8.3f}  " + "; ".join(notes))
        else:
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |  no row  "
                  + "; ".join(notes))

    common = ([f"# cases without psiBar.dat (no psi-bar term at all): {', '.join(noPsi)}"]
              if noPsi else [])
    ok1 = finish(tmp, out, [
        "# psibar: F1mean (X), F2mean (Y) without the psi-bar terms; Mzmean (N) includes"
        " the psi-bar yaw moment (collectResults.py)"] + common)
    ok2 = finish(tmpPsi, outPsi, [
        "# psibar: F1mean (X) includes the psi-bar surge force; F2mean (Y) does not;"
        " Mzmean (N) includes the psi-bar yaw moment (collectResults.py)"] + common)
    print(f"wrote {out}" if ok1 else "no rows for results.csv")
    print(f"wrote {outPsi}" if ok2 else "no rows for results_psi.csv")


if __name__ == "__main__":
    main()
