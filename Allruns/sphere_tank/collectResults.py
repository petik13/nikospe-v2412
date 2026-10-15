#!/usr/bin/env python3
"""
Collect the mean wave loads of all SPHERE cases in a directory into one
table, results.csv (MMG axes, the format of the sphere's meanLoads.py --csv:
Zhao et al. 1988 non-dimensional form, keyed on omega_e^2 R/g, U/sqrt(gD), mu).

The rows are exactly what meanLoads.py writes: X, Y and N include the psi-bar
terms of runs whose psiBar.dat was written with the divergence-form hull
forcing (middleFieldForm since 2026-10-10), X and Y with the mass correction
rho Uinf (Phi_C + S_C) that takes the momentum of the residual of the
numerical psibar's mass balance back out.  Older psiBar runs (gradient-form
hull forcing, spurious psi-bar surge force) have the psi-bar yaw moment in N
only; runs without psiBar.dat have no psi-bar term at all.  Both groups are
listed in comment lines of the file, since at forward speed their X lacks the
return-flow term (see meanLoads.py).

A case is a subdirectory with constant/waveConditions and
postProcessing/<fo>/*/force.dat; baseCase and cases that have not run are
skipped.  Each case is processed once by meanLoads.py --fo <fo> --csv, run in
the case directory: the current baseCase/meanLoads.py if there is one next to
this script, in DIR or in its parent (so not a case's own, possibly older,
copy), else each case's own meanLoads.py.

The file is written fresh; an existing one is kept as results.csv.bak.  Rows
are keyed on (w2R/g, U/sqrt(gD), heading): if two cases share a key, the later one
(alphabetically) wins and a warning names both.

usage: run it in (or point it at) a directory of cases; the table is written
there:
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
    """(omega_e^2 R/g, U/sqrt(gD), mu) of the table row the sphere's meanLoads.py
    writes for this case (Zhao et al. 1988 form; D = beam, R = D/2), computed as
    it does."""
    g = 9.81
    h = -math.degrees(wc["headingAngle"]) + 0.0
    mu = (360.0 - h) % 360.0
    if mu > 180.0 + 1e-9:
        mu = 360.0 - mu
    k = 2.0 * math.pi / wc["waveLength"]
    w0 = math.sqrt(g * k * math.tanh(k * wc["waterDepth"]))
    head, U0 = wc["headingAngle"], wc["currentSpeed"]
    Ux = U0 * math.cos(head) - wc.get("swaySpeed", 0.0) * math.sin(head)
    we = w0 + k * Ux
    D = bd.get("beam", 1.0)
    return (round(float(f"{we**2 * 0.5 * D / g:.6g}"), 4),
            round(float(f"{U0 / math.sqrt(g * D):.6g}"), 4), round(mu, 4))


def read_rows(path):
    """{(w2R/g, U/sqrt(gD), heading): row fields} of the sphere meanLoads.py csv."""
    rows = {}
    if not os.path.isfile(path):
        return rows
    for line in open(path):
        f = line.strip().split(",")
        try:
            k = tuple(round(float(v), 4) for v in f[:3])   # as table_key; the csv has 6 digits
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
    out = os.path.join(root, a.out)
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
    print(f"  {'case':38s} {'w2R/g':>6s} {'U/vgD':>7s} {'mu':>5s} | {'X':>8s} {'Y':>8s}"
          f" {'N':>8s}  psibar / notes")
    seen, noPsi, oldPsi = {}, [], []
    for name, c in cases:
        wc, bd = case_dicts(c)
        key = table_key(wc, bd)
        if key in seen:
            print(f"  WARNING: {name} and {seen[key]} are the same table row {key};"
                  f" {name} wins")
        seen[key] = name
        mlc = ml or os.path.join(c, "meanLoads.py")
        if not os.path.isfile(mlc):
            print(f"  skip {name}: no {mlc}; give --meanLoads")
            continue
        r = subprocess.run([sys.executable, mlc, "--fo", a.fo, "--csv", tmp],
                           cwd=c, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  {name:38s} meanLoads.py failed: "
                  + (r.stderr or r.stdout).strip().splitlines()[-1])
            continue
        if "psiBar.dat at t" not in r.stdout:
            psi = "none"
            noPsi.append(name)
        elif "psibar forces NOT added" in r.stdout:
            psi = "N only (old psiBar.dat)"
            oldPsi.append(name)
        else:
            psi = "X, Y, N"
        notes = [psi] + [l.strip() for l in r.stdout.splitlines() if "WARNING" in l]
        row = read_rows(tmp).get(key)
        if row:
            X, Y, N = float(row[3]), float(row[4]), float(row[5])
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |"
                  f" {X:+8.3f} {Y:+8.3f} {N:+8.3f}  " + "; ".join(notes))
        else:
            print(f"  {name:38s} {key[0]:6.3f} {key[1]:7.4f} {key[2]:5.0f} |  no row  "
                  + "; ".join(notes))

    if not os.path.isfile(tmp):
        sys.exit("no rows written")
    lines = open(tmp).read().splitlines()
    os.remove(tmp)

    note = ["# psibar: F1mean (X), F2mean (Y), Mzmean (N) include the psi-bar terms"
            " (divergence-form psiBar.dat; collectResults.py), X and Y with the mass"
            " correction rho Uinf (Phi_C + S_C)"]
    if oldPsi:
        note.append("# old psiBar.dat (gradient-form hull forcing): N includes psi-bar,"
                    f" X, Y do not: {', '.join(oldPsi)}")
    if noPsi:
        note.append(f"# no psiBar.dat (no psi-bar term at all): {', '.join(noPsi)}")
    write_table(lines, out, note)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
