#!/usr/bin/env python3
"""
Sweep of circular motion tests (runCMT.py) over a list of (u, v, r), and the
collection of all their results.csv into one manModel wave table.

Sweep file (e.g. sweep_lam07.dat): 'key value' lines and one 'run' line per
test; '#' starts a comment.

    lam       2.24                      wavelength [m], all runs
    straight  waveData07.dat            straight-course table (r = 0, MMG) added
                                        to the collected table; 'none' to leave out
    out       sweep_lam07/waveData.dat  collected table
    run  u v r  [runCMT.py options]     MMG at midship, model scale: u, v [m/s]
                                        (v > 0 starboard), r [rad/s] (> 0 starboard)

Paths are relative to the sweep file.  Each run is a copy of baseCase in
<sweep file name>/u<u>_v<v>_r<r>/ (e.g. sweep_lam07/u0.25_v-0.15_r0.05), run as

    python3 runCMT.py --lam <lam> --u <u> --v <v> --r <r> [options]

which meshes, runs manFlowPrescribed, removes the processor directories and
writes results.csv (cmtPost.py) and timeSeries.png (plotSeries.py).  Its
output goes to runSweep.log in the run directory.

A run is complete when its results.csv covers a full turn: head seas at the
yaw onset and again at the end (360 and 0 for a starboard turn), or every
10 deg row but the end one (runs whose endTime had no margin stop a hair
short of it).  Complete runs are skipped, so the script resumes after an
interruption; incomplete ones are run again in place.  --force reruns the
complete ones too.

Collection (after the runs, or alone with --collect):
  * every complete run's results.csv (MMG).  Its head-sea row at the yaw onset
    is transitional (the filtered loads there mix the straight run and the
    yaw-rate ramp, but carry the r of the steady turn), so it is replaced by
    the head-sea row after the full turn (--keepOnset keeps it).  Without an
    end-of-turn row the onset row is used for both (same heading).
  * the straight-course table, if given: its rows at the sweep's lam/L, with
    V = 0 and r = 0 where it has no such columns.  A half table (0..180,
    waves from starboard) is completed to 0..360 by the mirror image
    (360 - mu, with V, r, Y, N sign-flipped).  Headings missing at one speed
    are taken from the nearest speed that has them.
  * mu = 360 as a copy of mu = 0 (or the other way round) for every (U, V, r),
    so that each heading row closes periodically in manModel's grid.
  * written in the manModel waveData format with the MMG marker line, sorted
    by U, V, r, heading.  In the manModel case: modulesDict 'waves:
    waveLoads360', waveDict 'swaySpeedInterp: true' (yawRateInterp is on by
    default).  manModel interpolates linearly in U, V, r and the encounter
    angle; (U, V, r) states without data are filled along U, then V, then r.

usage (from Allruns/KVLCC2_CMT):
    python3 runSweep.py sweep_lam07.dat --dryRun     # the runs, their length, the commands
    python3 runSweep.py sweep_lam07.dat              # run what is missing, then collect
    python3 runSweep.py sweep_lam07.dat --only 1 3   # only these runs (numbers or names)
    python3 runSweep.py sweep_lam07.dat --collect    # collect only
Options runSweep.py does not know are passed on to every runCMT.py call, e.g.
    python3 runSweep.py sweep_lam07.dat --nproc 56 --procD 8 6 1
"""
import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE / "baseCase"

# Hull at the scale of runCMT.py (KVLCC2, 1/100); only for the run-length estimate
L_PP, DRAFT, G = 3.2, 0.208, 9.81

COLS = (["lam/L", "U", "V", "r", "heading", "F1mean", "F2mean", "Mzmean"]
        + [f"eta{i}" for i in range(1, 7)])
IU, IV, IR, IH = 1, 2, 3, 4
IY, IN = 6, 7
IGNORE = shutil.ignore_patterns("__pycache__", "*.pre*", "processor*",
                                "postProcessing", "log", "log.*")


# ------------------------------------------------------------- sweep file ----
def read_sweep(path):
    cfg = {"lam": None, "straight": None, "out": None}
    runs = []
    for ln, line in enumerate(open(path), 1):
        tok = line.split("#", 1)[0].split()
        if not tok:
            continue
        key = tok[0]
        if key == "run":
            if len(tok) < 4:
                sys.exit(f"{path}:{ln}: 'run u v r [options]'")
            u, v, r = (float(x) for x in tok[1:4])
            if r == 0.0:
                sys.exit(f"{path}:{ln}: r = 0 is a straight course, not a circle "
                         f"(it comes from the 'straight' table)")
            runs.append({"u": u, "v": v, "r": r, "extra": tok[4:],
                         "name": f"u{u:g}_v{v:g}_r{r:g}"})
        elif key in cfg and len(tok) == 2:
            cfg[key] = tok[1]
        else:
            sys.exit(f"{path}:{ln}: cannot read '{line.strip()}'")
    if cfg["lam"] is None:
        sys.exit(f"{path}: no 'lam' line")
    cfg["lam"] = float(cfg["lam"])
    names = [r["name"] for r in runs]
    if len(set(names)) != len(names):
        sys.exit(f"{path}: the same (u, v, r) twice")
    if not runs:
        sys.exit(f"{path}: no 'run' lines")
    return cfg, runs


def estimate(lam, u, v, r, onset=10.0, ramp=2.0, Co=0.2, discX=7):
    """endTime and number of time steps as runCMT.py sets them (its defaults)."""
    k = 2*np.pi/lam
    c = np.sqrt(G/k)
    zmin = min(-0.5*lam, -7*DRAFT)
    omega_h = np.sqrt(G*k*np.tanh(-k*zmin))
    Te0 = 2*np.pi/(omega_h + k*u)                        # head seas at the onset
    al = np.linspace(0.0, 2*np.pi, 721)
    we = np.abs(-omega_h + k*(np.cos(al)*u + np.sin(al)*v))
    Te_max = 2*np.pi/max(we.min(), 0.05*omega_h)
    end = (onset + 0.5*ramp)*Te0 + 2*np.pi/abs(r) + 3*Te_max
    xbody = 7*lam
    VS = lambda x, y: np.hypot(-u + r*y, v - r*(x - xbody))
    Nref = 2 + (lam >= 1.5) + (lam >= 3.0) + (lam >= 4.5)
    l3 = max(2.0*L_PP, lam)
    Vn = max(VS(xbody + sx*l3, sy*l3) for sx in (-1, 1) for sy in (-1, 1))
    Vf = max(VS(x, y) for x in (0.0, 14*lam) for y in (-7*lam, 7*lam))
    dt = min(Co*lam/((c + Vn)*discX)/(Nref + 1), Co*lam/((c + Vf)*discX))
    return end, end/dt


# ---------------------------------------------------------------- tables ----
def read_table(path):
    """(header, MMG marker line or None, rows as an array in the header's order)"""
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    head = [h.strip() for h in lines[0].split(",")]
    marker = next((l for l in lines[1:] if l.startswith("#") and "convention: MMG" in l), None)
    rows = [[float(x) for x in l.split(",")[:len(head)]]
            for l in lines[1:] if l.strip() and not l.lstrip().startswith("#")]
    return head, marker, np.array(rows, dtype=float).reshape(-1, len(head))


def to_cols(head, D, path):
    """Rows in COLS order; V = 0 and r = 0 where the table has no such column."""
    out = np.zeros((len(D), len(COLS)))
    for j, c in enumerate(COLS):
        if c in head:
            out[:, j] = D[:, head.index(c)]
        elif c not in ("V", "r"):
            sys.exit(f"{path}: no column '{c}'")
    return out


def states(T):
    """The distinct (U, V, r) of the rows, and a function selecting a state."""
    keys = np.unique(np.round(T[:, IU:IR + 1], 6) + 0.0, axis=0)
    sel = lambda key: np.all(np.isclose(T[:, IU:IR + 1], key), axis=1)
    return keys, sel


def run_status(case, run):
    """('complete' | description, rows in COLS order or None)"""
    p = case / "results.csv"
    if not p.is_file():
        return ("not run" if not case.exists() else "no results.csv"), None
    head, marker, D = read_table(p)
    if marker is None:
        return "results.csv without the MMG marker", None
    R = to_cols(head, D, p)
    mus = R[:, IH]
    if np.isclose(mus, 0.0).any() and np.isclose(mus, 360.0).any():
        return "complete", R
    # The turn stops a hair short of the end row (head seas again) when
    # endTime has no margin: complete if every other 10 deg row is there
    end = 0.0 if run["r"] > 0 else 360.0
    rest = [m for m in np.arange(0.0, 361.0, 10.0) if not np.isclose(m, end)]
    if all(np.isclose(mus, m).any() for m in rest):
        return "complete, no end-of-turn row", R
    return (f"partial ({len(np.unique(np.round(mus)))} headings, "
            f"{mus.min():g}..{mus.max():g} deg)"), R


def straight_rows(path, lamL):
    """Straight-course rows at lam/L, completed to 0..360, gaps filled across U."""
    head, marker, D = read_table(path)
    if marker is None:
        sys.exit(f"{path}: no '# convention: MMG' line (convert it with "
                 f"manModel/applications/convertWaveTable.py first)")
    T = to_cols(head, D, path)
    T = T[np.isclose(T[:, 0], lamL, atol=1e-3)]
    if not len(T):
        sys.exit(f"{path}: no rows at lam/L {lamL:g}")
    half = T[:, IH].max() <= 180.0 + 1e-6
    if half:                                             # waves from port: mirror image
        m = T[(T[:, IH] > 1e-6) & (T[:, IH] < 180.0 - 1e-6)].copy()
        m[:, IH] = 360.0 - m[:, IH]
        m[:, [IV, IR, IY, IN]] *= -1.0
        T = np.vstack([T, m])
    # headings missing at one speed: from the nearest speed (same V, r) that has them
    filled, add = [], []
    vr =np.unique(np.round(T[:, IV:IR + 1], 6) + 0.0, axis=0)
    for V, r in vr:
        S = T[np.isclose(T[:, IV], V) & np.isclose(T[:, IR], r)]
        Us = np.unique(S[:, IU])
        mus = np.unique(np.round(S[:, IH], 6))
        for U in Us:
            have = np.round(S[np.isclose(S[:, IU], U), IH], 6)
            for mu in mus:
                if np.isclose(have, mu).any():
                    continue
                cand = [u for u in Us if np.isclose(S[np.isclose(S[:, IU], u), IH], mu).any()]
                un = cand[int(np.argmin(np.abs(np.array(cand) - U)))]
                row = S[np.isclose(S[:, IU], un) & np.isclose(S[:, IH], mu)][0].copy()
                row[IU] = U
                add.append(row)
                filled.append((U, mu, un))
    if add:
        T = np.vstack([T] + [np.array(add)])
    return T, half, filled


def close_circle(T):
    """mu = 360 as a copy of mu = 0 (or the other way round) for every (U, V, r)."""
    add = []
    keys, sel = states(T)
    for key in keys:
        S = T[sel(key)]
        for a, b in ((0.0, 360.0), (360.0, 0.0)):
            ia = np.isclose(S[:, IH], a)
            if ia.any() and not np.isclose(S[:, IH], b).any():
                row = S[ia][0].copy()
                row[IH] = b
                add.append(row)
    return np.vstack([T] + [np.array(add)]) if add else T


# ------------------------------------------------------------------ runs ----
def run_one(i, run, case, lam, passthru):
    cmd = (["python3", "runCMT.py", f"--lam={lam:g}", f"--u={run['u']:g}",
            f"--v={run['v']:g}", f"--r={run['r']:g}"] + run["extra"] + passthru)
    if not case.exists():
        shutil.copytree(BASE, case, ignore=IGNORE)
    print(f"[{i}] {run['name']}: {' '.join(cmd)}\n      in {case} "
          f"(output in runSweep.log), started {time.strftime('%Y-%m-%d %H:%M')}", flush=True)
    t0 = time.time()
    with open(case / "runSweep.log", "w") as f:
        rc = subprocess.run(cmd, cwd=case, stdout=f, stderr=subprocess.STDOUT).returncode
    st, _ = run_status(case, run)
    print(f"[{i}] {run['name']}: runCMT.py rc {rc}, {st}, "
          f"{(time.time() - t0)/3600:.2f} h", flush=True)


# --------------------------------------------------------------- collect ----
def collect(cfg, runs, sweep_dir, root, keep_onset):
    blocks, used, marker, lamL = [], [], None, None
    for i, run in enumerate(runs, 1):
        case = sweep_dir / run["name"]
        st, R = run_status(case, run)
        if not st.startswith("complete"):
            print(f"  [{i}] {run['name']}: {st} -- left out")
            continue
        head, mk, _ = read_table(case / "results.csv")
        marker = marker or mk
        if not (np.allclose(R[:, IU], run["u"]) and np.allclose(R[:, IV], run["v"])
                and np.allclose(R[:, IR], run["r"])):
            print(f"  [{i}] {run['name']}: WARNING results.csv has U, V, r "
                  f"{R[0, IU]:g}, {R[0, IV]:g}, {R[0, IR]:g}")
        if lamL is None:
            lamL = R[0, 0]
        if not np.allclose(R[:, 0], lamL, atol=1e-3):
            sys.exit(f"{case}/results.csv: lam/L {R[0, 0]:g}, the others {lamL:g}")
        onset, end = (360.0, 0.0) if run["r"] > 0 else (0.0, 360.0)
        io, ie = np.isclose(R[:, IH], onset), np.isclose(R[:, IH], end)
        if not ie.any():
            # same heading: the end row is the onset row (close_circle copies it)
            print(f"  [{i}] {run['name']}: {len(R)} headings; no end-of-turn row, "
                  f"head seas mu {end:g} = the onset row (mu {onset:g})")
            run["onsetOnly"] = True
        elif not keep_onset:
            d = R[io][0, IH + 1:] - R[ie][0, IH + 1:]
            R[io, IH + 1:] = R[ie][0, IH + 1:]
            print(f"  [{i}] {run['name']}: {len(R)} headings; onset row (mu {onset:g}) replaced "
                  f"by the end of the turn (changed X, Y, N by {d[0]:+.3f}, {d[1]:+.3f}, {d[2]:+.3f})")
        else:
            print(f"  [{i}] {run['name']}: {len(R)} headings")
        blocks.append(R)
        used.append(run)
    if not blocks:
        print("collect: no complete run yet")
        return

    contents = (f"# contents: circular motion tests at lam/L {lamL:g} (runSweep.py "
                f"{cfg['file']}), (u, v, r) = "
                + ", ".join(f"({r['u']:g}, {r['v']:g}, {r['r']:g})" for r in used))
    full = [r for r in used if not r.get("onsetOnly")]
    part = [r for r in used if r.get("onsetOnly")]
    if full and not keep_onset:
        contents += "; head seas at the yaw onset replaced by the end of the turn"
    if part:
        contents += ("; head seas = the yaw onset (no end-of-turn row) for "
                     + ", ".join(f"({r['u']:g}, {r['v']:g}, {r['r']:g})" for r in part))
    if cfg["straight"] and cfg["straight"].lower() != "none":
        sp = (root / cfg["straight"]).resolve()
        S, half, filled = straight_rows(sp, lamL)
        blocks.append(S)
        keys, _ = states(S)
        print(f"  straight course {sp}: (U, V, r) "
              + ", ".join(f"({k[0]:g}, {k[1]:g}, {k[2]:g})" for k in keys)
              + (", completed to 0..360 by symmetry" if half else ""))
        if filled:
            print("    filled from the nearest speed: "
                  + ", ".join(f"U {u:g} mu {m:g} <- U {un:g}" for u, m, un in filled))
        contents += (f"; straight course from {sp.name} ("
                     + ("half table completed to 0..360 by symmetry, " if half else "")
                     + "U " + ", ".join(f"{u:g}" for u in np.unique(S[:, IU]))
                     + (f"; {len(filled)} headings filled from the nearest speed" if filled else "")
                     + ")")
    T = close_circle(np.vstack(blocks))
    contents += "; mu 360 = mu 0"
    T = T[np.lexsort([T[:, IH], T[:, IR], T[:, IV], T[:, IU], T[:, 0]])]

    out = (root / cfg["out"]) if cfg["out"] else (sweep_dir / "waveData.dat")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(",".join(COLS) + "\n")
        f.write(marker + "\n")
        f.write(contents + "\n")
        for row in T:
            f.write(",".join(f"{x + 0.0:.6g}" for x in row) + "\n")

    keys, sel = states(T)
    Us, Vs, Rs = (np.unique(keys[:, j]) for j in range(3))
    print(f"wrote {out}: {len(T)} rows")
    print(f"  grid U {list(Us)}, V {list(Vs)}, r {list(Rs)}: {len(keys)} of "
          f"{len(Us)*len(Vs)*len(Rs)} (U, V, r) states have data (manModel fills the "
          f"others along U, then V, then r)")
    full = np.arange(0.0, 361.0, 10.0)
    for key in keys:
        mus = T[sel(key), IH]
        miss = [m for m in full if not np.isclose(mus, m).any()]
        if miss:
            print(f"  (U, V, r) = ({key[0]:g}, {key[1]:g}, {key[2]:g}): no data at mu "
                  + ", ".join(f"{m:g}" for m in miss))


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep", help="sweep file, e.g. sweep_lam07.dat")
    ap.add_argument("--dryRun", action="store_true", help="list the runs, do nothing")
    ap.add_argument("--only", nargs="+", default=None, metavar="RUN",
                    help="only these runs (numbers from --dryRun, or names)")
    ap.add_argument("--collect", action="store_true", help="only collect the results")
    ap.add_argument("--noCollect", action="store_true", help="do not collect after the runs")
    ap.add_argument("--force", action="store_true", help="rerun complete runs too")
    ap.add_argument("--keepOnset", action="store_true",
                    help="keep the transitional head-sea row at the yaw onset")
    a, passthru = ap.parse_known_args()

    sweep = Path(a.sweep).resolve()
    root = sweep.parent
    cfg, runs = read_sweep(sweep)
    cfg["file"] = sweep.name
    sweep_dir = root / sweep.stem
    lam = cfg["lam"]

    pick = list(enumerate(runs, 1))
    if a.only:
        pick = [(i, r) for i, r in pick if str(i) in a.only or r["name"] in a.only]
        unknown = set(a.only) - {str(i) for i, _ in pick} - {r["name"] for _, r in pick}
        if unknown:
            sys.exit(f"--only: no run {', '.join(sorted(unknown))}")

    print(f"sweep {sweep.name}: lambda {lam:g} m (lam/L {lam/L_PP:.3g}), {len(runs)} runs in {sweep_dir}")
    if passthru:
        print(f"  passed on to runCMT.py: {' '.join(passthru)}")
    if not a.collect or a.dryRun:
        tot = 0.0
        for i, run in pick:
            end, steps = estimate(lam, run["u"], run["v"], run["r"])
            tot += steps
            st, _ = run_status(sweep_dir / run["name"], run)
            beta = np.degrees(np.arctan2(-run["v"], run["u"]))
            print(f"  [{i}] {run['name']:22s} beta {beta:5.1f} deg, radius "
                  f"{np.hypot(run['u'], run['v'])/abs(run['r'])/L_PP:4.2f} L, "
                  f"endTime ~{end:6.1f} s, ~{steps/1e3:5.1f}k steps  [{st}]"
                  + (f"  options {' '.join(run['extra'])}" if run["extra"] else ""))
        print(f"  total ~{tot/1e3:.0f}k time steps (runCMT.py defaults; same mesh for every run)")
    if a.dryRun:
        return

    if not a.collect:
        if not BASE.is_dir():
            sys.exit(f"no {BASE}")
        sweep_dir.mkdir(exist_ok=True)
        for i, run in pick:
            case = sweep_dir / run["name"]
            st, _ = run_status(case, run)
            if st.startswith("complete") and not a.force:
                print(f"[{i}] {run['name']}: complete, skipped")
                continue
            if case.exists():
                print(f"[{i}] {run['name']}: {st} -- running again in place")
            run_one(i, run, case, lam, passthru)
        if a.noCollect:
            return

    print("collecting")
    collect(cfg, runs, sweep_dir, root, a.keepOnset)


if __name__ == "__main__":
    main()
