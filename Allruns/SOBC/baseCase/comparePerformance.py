#!/usr/bin/env python3
"""Compare manFlow execution logs without running or editing either case."""
import argparse
import math
from pathlib import Path
import re
import statistics

NUMBER = r"[\d.eE+-]+"


def read_steps(path):
    text = path.read_text()
    rows = []
    for block in re.split(r"^Time = ", text, flags=re.MULTILINE)[1:]:
        step = re.match(rf"({NUMBER})\s+deltaT = ({NUMBER})", block)
        timing = re.search(
            rf"ExecutionTime = ({NUMBER}) s\s+ClockTime = ({NUMBER}) s", block
        )
        solves = re.findall(
            r"Solving for PhiD,.*?No Iterations (\d+)", block
        )
        if step and timing:
            rows.append(dict(t=float(step[1]), dt=float(step[2]),
                             cpu=float(timing[1]), wall=float(timing[2]),
                             iterations=[int(n) for n in solves]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", type=Path, nargs="+", help="case directories or log files")
    parser.add_argument("--from-time", type=float, default=20.0)
    parser.add_argument("--to-time", type=float, default=30.0)
    args = parser.parse_args()
    if not (math.isfinite(args.from_time) and math.isfinite(args.to_time)
            and args.from_time < args.to_time):
        parser.error("require finite from-time < to-time")

    for case in args.cases:
        path = case / "log" if case.is_dir() else case
        try:
            rows = read_steps(path)
        except OSError as exc:
            parser.error(str(exc))
        if any(b["t"] <= a["t"] or b["cpu"] < a["cpu"] or b["wall"] < a["wall"]
               for a, b in zip(rows, rows[1:])):
            parser.error(f"{path}: restarted/concatenated log; supply one continuous run")
        window = [r for r in rows if args.from_time <= r["t"] <= args.to_time]
        if len(window) < 2:
            parser.error(f"{path}: need at least two completed steps in the chosen window")
        # Subtract the first timestamp: neither initialization nor that first
        # step belongs to the measured interval.
        measured = window[1:]
        cpu = window[-1]["cpu"] - window[0]["cpu"]
        wall = window[-1]["wall"] - window[0]["wall"]
        simulated = window[-1]["t"] - window[0]["t"]
        print(f"\n{case}")
        print(f"  Window: {window[0]['t']:g} .. {window[-1]['t']:g} s; {len(measured)} steps")
        print(f"  CPU / step: {cpu / len(measured):.6f} s")
        print(f"  Wall / step: {wall / len(measured):.6f} s")
        print(f"  Wall / simulated second: {wall / simulated:.4f} s")
        print("  dt in whole log:", ", ".join(f"{v:g}" for v in sorted({r['dt'] for r in rows})))
        print("  dt in measured window:", ", ".join(f"{v:g}" for v in sorted({r['dt'] for r in measured})))
        counts = {len(r["iterations"]) for r in measured}
        if len(counts) == 1 and next(iter(counts)):
            n = next(iter(counts))
            means = [statistics.mean(r["iterations"][j] for r in measured) for j in range(n)]
            print("  Mean iterations by solve:", ", ".join(f"{v:.3f}" for v in means))
        else:
            print("  Solve counts vary or iteration logging is unavailable:", sorted(counts))
        print(f"  Whole-log final wall time (includes startup): {rows[-1]['wall']:g} s")
    print("\nCompare identical time windows, dt, mesh, ranks, output settings and hardware load.")
    print("Repeat timings; fewer iterations alone do not establish a speedup.")


if __name__ == "__main__":
    main()
