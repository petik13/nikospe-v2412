import shutil
import subprocess
from pathlib import Path

BASE = Path("baseCase")
RUNS = Path("runs")
RUNS.mkdir(exist_ok=True)
CASE_SCRIPT = "runsim.py"   

def run_one_case(lam: float, Ucur: float, heading: float) -> int:
    case_name = f"lam{lam:g}_U{Ucur:g}_head{heading:g}"
    case_dir = RUNS / case_name

    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(BASE, case_dir)

    # run the per-case script inside the case directory
    cmd = [
        "python3", CASE_SCRIPT,
        "--lam", str(lam),
        "--Ucur", str(Ucur),
        "--heading", str(heading),
    ]

    with open(case_dir / "log.driver", "w") as log:
        p = subprocess.Popen(cmd, cwd=case_dir, stdout=log, stderr=subprocess.STDOUT)
        rc = p.wait()

    return rc

def main():
    # define your sweep grid
    lams = [7.0, 6.0, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0]
    Ucurs = [0.0]
    headings = [60.0]

    for lam in lams:
        for Ucur in Ucurs:
            for heading in headings:
                print(f"=== Running lam={lam}, Ucur={Ucur}, heading={heading} ===")
                rc = run_one_case(lam, Ucur, heading)
                print(f"=== Done rc={rc} ===")

                if rc != 0:
                    print("Stopping sweep because this case failed.")
                    return

if __name__ == "__main__":
    main()