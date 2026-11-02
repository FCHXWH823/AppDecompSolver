#!/usr/bin/env python3
"""Run the six 9-input benchmarks used by Table 1.

The script intentionally records raw stdout, parsed MED/runtime, command lines,
and failures. Long-running Gurobi methods are opt-in.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = {
    "cos": "cos-i-9-o-10-a-0-b-1_57079632-acc.txt",
    "tan": "tan-i-9-o-10-a-0-b-1_256637061-acc.txt",
    "exp": "exp-i-9-o-10-a-0-b-3-acc.txt",
    "ln": "ln-i-9-o-10-a-1-b-10-acc.txt",
    "erf": "erf-i-9-o-10-a-0-b-3-acc.txt",
    "denoise": "denoise-i-9-o-10-a-0-b-3-acc.txt",
}
PAPER = {
    "dalta_joint": {
        "cos": [2.96, 3.06], "tan": [3.24, 2.83], "exp": [4.22, 2.72],
        "ln": [4.69, 6.77], "erf": [1.85, 2.76], "denoise": [4.75, 2.81],
    },
    "dalta_joint_ilp": {
        "cos": [2.48, 3600.0], "tan": [2.62, 3600.0],
        "exp": [3.55, 3600.0], "ln": [2.55, 3600.0],
        "erf": [2.66, 3600.0], "denoise": [3.38, 3600.0],
    },
    "dalta_separate_ilp": {
        "cos": [11.64, 258.37], "tan": [10.91, 236.32],
        "exp": [9.26, 242.58], "ln": [8.32, 224.68],
        "erf": [5.07, 139.60], "denoise": [10.91, 229.25],
    },
    "ising_separate": {
        "cos": [8.33, 0.56], "tan": [10.45, 0.56], "exp": [7.07, 0.74],
        "ln": [6.57, 0.49], "erf": [4.61, 0.42], "denoise": [9.69, 0.46],
    },
    "ising_joint": {
        "cos": [2.50, 1.75], "tan": [2.50, 1.87], "exp": [2.66, 1.92],
        "ln": [2.72, 2.77], "erf": [1.90, 1.55], "denoise": [2.80, 1.51],
    },
}

def command_for(
    method: str, benchmark: Path, seed: int, mip_time: float, mip_gap: float,
) -> tuple[list[str], Path, dict[str, str]]:
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    if method.startswith("ising_"):
        paper_plus = method.endswith("_paper_plus")
        paper = method.endswith("_paper")
        if not paper and not paper_plus:
            raise ValueError(
                "Ising method must end in _paper or _paper_plus"
            )
        mode = method.removeprefix("ising_")
        mode = mode.removesuffix("_paper_plus")
        mode = mode.removesuffix("_paper")
        command = [
            "build/ising-decomp",
            "--input", str(benchmark), "--mode", mode,
            "--n-free", "4", "--rounds", "5" if mode == "joint" else "1",
            "--partitions", "1000", "--seed", str(seed),
            "--max-steps", "2000", "--quiet",
        ]
        if paper_plus:
            command.extend(["--preset", "paper-plus"])
        else:
            command.extend(["--preset", "paper"])
        return command, ROOT, env

    solver = {
        "dalta_joint": "joint-heuristic",
        "dalta_joint_ilp": "joint-ilp",
        "dalta_separate_legacy": "separate-legacy",
        "dalta_separate_ilp": "separate-ilp",
    }[method]
    command = [
        "build/dalta-baseline", "--input", str(benchmark),
        "--output", "results/table1/artifacts",
        "--nFree", "4", "--startBit", "8", "--nTry", "5",
        "--nComb", "1000", "--admmReInitTime", "30",
        "--admmStuckTol", "3", "--seed", str(seed), "--solver", solver,
        "--MIPTime", str(mip_time), "--MIPGap", str(mip_gap),
    ]
    return command, ROOT, env


def parse_result(method: str, output: str) -> tuple[float, float]:
    if method.startswith("ising_"):
        match = re.search(r"RESULT .*?med=([0-9.eE+-]+) time_s=([0-9.eE+-]+)", output)
        if not match:
            raise ValueError("missing Ising RESULT line")
        return float(match.group(1)), float(match.group(2))
    med_values = re.findall(r"^med = ([0-9.eE+-]+)$", output, re.MULTILINE)
    time_values = re.findall(r"^time = ([0-9.eE+-]+) s$", output, re.MULTILINE)
    if not med_values or not time_values:
        raise ValueError("missing DALTA MED/runtime lines")
    return float(med_values[-1]), float(time_values[-1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--methods", nargs="+",
        default=[
            "ising_separate_paper", "ising_joint_paper",
            "ising_separate_paper_plus", "ising_joint_paper_plus",
        ],
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--benchmarks", nargs="+", choices=BENCHMARKS, default=list(BENCHMARKS))
    parser.add_argument("--mip-time", type=float, default=3600.0,
                        help="Gurobi time limit for each model")
    parser.add_argument("--mip-gap", type=float, default=0.01,
                        help="Gurobi relative MIP gap")
    parser.add_argument("--output", type=Path, default=ROOT / "results/table1/results.json")
    args = parser.parse_args()

    if any(not method.startswith("ising_") for method in args.methods):
        (ROOT / "results/table1/artifacts").mkdir(parents=True, exist_ok=True)

    records: dict[str, object] = {"paper": PAPER, "seed": args.seed, "runs": {}}
    for method in args.methods:
        method_records: dict[str, object] = {}
        records["runs"][method] = method_records  # type: ignore[index]
        for name in args.benchmarks:
            filename = BENCHMARKS[name]
            benchmark = Path("benchmarks/numerics") / filename
            command, cwd, env = command_for(
                method, benchmark, args.seed, args.mip_time, args.mip_gap,
            )
            started = time.time()
            print(f"[{method}] {name}", flush=True)
            try:
                process = subprocess.run(
                    command, cwd=cwd, env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=args.timeout, check=False,
                )
                med, runtime = parse_result(method, process.stdout)
                method_records[name] = {
                    "status": "ok" if process.returncode == 0 else "failed",
                    "med": med, "reported_time_s": runtime,
                    "wall_time_s": time.time() - started,
                    "returncode": process.returncode, "command": command,
                    "stdout": process.stdout,
                }
            except subprocess.TimeoutExpired as error:
                method_records[name] = {
                    "status": "timeout", "wall_time_s": time.time() - started,
                    "command": command, "stdout": error.stdout or "",
                }
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(records, indent=2), encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
