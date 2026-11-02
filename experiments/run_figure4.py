#!/usr/bin/env python3
"""Run and incrementally record all ten 16-input Figure 4 benchmarks."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT_RE = re.compile(
    r"RESULT mode=joint med=([0-9.eE+-]+) time_s=([0-9.eE+-]+)"
)
CASES = {
    "cos": {
        "name": "cos(x)",
        "input": "benchmarks/numerics/cos-i-16-o-16-a-0-b-1_57079632-acc.txt",
        "dalta_med": 10.18,
        "dalta_time_s": 4477.86,
    },
    "tan": {
        "name": "tan(x)",
        "input": "benchmarks/numerics/tan-i-16-o-16-a-0-b-1_256637061-acc.txt",
        "dalta_med": 3.14,
        "dalta_time_s": 3389.72,
    },
    "exp": {
        "name": "exp(x)",
        "input": "benchmarks/numerics/exp-i-16-o-16-a-0-b-3-acc.txt",
        "dalta_med": 10.30,
        "dalta_time_s": 3865.39,
    },
    "ln": {
        "name": "ln(x)",
        "input": "benchmarks/numerics/ln-i-16-o-16-a-1-b-10-acc.txt",
        "dalta_med": 10.65,
        "dalta_time_s": 5312.16,
    },
    "erf": {
        "name": "erf(x)",
        "input": "benchmarks/numerics/erf-i-16-o-16-a-0-b-3-acc.txt",
        "dalta_med": 15.01,
        "dalta_time_s": 3983.88,
    },
    "denoise": {
        "name": "denoise(x)",
        "input": "benchmarks/numerics/denoise-i-16-o-16-a-0-b-3-acc.txt",
        "dalta_med": 33.94,
        "dalta_time_s": 4503.74,
    },
    "brent-kung": {
        "name": "Brent-Kung",
        "input": "benchmarks/axbench-verilog/adder_8.txt",
        "dalta_med": 1.11,
        "dalta_time_s": 1886.98,
    },
    "forwardk2j": {
        "name": "Forwardk2j",
        "input": "benchmarks/axbench-verilog/kin_16.txt",
        "dalta_med": 1168.18,
        "dalta_time_s": 8261.07,
    },
    "inversek2j": {
        "name": "Inversek2j",
        "input": "benchmarks/axbench-verilog/kin_inv.txt",
        "dalta_med": 379.68,
        "dalta_time_s": 6970.90,
    },
    "multiplier": {
        "name": "Multiplier",
        "input": "benchmarks/axbench-verilog/mult_8.txt",
        "dalta_med": 481.47,
        "dalta_time_s": 8506.03,
    },
}


def command_for(case: dict[str, object], seed: int) -> list[str]:
    return [
        "build/ising-decomp",
        "--input", str(case["input"]),
        "--mode", "joint",
        "--n-free", "7",
        "--rounds", "5",
        "--partitions", "1000",
        "--seed", str(seed),
        "--max-steps", "2000",
        "--preset", "paper",
        "--quiet",
    ]


def run_case(slug: str, seed: int, timeout: float) -> tuple[str, dict[str, object]]:
    case = CASES[slug]
    command = command_for(case, seed)
    env = os.environ.copy()
    env.update({
        "LC_ALL": "C",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    })
    started = time.time()
    record: dict[str, object] = {
        "name": case["name"],
        "input": case["input"],
        "command": command,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "figure4_dalta": {
            "med": case["dalta_med"],
            "time_s": case["dalta_time_s"],
        },
    }
    try:
        process = subprocess.run(
            command,
            cwd=ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
        record["wall_time_s"] = time.time() - started
        record["returncode"] = process.returncode
        record["stdout"] = process.stdout
        match = RESULT_RE.search(process.stdout)
        if process.returncode != 0 or not match:
            record["status"] = "failed"
        else:
            med = float(match.group(1))
            solver_time = float(match.group(2))
            record["status"] = "ok"
            record["result"] = {"med": med, "solver_time_s": solver_time}
            record["ratio_to_figure4_dalta"] = {
                "med": med / float(case["dalta_med"]),
                "time": solver_time / float(case["dalta_time_s"]),
                "speedup": float(case["dalta_time_s"]) / solver_time,
            }
    except subprocess.TimeoutExpired as error:
        record.update({
            "status": "timeout",
            "wall_time_s": time.time() - started,
            "stdout": error.stdout or "",
        })
    return slug, record


def import_brent_kung(records: dict[str, object]) -> None:
    path = ROOT / "results/figure4/brent-kung-paper.json"
    if not path.exists() or "brent-kung" in records["runs"]:  # type: ignore[operator]
        return
    source = json.loads(path.read_text(encoding="utf-8"))
    if source.get("configuration", {}).get("partition_strategy") != "dalta-draw":
        return
    records["runs"]["brent-kung"] = {  # type: ignore[index]
        "name": source["benchmark"],
        "input": source["input"],
        "command": source["command"],
        "status": "ok",
        "wall_time_s": source["result"]["solver_time_s"],
        "result": source["result"],
        "figure4_dalta": source["figure4_dalta"],
        "ratio_to_figure4_dalta": source["ratio_to_figure4_dalta"],
        "source": "imported from brent-kung-paper.json",
    }


def update_summary(records: dict[str, object]) -> None:
    successful = [
        record for record in records["runs"].values()  # type: ignore[union-attr]
        if record.get("status") == "ok"
    ]
    summary: dict[str, object] = {
        "completed": len(successful),
        "total": len(CASES),
    }
    if successful:
        med_ratios = [record["ratio_to_figure4_dalta"]["med"] for record in successful]
        time_ratios = [record["ratio_to_figure4_dalta"]["time"] for record in successful]
        solver_meds = [record["result"]["med"] for record in successful]
        solver_times = [record["result"]["solver_time_s"] for record in successful]
        dalta_meds = [record["figure4_dalta"]["med"] for record in successful]
        dalta_times = [record["figure4_dalta"]["time_s"] for record in successful]
        summary.update({
            "mean_case_med_ratio": sum(med_ratios) / len(med_ratios),
            "mean_case_time_ratio": sum(time_ratios) / len(time_ratios),
            "speedup_from_mean_case_time_ratio": len(time_ratios) / sum(time_ratios),
            "aggregate_med_ratio": sum(solver_meds) / sum(dalta_meds),
            "aggregate_time_ratio": sum(solver_times) / sum(dalta_times),
        })
    records["summary"] = summary


def write_records(path: Path, records: dict[str, object]) -> None:
    update_summary(records)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=14400.0)
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results/figure4/results.json",
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if not (ROOT / "build/ising-decomp").exists():
        parser.error("build/ising-decomp is missing; build the solver first")
    for slug in args.cases:
        if not (ROOT / str(CASES[slug]["input"])).exists():
            parser.error(f"missing input for {slug}: {CASES[slug]['input']}")

    if args.output.exists():
        records = json.loads(args.output.read_text(encoding="utf-8"))
    else:
        records = {
            "configuration": {
                "preset": "paper",
                "mode": "joint",
                "inputs": 16,
                "n_free": 7,
                "n_bound": 9,
                "partitions": 1000,
                "partition_strategy": "dalta-draw",
                "rounds": 5,
                "seed": args.seed,
                "max_steps": 2000,
                "sample_period": 10,
                "energy_window": 10,
                "epsilon": 1e-8,
            },
            "runs": {},
        }
    import_brent_kung(records)
    write_records(args.output, records)

    pending = [
        slug for slug in args.cases
        if args.rerun or records["runs"].get(slug, {}).get("status") != "ok"
    ]
    skipped = len(args.cases) - len(pending)
    print(f"Figure 4: pending={len(pending)} skipped={skipped} jobs={args.jobs}", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_case, slug, args.seed, args.timeout): slug
            for slug in pending
        }
        for future in concurrent.futures.as_completed(futures):
            slug, record = future.result()
            records["runs"][slug] = record
            write_records(args.output, records)
            if record["status"] == "ok":
                result = record["result"]
                print(
                    f"[{slug}] med={result['med']:.10g} "
                    f"time_s={result['solver_time_s']:.3f}",
                    flush=True,
                )
            else:
                print(f"[{slug}] status={record['status']}", flush=True)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
