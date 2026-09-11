#!/usr/bin/env python3
"""Run the four-level SSPRK3 self-convergence gate for stage S."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path


def checked(command: list[str], log: Path) -> str:
    result = subprocess.run(command, capture_output=True, text=True)
    text = result.stdout + result.stderr
    log.write_text(text, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}")
    return text


def density_l2(text: str) -> float:
    match = re.search(r"field=Density\s+samples=\d+\s+l1=\S+\s+l2=(\S+)", text)
    if match is None:
        raise RuntimeError("field-error output has no Density L2")
    value = float(match.group(1))
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError("Density L2 difference is not finite and positive")
    return value


def representative_dt(log: Path) -> tuple[int, float]:
    text = log.read_text(encoding="utf-8")
    samples = [
        float(value) for value in re.findall(r"^step=\d+\s+time=\S+\s+dt=(\S+)", text, re.MULTILINE)
    ]
    stopped = re.search(r"reason=physical_time_reached step=(\d+) time=", text)
    if stopped is None or len(samples) < 2:
        raise RuntimeError(f"time run did not provide enough accepted steps: {log}")
    # The last step is normally clipped to t_end and is not the nominal CFL step.
    return int(stopped.group(1)), statistics.median(samples[:-1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--resolution", type=int, default=48)
    parser.add_argument("--end-time", type=float, default=0.1)
    parser.add_argument("--cfl-values", default="0.4,0.2,0.1,0.05")
    parser.add_argument("--minimum-order", type=float, default=2.8)
    args = parser.parse_args()
    cfl_values = [float(value) for value in args.cfl_values.split(",")]
    if len(cfl_values) < 4 or any(value <= 0.0 for value in cfl_values):
        raise RuntimeError("stage S time convergence requires four positive CFL values")
    for coarse, fine in zip(cfl_values, cfl_values[1:]):
        if abs(coarse / fine - 2.0) > 1.0e-12:
            raise RuntimeError("stage S CFL values must form a factor-two refinement")
    root = args.work_dir.resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    vortex_matrix = Path(__file__).with_name("run_vortex_matrix.py")
    runs: list[dict[str, object]] = []
    fields: list[Path] = []
    for index, cfl in enumerate(cfl_values):
        run_root = root / f"level-{index}"
        command = [
            sys.executable,
            str(vortex_matrix),
            "--run",
            str(args.run),
            "--generator",
            str(args.generator),
            "--validator",
            str(args.validator),
            "--template",
            str(args.template),
            "--work-dir",
            str(run_root),
            "--ranks",
            "1",
            "--resolutions",
            str(args.resolution),
            "--profile",
            args.profile,
            "--end-time",
            str(args.end_time),
            "--cfl",
            str(cfl),
            "--finest-l1",
            "1",
        ]
        root_log = root / f"level-{index}.log"
        checked(command, root_log)
        case_root = run_root / f"n{args.resolution}"
        field_matches = sorted((case_root / "output-r1").glob("*.field.*.cgns"))
        if len(field_matches) != 1:
            raise RuntimeError(f"level {index} did not produce exactly one final field")
        fields.append(field_matches[0])
        steps, dt = representative_dt(case_root / "run-r1.log")
        runs.append(
            {
                "level": index,
                "cfl": cfl,
                "steps": steps,
                "representative_dt": dt,
                "field": str(fields[-1]),
            }
        )
    differences: list[float] = []
    for index in range(len(fields) - 1):
        log = root / f"difference-{index}-{index + 1}.log"
        text = checked(
            [str(args.validator), "field-error", str(fields[index + 1]), str(fields[index])],
            log,
        )
        differences.append(density_l2(text))
    orders: list[float] = []
    for coarse, fine in zip(differences, differences[1:]):
        order = math.log(coarse / fine) / math.log(2.0)
        orders.append(order)
        if order < args.minimum_order:
            raise RuntimeError(
                f"SSPRK3 self-convergence order {order} is below {args.minimum_order}"
            )
    summary = {
        "matrix_version": 1,
        "stage": "S",
        "status": "passed",
        "method": "SSPRK3 factor-two self-convergence",
        "profile": args.profile,
        "resolution": args.resolution,
        "end_time": args.end_time,
        "runs": runs,
        "density_l2_successive_differences": differences,
        "observed_orders": orders,
        "minimum_order": args.minimum_order,
    }
    path = root / "matrix-summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"stage S SSPRK3 convergence passed: {path}; orders={orders}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"stage S time convergence failed: {error}", file=sys.stderr)
        raise SystemExit(1)
