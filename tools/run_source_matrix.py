#!/usr/bin/env python3
"""Run the uniform conservative-source production validation matrix."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from run_release_matrix import clean_work_directory, one_field, one_file, render, run


INITIAL = (1.0, 0.1, -0.05, 0.0, 2.50625)
SOURCE = (0.01, 0.02, -0.01, 0.0, 0.03)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--resolution", type=int, default=16)
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--reconstruction", default="weno_z")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--end-time", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.1)
    parser.add_argument("--tolerance", type=float, default=2.0e-12)
    parser.add_argument("--min-cells", type=int, default=4)
    return parser.parse_args()


def validate_series(path: Path, end_time: float, tolerance: float) -> float:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if len(rows) < 2 or any(len(row) != 6 for row in rows):
        raise RuntimeError(f"uniform-source statistics have invalid rows: {path}")
    selected_initial = (INITIAL[0], INITIAL[1], INITIAL[2], INITIAL[4])
    selected_source = (SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[4])
    maximum = 0.0
    for row in rows:
        if not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"uniform-source statistics are non-finite: {path}")
        time = row[1]
        for actual, initial, source in zip(row[2:], selected_initial, selected_source):
            maximum = max(maximum, abs(actual - (initial + time * source)))
    if abs(rows[-1][1] - end_time) > 1.0e-14 or maximum > tolerance:
        raise RuntimeError(
            f"uniform-source statistics error {maximum} exceeds {tolerance}: {path}"
        )
    return maximum


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if any(value < 1 for value in ranks) or (
        any(value > 1 for value in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("source rank configuration is invalid")
    if min(args.resolution, args.min_cells) < 1:
        raise RuntimeError("source grid parameters must be positive")
    if min(args.end_time, args.cfl, args.tolerance) <= 0.0:
        raise RuntimeError("source time and tolerances must be positive")

    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / "periodic-square.cgns"
    records: list[dict[str, object]] = []
    records.append(
        run(
            [
                str(args.generator), "periodic-square", str(mesh),
                str(args.resolution), str(args.resolution), "1.0",
            ],
            root / "generate.log",
        )
    )
    template = args.template.read_text(encoding="utf-8")
    fields: dict[int, Path] = {}
    for rank_count in ranks:
        case_name = f"uniform-source-{args.profile}-r{rank_count}"
        output = root / f"output-r{rank_count}"
        config = root / f"{case_name}.wcns"
        config.write_text(
            render(
                template,
                {
                    "CASE_NAME": case_name,
                    "MESH_PATH": mesh.as_posix(),
                    "ALGORITHM_PROFILE": args.profile,
                    "RECONSTRUCTION": args.reconstruction,
                    "RIEMANN": args.riemann,
                    "MIN_CELLS": str(args.min_cells),
                    "CFL": str(args.cfl),
                    "END_TIME": str(args.end_time),
                    "OUTPUT_DIRECTORY": output.as_posix(),
                },
            ),
            encoding="utf-8",
        )
        command = [str(args.run), "--config", str(config)]
        if rank_count > 1:
            command = [str(args.mpi_exec), "-n", str(rank_count)] + command
        records.append(run(command, root / f"run-r{rank_count}.log"))
        fields[rank_count] = one_field(output)
        records.append(
            run(
                [
                    str(args.validator), "uniform-source", str(fields[rank_count]),
                    str(args.end_time),
                    *(str(value) for value in INITIAL),
                    *(str(value) for value in SOURCE),
                    str(args.tolerance),
                ],
                root / f"field-r{rank_count}.log",
            )
        )
        series = one_file(output, f"*.statistics.r{rank_count}.txt")
        maximum = validate_series(series, args.end_time, args.tolerance)
        series_log = root / f"series-r{rank_count}.log"
        series_log.write_text(
            f"check=uniform_source_series max_abs={maximum:.17g} "
            f"tolerance={args.tolerance:.17g}\n",
            encoding="utf-8",
        )
        records.append({"check": "uniform_source_series", "log": str(series_log)})
    for rank_count in ranks[1:]:
        records.append(
            run(
                [
                    str(args.validator), "compare", str(fields[ranks[0]]),
                    str(fields[rank_count]), "2e-11",
                ],
                root / f"compare-r{ranks[0]}-r{rank_count}.log",
            )
        )
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "profile": args.profile,
        "resolution": args.resolution,
        "end_time": args.end_time,
        "ranks": ranks,
        "initial": INITIAL,
        "source": SOURCE,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"source matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"source matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)
