#!/usr/bin/env python3
"""Run production Sod or four-quadrant Riemann validation matrices."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from run_release_matrix import clean_work_directory, one_field, render, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=("sod", "quadrant"))
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--cells-i", type=int, default=80)
    parser.add_argument("--cells-j", type=int, default=8)
    parser.add_argument("--zones-i", type=int, default=1)
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--reconstruction", default="weno_z")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--end-time", type=float, default=0.2)
    parser.add_argument("--cfl", type=float, default=0.1)
    parser.add_argument("--min-cells", type=int, default=8)
    parser.add_argument("--density-l1", type=float, default=2.5e-2)
    parser.add_argument("--position-cells", type=float, default=3.0)
    parser.add_argument("--symmetry-l1", type=float, default=5.0e-4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if any(value < 1 for value in ranks) or (
        any(value > 1 for value in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("shock rank configuration is invalid")
    if min(args.cells_i, args.cells_j, args.zones_i, args.min_cells) < 1:
        raise RuntimeError("shock grid parameters must be positive")
    if args.cells_i % args.zones_i != 0:
        raise RuntimeError("shock cells_i must be divisible by zones_i")
    if min(args.end_time, args.cfl, args.density_l1, args.position_cells, args.symmetry_l1) <= 0.0:
        raise RuntimeError("shock times and tolerances must be positive")

    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / f"{args.case}.cgns"
    length_y = 0.1 if args.case == "sod" else 1.0
    records: list[dict[str, object]] = []
    records.append(
        run(
            [
                str(args.generator),
                "rectangle",
                str(mesh),
                str(args.cells_i),
                str(args.cells_j),
                str(args.zones_i),
                "1.0",
                str(length_y),
                "false",
            ],
            root / "generate.log",
        )
    )
    template = args.template.read_text(encoding="utf-8")
    fields: dict[int, Path] = {}
    for rank_count in ranks:
        case_name = (
            f"{args.case}-{args.profile}-{args.reconstruction}-"
            f"{args.riemann}-r{rank_count}"
        )
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
                [str(args.validator), "finite", str(fields[rank_count])],
                root / f"finite-r{rank_count}.log",
            )
        )
        if args.case == "sod":
            validation = [
                str(args.validator),
                "sod",
                str(fields[rank_count]),
                str(args.end_time),
                "0.5",
                "1.4",
                str(args.density_l1),
                str(args.position_cells),
            ]
        else:
            validation = [
                str(args.validator),
                "diagonal-symmetry",
                str(fields[rank_count]),
                str(args.symmetry_l1),
            ]
        records.append(run(validation, root / f"analytic-r{rank_count}.log"))

    for rank_count in ranks[1:]:
        records.append(
            run(
                [
                    str(args.validator),
                    "compare",
                    str(fields[ranks[0]]),
                    str(fields[rank_count]),
                    "2e-11",
                ],
                root / f"compare-r{ranks[0]}-r{rank_count}.log",
            )
        )
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "case": args.case,
        "grid": [args.cells_i, args.cells_j],
        "zones_i": args.zones_i,
        "profile": args.profile,
        "reconstruction": args.reconstruction,
        "riemann": args.riemann,
        "end_time": args.end_time,
        "cfl": args.cfl,
        "ranks": ranks,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"shock matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"shock matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)
