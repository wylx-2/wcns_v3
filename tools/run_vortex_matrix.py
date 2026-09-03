#!/usr/bin/env python3
"""Run the production isentropic-vortex accuracy and MPI matrix."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

from run_release_matrix import clean_work_directory, one_field, one_file, render, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--resolutions", default="16")
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--reconstruction", default="weno_z")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--end-time", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.1)
    parser.add_argument("--minimum-order", type=float, default=0.0)
    parser.add_argument("--finest-l1", type=float, default=1.0)
    parser.add_argument("--min-cells", type=int, default=8)
    parser.add_argument("--reference-velocity", type=float, default=340.0)
    parser.add_argument("--reference-temperature", type=float, default=288.15)
    parser.add_argument("--molar-mass", type=float, default=0.029)
    return parser.parse_args()


def vortex_error(log_path: Path) -> float:
    match = re.search(r"\brho_l1=([^ ]+)", log_path.read_text(encoding="utf-8"))
    if match is None:
        raise RuntimeError(f"vortex validator did not report rho_l1: {log_path}")
    return float(match.group(1))


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    resolutions = [int(value) for value in args.resolutions.split(",")]
    if any(value < 1 for value in ranks) or any(value < 16 for value in resolutions):
        raise RuntimeError("vortex ranks/resolutions are invalid")
    if any(value > 1 for value in ranks) and args.mpi_exec is None:
        raise RuntimeError("MPI vortex ranks require --mpi-exec")
    if args.end_time <= 0.0 or args.cfl <= 0.0 or args.finest_l1 <= 0.0:
        raise RuntimeError("vortex time, CFL and finest tolerance must be positive")
    if min(args.reference_velocity, args.reference_temperature, args.molar_mass) <= 0.0:
        raise RuntimeError("vortex reference scales must be positive")
    gamma = 1.4
    mach = args.reference_velocity / math.sqrt(
        gamma * (8.314 / args.molar_mass) * args.reference_temperature
    )
    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    template = args.template.read_text(encoding="utf-8")
    records: list[dict[str, object]] = []
    errors: dict[int, float] = {}
    for resolution in resolutions:
        resolution_root = root / f"n{resolution}"
        resolution_root.mkdir()
        mesh = resolution_root / "periodic-square.cgns"
        records.append(
            run(
                [
                    str(args.generator),
                    "periodic-square",
                    str(mesh),
                    str(resolution),
                    str(resolution),
                    "10.0",
                ],
                resolution_root / "generate.log",
            )
        )
        fields: dict[int, Path] = {}
        for rank_count in ranks:
            case_name = f"vortex-{args.profile}-n{resolution}-r{rank_count}"
            output = resolution_root / f"output-r{rank_count}"
            config = resolution_root / f"{case_name}.wcns"
            config.write_text(
                render(
                    template,
                    {
                        "CASE_NAME": case_name,
                        "MESH_PATH": mesh.as_posix(),
                        "ALGORITHM_PROFILE": args.profile,
                        "RECONSTRUCTION": args.reconstruction,
                        "RIEMANN": args.riemann,
                        "REFERENCE_VELOCITY": str(args.reference_velocity),
                        "REFERENCE_TEMPERATURE": str(args.reference_temperature),
                        "MOLAR_MASS": str(args.molar_mass),
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
            records.append(
                run(command, resolution_root / f"run-r{rank_count}.log")
            )
            fields[rank_count] = one_field(output)
            validation_log = resolution_root / f"vortex-r{rank_count}.log"
            records.append(
                run(
                    [
                        str(args.validator),
                        "vortex",
                        str(fields[rank_count]),
                        str(args.end_time),
                        "10.0",
                        "5.0",
                        "5.0",
                        "5.0",
                        "1.0",
                        "1.0",
                        str(gamma),
                        str(mach),
                        str(max(1.0, args.finest_l1)),
                    ],
                    validation_log,
                )
            )
            records.append(
                run(
                    [
                        str(args.validator),
                        "series-constant",
                        str(one_file(output, f"*.statistics.r{rank_count}.txt")),
                        "2e-11",
                    ],
                    resolution_root / f"conservation-r{rank_count}.log",
                )
            )
            if rank_count == ranks[0]:
                errors[resolution] = vortex_error(validation_log)
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
                    resolution_root / f"compare-r{ranks[0]}-r{rank_count}.log",
                )
            )
    ordered = sorted(errors)
    orders: list[dict[str, float]] = []
    for coarse, fine in zip(ordered, ordered[1:]):
        order = math.log(errors[coarse] / errors[fine]) / math.log(fine / coarse)
        orders.append({"coarse": coarse, "fine": fine, "rho_l1_order": order})
        if order < args.minimum_order:
            raise RuntimeError(
                f"vortex order {order} for {coarse}->{fine} is below {args.minimum_order}"
            )
    if errors[ordered[-1]] > args.finest_l1:
        raise RuntimeError("vortex finest-grid L1 exceeds tolerance")
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "profile": args.profile,
        "reconstruction": args.reconstruction,
        "riemann": args.riemann,
        "end_time": args.end_time,
        "cfl": args.cfl,
        "mach": mach,
        "ranks": ranks,
        "rho_l1": errors,
        "orders": orders,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"vortex matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"vortex matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)
