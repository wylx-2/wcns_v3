#!/usr/bin/env python3
"""Run deterministic WCNS release cases through the production executable."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path


MARKER = ".wcns-release-matrix"


def clean_work_directory(path: Path) -> None:
    resolved = path.resolve()
    if resolved.parent == resolved or len(resolved.parts) < 3:
        raise RuntimeError(f"refusing broad release work directory: {resolved}")
    if resolved.exists():
        if not (resolved / MARKER).is_file():
            raise RuntimeError(
                f"refusing to replace unmarked release work directory: {resolved}"
            )
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    (resolved / MARKER).write_text("wcns release matrix\n", encoding="utf-8")


def run(command: list[str], log_path: Path) -> dict[str, object]:
    start = time.perf_counter()
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.perf_counter() - start
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"see {log_path}"
        )
    return {
        "command": command,
        "return_code": completed.returncode,
        "wall_seconds": elapsed,
        "log": str(log_path),
        "last_line": completed.stdout.strip().splitlines()[-1]
        if completed.stdout.strip()
        else "",
    }


def render(template: str, replacements: dict[str, str]) -> str:
    result = template
    for key, value in replacements.items():
        result = result.replace(f"@{key}@", value)
    if "@" in result:
        unresolved = sorted({part.split("@", 1)[0] for part in result.split("@")[1::2]})
        raise RuntimeError(f"unresolved configuration placeholders: {unresolved}")
    return result


def one_field(directory: Path) -> Path:
    fields = sorted(directory.glob("*.field.*.cgns"))
    if len(fields) != 1:
        raise RuntimeError(
            f"expected exactly one final CGNS field in {directory}, found {len(fields)}"
        )
    return fields[0]


def one_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one {pattern} in {directory}, found {len(matches)}"
        )
    return matches[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--cells-i", type=int, default=16)
    parser.add_argument("--cells-j", type=int, default=8)
    parser.add_argument("--cells-k", type=int, default=1)
    parser.add_argument("--dimension", type=int, default=2)
    parser.add_argument("--zones-i", type=int, default=1)
    parser.add_argument("--warp", type=float, default=0.0)
    parser.add_argument("--periodic-x", action="store_true")
    parser.add_argument("--profile", default="phenglei_wcns")
    parser.add_argument("--reconstruction", default="weno_z")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--uniform-tolerance", type=float, default=1.0e-11)
    parser.add_argument("--case-prefix", default="release-freestream")
    parser.add_argument("--reference-directory", type=Path)
    parser.add_argument("--min-cells", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if not ranks or any(value < 1 for value in ranks) or len(set(ranks)) != len(ranks):
        raise RuntimeError("--ranks must contain distinct positive integers")
    if any(value > 1 for value in ranks) and args.mpi_exec is None:
        raise RuntimeError("MPI ranks require --mpi-exec")
    if args.steps < 1 or args.min_cells < 1 or args.uniform_tolerance < 0.0:
        raise RuntimeError("steps/min-cells must be positive and tolerance nonnegative")
    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    mesh = root / "release-grid.cgns"
    records: list[dict[str, object]] = []
    records.append(
        run(
            [
                str(args.generator),
                str(mesh),
                str(args.dimension),
                str(args.cells_i),
                str(args.cells_j),
                str(args.cells_k),
                str(args.zones_i),
                str(args.warp),
                "true" if args.periodic_x else "false",
            ],
            root / "generate.log",
        )
    )
    template = args.template.read_text(encoding="utf-8")
    fields: dict[int, Path] = {}
    for rank_count in ranks:
        case_name = f"{args.case_prefix}-r{rank_count}"
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
                    "MAX_STEPS": str(args.steps),
                    "MIN_STEPS": str(args.steps),
                    "MIN_CELLS": str(args.min_cells),
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
                    str(args.validator),
                    "uniform",
                    str(fields[rank_count]),
                    "1.0",
                    "0.2",
                    "-0.1",
                    "0.0",
                    "1.0",
                    str(args.uniform_tolerance),
                ],
                root / f"uniform-r{rank_count}.log",
            )
        )
        records.append(
            run(
                [str(args.validator), "finite", str(fields[rank_count])],
                root / f"finite-r{rank_count}.log",
            )
        )
        statistics = one_file(output, f"*.statistics.r{rank_count}.txt")
        records.append(
            run(
                [
                    str(args.validator),
                    "series-constant",
                    str(statistics),
                    "1e-12",
                ],
                root / f"conservation-r{rank_count}.log",
            )
        )
    reference_rank = ranks[0]
    for rank_count in ranks[1:]:
        records.append(
            run(
                [
                    str(args.validator),
                    "compare",
                    str(fields[reference_rank]),
                    str(fields[rank_count]),
                    "2e-11",
                ],
                root / f"compare-r{reference_rank}-r{rank_count}.log",
            )
        )
    if args.reference_directory is not None:
        reference_field = one_field(args.reference_directory.resolve())
        for rank_count in ranks:
            records.append(
                run(
                    [
                        str(args.validator),
                        "compare-spatial",
                        str(reference_field),
                        str(fields[rank_count]),
                        "2e-11",
                        "1e-13",
                    ],
                    root / f"compare-spatial-reference-r{rank_count}.log",
                )
            )
    summary = {
        "matrix_version": 1,
        "status": "passed",
        "grid": {
            "dimension": args.dimension,
            "cells": [args.cells_i, args.cells_j, args.cells_k],
            "zones_i": args.zones_i,
            "warp": args.warp,
            "periodic_x": args.periodic_x,
            "profile": args.profile,
            "reconstruction": args.reconstruction,
            "riemann": args.riemann,
            "steps": args.steps,
        },
        "ranks": ranks,
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"release matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:  # concise, deterministic CI diagnostic
        print(f"release matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)
