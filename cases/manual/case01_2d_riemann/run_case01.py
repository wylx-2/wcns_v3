#!/usr/bin/env python3
"""Generate, run and validate manual acceptance case01."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY = CASE_DIR.parents[2]
sys.path.insert(0, str(REPOSITORY / "tools"))

from process_metrics import run_measured  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=REPOSITORY / "build-rc-mpi/wcns_run.exe")
    parser.add_argument(
        "--generator",
        type=Path,
        default=REPOSITORY / "build-rc-mpi/wcns_generate_release_cgns.exe",
    )
    parser.add_argument(
        "--validator",
        type=Path,
        default=REPOSITORY / "build-rc-mpi/wcns_validate_release_case.exe",
    )
    parser.add_argument(
        "--mpi-exec",
        type=Path,
        default=Path("C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe"),
    )
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--clean", action="store_true")
    return parser.parse_args()


def remove_generated() -> None:
    for name in ("grids", "results", "logs", "validation"):
        target = (CASE_DIR / name).resolve()
        if target.parent != CASE_DIR:
            raise RuntimeError(f"refusing unsafe cleanup target: {target}")
        if target.exists():
            shutil.rmtree(target)
    for name in ("case01-summary.json", "files.sha256"):
        target = CASE_DIR / name
        if target.exists():
            target.unlink()


def execute(command: list[str], log: Path) -> dict[str, object]:
    print(f"running: {' '.join(command)}", flush=True)
    code, output, wall, peak_rss = run_measured(command, log)
    if code != 0:
        raise RuntimeError(f"command failed ({code}); see {log}")
    return {
        "command": command,
        "return_code": code,
        "wall_seconds": wall,
        "peak_process_tree_rss_bytes": peak_rss,
        "log": str(log.relative_to(CASE_DIR)),
        "last_line": output.strip().splitlines()[-1] if output.strip() else "",
    }


def clustered_coordinate(logical: float, center: float, strength: float) -> float:
    denominator = math.sinh(strength)
    if logical <= center:
        return center * (1.0 - math.sinh(strength * (center - logical) / center) / denominator)
    return (
        center
        + (1.0 - center) * math.sinh(strength * (logical - center) / (1.0 - center)) / denominator
    )


def spacing_metadata(clustered: bool) -> dict[str, object]:
    coordinates = [
        clustered_coordinate(index / 256.0, 0.5, 2.5) if clustered else index / 256.0
        for index in range(257)
    ]
    spacing = [right - left for left, right in zip(coordinates, coordinates[1:])]
    discontinuity_cell = next(
        index
        for index, (left, right) in enumerate(zip(coordinates, coordinates[1:]))
        if left <= 0.8 < right
    )
    return {
        "cells": [256, 256],
        "vertices": [257, 257],
        "zones_i": 4,
        "length": [1.0, 1.0],
        "cluster_center": [0.5, 0.5] if clustered else None,
        "cluster_strength": 2.5 if clustered else 0.0,
        "minimum_spacing": min(spacing),
        "maximum_spacing": max(spacing),
        "center_vertex": coordinates[128],
        "cell_containing_discontinuity": discontinuity_cell,
        "discontinuity_cell_bounds": [
            coordinates[discontinuity_cell],
            coordinates[discontinuity_cell + 1],
        ],
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    if args.ranks < 1:
        raise RuntimeError("rank count must be positive")
    executables = [args.run, args.generator, args.validator]
    if args.ranks > 1:
        executables.append(args.mpi_exec)
    executables = [path.resolve() for path in executables]
    for path in executables:
        if not path.is_file():
            raise RuntimeError(f"required executable does not exist: {path}")
    run_executable, generator, validator = executables[:3]
    mpi_executable = executables[3] if args.ranks > 1 else None

    if args.clean:
        remove_generated()
    for name in ("grids", "results", "logs", "validation"):
        (CASE_DIR / name).mkdir(parents=True, exist_ok=True)
    if any((CASE_DIR / "results" / name).exists() for name in ("uniform", "clustered")):
        raise RuntimeError("result directories already exist; rerun with --clean to replace them")

    old_cwd = Path.cwd()
    os.chdir(CASE_DIR)
    try:
        records: list[dict[str, object]] = []
        records.append(
            execute(
                [
                    str(generator),
                    "rectangle",
                    "grids/uniform_256x256.cgns",
                    "256",
                    "256",
                    "4",
                    "1.0",
                    "1.0",
                    "false",
                ],
                CASE_DIR / "logs/generate-uniform.log",
            )
        )
        records.append(
            execute(
                [
                    str(generator),
                    "clustered-rectangle",
                    "grids/clustered_256x256_x0p5_y0p5.cgns",
                    "256",
                    "256",
                    "4",
                    "1.0",
                    "1.0",
                    "0.5",
                    "0.5",
                    "2.5",
                    "false",
                ],
                CASE_DIR / "logs/generate-clustered.log",
            )
        )

        final_fields: dict[str, Path] = {}
        for name, config in (
            ("uniform", "uniform_256x256.wcns"),
            ("clustered", "clustered_256x256.wcns"),
        ):
            command = [str(run_executable), "--config", config]
            if mpi_executable is not None:
                command = [str(mpi_executable), "-n", str(args.ranks)] + command
            records.append(execute(command, CASE_DIR / f"logs/run-{name}.log"))
            fields = sorted((CASE_DIR / "results" / name).glob("*.field.*.cgns"))
            if len(fields) != 2:
                raise RuntimeError(f"{name} expected initial and final fields, found {len(fields)}")
            final = fields[-1]
            final_fields[name] = final
            for label, field in (("initial", fields[0]), ("final", final)):
                records.append(
                    execute(
                        [str(validator), "finite", str(field)],
                        CASE_DIR / f"validation/{name}-{label}-finite.txt",
                    )
                )
            records.append(
                execute(
                    [str(validator), "diagonal-symmetry", str(final), "5e-4"],
                    CASE_DIR / f"validation/{name}-final-symmetry.txt",
                )
            )

        summary = {
            "case": "case01_2d_riemann",
            "status": "passed",
            "ranks": args.ranks,
            "initial_discontinuity": [0.8, 0.8],
            "end_time": 0.3,
            "algorithm": {
                "profile": "phenglei_wcns",
                "reconstruction": "weno_z",
                "reconstruction_variables": "characteristic",
                "riemann": "hllc",
                "time_integrator": "ssprk3",
                "cfl": 0.4,
            },
            "grids": {
                "uniform": spacing_metadata(False),
                "clustered": spacing_metadata(True),
            },
            "final_fields": {
                name: str(path.relative_to(CASE_DIR)) for name, path in final_fields.items()
            },
            "records": records,
        }
        summary_path = CASE_DIR / "case01-summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

        generated = [summary_path]
        for directory in ("grids", "results", "logs", "validation"):
            generated.extend(path for path in (CASE_DIR / directory).rglob("*") if path.is_file())
        checksum_lines = [
            f"{sha256(path)}  {path.relative_to(CASE_DIR).as_posix()}" for path in sorted(generated)
        ]
        (CASE_DIR / "files.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    finally:
        os.chdir(old_cwd)
    print(f"case01 passed: {CASE_DIR / 'case01-summary.json'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"case01 failed: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1)
