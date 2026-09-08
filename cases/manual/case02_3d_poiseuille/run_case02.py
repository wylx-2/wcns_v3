#!/usr/bin/env python3
"""Generate, run and validate manual acceptance case02."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
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
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse a completed case only after checking its manifest and current config digest",
    )
    parser.add_argument("--velocity-l2-tolerance", type=float, default=5.0e-3)
    parser.add_argument("--crossflow-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--pressure-span-tolerance", type=float, default=1.0e-2)
    parser.add_argument("--homogeneity-tolerance", type=float, default=1.0e-7)
    return parser.parse_args()


def remove_generated() -> None:
    for name in ("grids", "results", "logs", "validation"):
        target = (CASE_DIR / name).resolve()
        if target.parent != CASE_DIR:
            raise RuntimeError(f"refusing unsafe cleanup target: {target}")
        if target.exists():
            shutil.rmtree(target)
    for name in ("case02-summary.json", "files.sha256"):
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


def wall_coordinate(logical: float, strength: float) -> float:
    if strength == 0.0:
        return logical
    return 0.5 * (
        1.0 + math.tanh(strength * (2.0 * logical - 1.0)) / math.tanh(strength)
    )


def grid_metadata(strength: float) -> dict[str, object]:
    y = [wall_coordinate(index / 48.0, strength) for index in range(49)]
    spacing = [right - left for left, right in zip(y, y[1:])]
    return {
        "cells": [36, 48, 36],
        "vertices": [37, 49, 37],
        "zones": [2, 1, 2],
        "length": [2.0 * math.pi, 1.0, math.pi],
        "periodic_directions": ["x", "z"],
        "wall_cluster_strength": strength,
        "minimum_dy": min(spacing),
        "maximum_dy": max(spacing),
        "maximum_to_minimum_dy": max(spacing) / min(spacing),
        "center_vertex": y[24],
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def last_nonempty_line(path: Path) -> str:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    return lines[-1] if lines else ""


def poiseuille_validation_command(
    validator: Path,
    field: Path,
    velocity_l2: float,
    crossflow: float,
    pressure_span: float,
    homogeneity: float,
) -> list[str]:
    return [
        str(validator),
        "poiseuille-profile",
        str(field),
        "0.0",
        "1.0",
        "1.0",
        str(velocity_l2),
        str(crossflow),
        str(pressure_span),
        str(homogeneity),
    ]


def main() -> int:
    args = parse_args()
    if args.clean and args.resume:
        raise RuntimeError("--clean and --resume are mutually exclusive")
    if args.ranks < 1:
        raise RuntimeError("rank count must be positive")
    for value, label in (
        (args.velocity_l2_tolerance, "velocity L2 tolerance"),
        (args.crossflow_tolerance, "crossflow tolerance"),
        (args.pressure_span_tolerance, "pressure span tolerance"),
        (args.homogeneity_tolerance, "homogeneity tolerance"),
    ):
        if not math.isfinite(value) or value < 0.0:
            raise RuntimeError(f"{label} must be finite and nonnegative")
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
    if not args.resume and any(
        (CASE_DIR / "results" / name).exists() for name in ("uniform", "wall_clustered")
    ):
        raise RuntimeError("result directories already exist; rerun with --clean to replace them")

    old_cwd = Path.cwd()
    os.chdir(CASE_DIR)
    try:
        records: list[dict[str, object]] = []
        common_grid = [
            "36", "48", "36", "2", "2",
            str(2.0 * math.pi), "1.0", str(math.pi),
        ]
        grid_commands = (
            ("uniform", "grids/uniform_36x48x36.cgns", "0.0"),
            ("wall-clustered", "grids/wall_clustered_36x48x36.cgns", "1.0"),
        )
        for label, relative_path, strength in grid_commands:
            grid_path = CASE_DIR / relative_path
            command = [str(generator), "periodic-channel", relative_path] + common_grid + [strength]
            if args.resume:
                if not grid_path.is_file():
                    raise RuntimeError(f"resume grid is missing: {grid_path}")
                records.append({
                    "command": command,
                    "return_code": 0,
                    "wall_seconds": None,
                    "peak_process_tree_rss_bytes": None,
                    "log": f"logs/generate-{label}.log",
                    "last_line": "existing grid reused",
                    "reused_existing": True,
                    "sha256": sha256(grid_path),
                })
            else:
                records.append(execute(
                    command,
                    CASE_DIR / f"logs/generate-{label}.log",
                ))

        final_fields: dict[str, Path] = {}
        for name, config in (
            ("uniform", "uniform_36x48x36.wcns"),
            ("wall_clustered", "wall_clustered_36x48x36.wcns"),
        ):
            command = [str(run_executable), "--config", config]
            if mpi_executable is not None:
                command = [str(mpi_executable), "-n", str(args.ranks)] + command
            result_directory = CASE_DIR / "results" / name
            if args.resume and result_directory.exists():
                manifests = sorted(result_directory.glob(f"*.manifest.r{args.ranks}.txt"))
                if len(manifests) != 1:
                    raise RuntimeError(
                        f"{name} resume expected one rank-matched manifest, found {len(manifests)}"
                    )
                manifest = read_manifest(manifests[0])
                if manifest.get("stop_reason") != "steady_converged":
                    raise RuntimeError(f"{name} resume manifest is not steady_converged")
                dry_command = command + ["--dry-run"]
                dry_log = CASE_DIR / f"logs/resume-check-{name}.log"
                dry_record = execute(dry_command, dry_log)
                digest_match = re.search(
                    r"digest=0x([0-9a-fA-F]+)", dry_log.read_text(encoding="utf-8")
                )
                if digest_match is None:
                    raise RuntimeError(f"{name} resume could not read current config digest")
                current_digest = int(digest_match.group(1), 16)
                if current_digest != int(manifest.get("config_digest", "-1")):
                    raise RuntimeError(f"{name} resume config digest does not match manifest")
                records.append(dry_record | {"purpose": "resume configuration check"})
                run_log = CASE_DIR / f"logs/run-{name}.log"
                run_record = {
                    "command": command,
                    "return_code": 0,
                    "wall_seconds": float(manifest["wall_time"]),
                    "peak_process_tree_rss_bytes": None,
                    "log": str(run_log.relative_to(CASE_DIR)),
                    "last_line": last_nonempty_line(run_log),
                    "reused_existing": True,
                    "manifest": str(manifests[0].relative_to(CASE_DIR)),
                }
            else:
                run_record = execute(command, CASE_DIR / f"logs/run-{name}.log")
            if "reason=steady_converged" not in str(run_record["last_line"]):
                raise RuntimeError(f"{name} did not stop by steady residual convergence")
            records.append(run_record)
            fields = sorted((CASE_DIR / "results" / name).glob("*.field.*.cgns"))
            if len(fields) != 2:
                raise RuntimeError(f"{name} expected initial and final fields, found {len(fields)}")
            final_fields[name] = fields[-1]
            for label, field in (("initial", fields[0]), ("final", fields[-1])):
                records.append(execute(
                    [str(validator), "finite", str(field)],
                    CASE_DIR / f"validation/{name}-{label}-finite.txt",
                ))
            records.append(execute(
                poiseuille_validation_command(
                    validator, fields[0], 1.0e-12, 1.0e-12, 1.0e-10, 1.0e-10
                ),
                CASE_DIR / f"validation/{name}-initial-profile.txt",
            ))
            records.append(execute(
                poiseuille_validation_command(
                    validator,
                    fields[-1],
                    args.velocity_l2_tolerance,
                    args.crossflow_tolerance,
                    args.pressure_span_tolerance,
                    args.homogeneity_tolerance,
                ),
                CASE_DIR / f"validation/{name}-final-profile.txt",
            ))

        summary = {
            "case": "case02_3d_poiseuille",
            "status": "passed",
            "ranks": args.ranks,
            "domain": [2.0 * math.pi, 1.0, math.pi],
            "end_condition": "steady residual convergence",
            "steady_convergence": {
                "check_interval_steps": 20,
                "consecutive_checks": 3,
                "l2_absolute": 5.0e-5,
                "l2_relative": 1.0e-6,
                "linf_absolute": 2.0e-4,
                "linf_relative": 1.0e-5,
            },
            "physics": {
                "equations": "compressible Navier-Stokes",
                "reynolds": 10.0,
                "mach": 0.2,
                "prandtl": 0.72,
                "centerline_velocity": 1.0,
                "pressure_gradient_force": [0.8, 0.0, 0.0],
                "wall_temperature": 1.0,
            },
            "algorithm": {
                "profile": "phenglei_wcns",
                "reconstruction": "mdcd_linear",
                "reconstruction_variables": "primitive",
                "riemann": "hllc",
                "time_integrator": "ssprk3",
                "cfl": 0.5,
            },
            "grids": {
                "uniform": grid_metadata(0.0),
                "wall_clustered": grid_metadata(1.0),
            },
            "validation_tolerances": {
                "velocity_l2": args.velocity_l2_tolerance,
                "crossflow_linf": args.crossflow_tolerance,
                "pressure_span": args.pressure_span_tolerance,
                "homogeneity_linf": args.homogeneity_tolerance,
            },
            "final_fields": {
                name: str(path.relative_to(CASE_DIR)) for name, path in final_fields.items()
            },
            "records": records,
        }
        summary_path = CASE_DIR / "case02-summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

        generated = [summary_path]
        for directory in ("grids", "results", "logs", "validation"):
            generated.extend(path for path in (CASE_DIR / directory).rglob("*") if path.is_file())
        checksum_lines = [
            f"{sha256(path)}  {path.relative_to(CASE_DIR).as_posix()}"
            for path in sorted(generated)
        ]
        (CASE_DIR / "files.sha256").write_text(
            "\n".join(checksum_lines) + "\n", encoding="utf-8"
        )
    finally:
        os.chdir(old_cwd)
    print(f"case02 passed: {CASE_DIR / 'case02-summary.json'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"case02 failed: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1)
