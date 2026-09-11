#!/usr/bin/env python3
"""Generate, run, validate and summarize manual acceptance case03."""

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
        "--metric-comparator",
        type=Path,
        default=REPOSITORY / "build-rc-mpi/wcns_compare_metric_profiles.exe",
    )
    parser.add_argument(
        "--mpi-exec",
        type=Path,
        default=Path("C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe"),
    )
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--vortex-density-l1-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--conservation-tolerance", type=float, default=2.0e-10)
    return parser.parse_args()


def remove_generated() -> None:
    for name in ("grids", "results", "logs", "validation"):
        target = (CASE_DIR / name).resolve()
        if target.parent != CASE_DIR:
            raise RuntimeError(f"refusing unsafe cleanup target: {target}")
        if target.exists():
            shutil.rmtree(target)
    for name in ("case03-summary.json", "files.sha256"):
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_value(value: str) -> object:
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    try:
        numeric = float(value)
    except ValueError:
        return value
    return numeric if math.isfinite(numeric) else value


def parse_records(path: Path, selector: str) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(selector):
            continue
        record: dict[str, object] = {}
        for token in line.split():
            if "=" in token:
                key, value = token.split("=", 1)
                record[key] = parse_value(value)
        result.append(record)
    return result


def manifest(path: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = parse_value(value)
    return result


def main() -> int:
    args = parse_args()
    mach = 340.0 / math.sqrt(1.4 * (8.314 / 0.029) * 288.15)
    if args.ranks != 4:
        raise RuntimeError("case03 uses four native zones and requires --ranks 4")
    for value, label in (
        (args.vortex_density_l1_tolerance, "vortex density L1 tolerance"),
        (args.conservation_tolerance, "conservation tolerance"),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError(f"{label} must be positive and finite")

    executables = [
        args.run,
        args.generator,
        args.validator,
        args.metric_comparator,
        args.mpi_exec,
    ]
    executables = [path.resolve() for path in executables]
    for path in executables:
        if not path.is_file():
            raise RuntimeError(f"required executable does not exist: {path}")
    run_executable, generator, validator, metric_comparator, mpi_executable = executables

    if args.clean:
        remove_generated()
    for name in ("grids", "logs", "validation"):
        (CASE_DIR / name).mkdir(parents=True, exist_ok=True)
    if not args.dry_run and (CASE_DIR / "results").exists():
        raise RuntimeError("results already exist; rerun with --clean to replace them")

    old_cwd = Path.cwd()
    os.chdir(CASE_DIR)
    try:
        records: list[dict[str, object]] = []
        grid_command = [
            str(generator),
            "warped-periodic-square",
            "grids/warped_100x100.cgns",
            "100",
            "100",
            "10.0",
            "0.5",
            "0.5",
        ]
        records.append(execute(grid_command, CASE_DIR / "logs/generate-grid.log"))
        records.append(
            execute(
                [str(metric_comparator), "grids/warped_100x100.cgns"],
                CASE_DIR / "validation/metric-profile-comparison.txt",
            )
        )

        profiles: dict[str, dict[str, object]] = {}
        for name, config_name in (
            ("phenglei", "phenglei_100x100.wcns"),
            ("scmm6", "scmm6_100x100.wcns"),
        ):
            command = [
                str(mpi_executable),
                "-n",
                str(args.ranks),
                str(run_executable),
                "--config",
                config_name,
            ]
            if args.dry_run:
                command.append("--dry-run")
                records.append(execute(command, CASE_DIR / f"logs/dry-run-{name}.log"))
                continue

            records.append(execute(command, CASE_DIR / f"logs/run-{name}.log"))
            fields = sorted((CASE_DIR / "results" / name).glob("*.field.*.cgns"))
            manifests = sorted((CASE_DIR / "results" / name).glob("*.manifest.r4.txt"))
            statistics = sorted((CASE_DIR / "results" / name).glob("*.statistics.r4.txt"))
            if len(fields) != 2 or len(manifests) != 1 or len(statistics) != 1:
                raise RuntimeError(
                    f"{name} expected two fields, one manifest and one statistics file"
                )
            initial, final = fields
            run_manifest = manifest(manifests[0])
            if run_manifest.get("stop_reason") != "physical_time_reached":
                raise RuntimeError(f"{name} did not stop at the requested physical time")

            for label, field in (("initial", initial), ("final", final)):
                records.append(
                    execute(
                        [str(validator), "finite", str(field)],
                        CASE_DIR / f"validation/{name}-{label}-finite.txt",
                    )
                )
            analytic_errors: dict[str, object] = {}
            for label, field, time in (
                ("initial", initial, 0.0),
                ("final", final, 10.0),
            ):
                analytic_log = CASE_DIR / f"validation/{name}-{label}-analytic-error.txt"
                records.append(
                    execute(
                        [
                            str(validator),
                            "vortex",
                            str(field),
                            str(time),
                            "10.0",
                            "5.0",
                            "5.0",
                            "5.0",
                            "1.0",
                            "1.0",
                            "1.4",
                            str(mach),
                            str(args.vortex_density_l1_tolerance),
                        ],
                        analytic_log,
                    )
                )
                parsed = parse_records(analytic_log, "check=isentropic_vortex")
                analytic_errors[label] = parsed[0] if len(parsed) == 1 else parsed
            field_error_log = CASE_DIR / f"validation/{name}-initial-final-error.txt"
            records.append(
                execute(
                    [str(validator), "field-error", str(initial), str(final)],
                    field_error_log,
                )
            )
            conservation_log = CASE_DIR / f"validation/{name}-conservation.txt"
            records.append(
                execute(
                    [
                        str(validator),
                        "series-constant",
                        str(statistics[0]),
                        str(args.conservation_tolerance),
                    ],
                    conservation_log,
                )
            )
            profiles[name] = {
                "manifest": run_manifest,
                "initial_field": str(initial.relative_to(CASE_DIR)),
                "final_field": str(final.relative_to(CASE_DIR)),
                "analytic_errors": analytic_errors,
                "initial_final_errors": parse_records(field_error_log, "check=field_error field="),
                "conservation": parse_records(conservation_log, "check=series_constant"),
            }

        if args.dry_run:
            print("case03 dry-run passed", flush=True)
            return 0

        metric_log = CASE_DIR / "validation/metric-profile-comparison.txt"
        metric_profiles = parse_records(metric_log, "profile=")
        metric_differences = parse_records(metric_log, "comparison=")
        if len(metric_profiles) != 2:
            raise RuntimeError("metric comparison did not report both profiles")
        for metric in metric_profiles:
            if metric["fallback_cells"] != 0:
                raise RuntimeError("metric construction used a forbidden fallback")
            if abs(float(metric["volume_sum"]) - 100.0) > 1.0e-10:
                raise RuntimeError("metric volume sum is inconsistent with the domain area")
            if float(metric["max_reference_relative_difference"]) > 1.0e-3:
                raise RuntimeError("high-order Jacobian differs excessively from cell volume")
            if float(metric["gcl_closure_linf"]) > 1.0e-10:
                raise RuntimeError("metric geometric-conservation closure exceeds tolerance")
        for metric in metric_differences:
            if metric["comparison"] == "jacobian" and float(metric["relative_linf"]) > 1.0e-3:
                raise RuntimeError("profile Jacobians differ beyond the case03 tolerance")

        summary = {
            "case": "case03_2d_isentropic_vortex_warped_grid",
            "status": "passed",
            "ranks": args.ranks,
            "domain": [10.0, 10.0],
            "cells": [100, 100],
            "native_zones": [2, 2],
            "warp": {
                "x_amplitude": 0.5,
                "x_eta_wave_count": 1,
                "y_amplitude": 0.5,
                "y_xi_wave_count": 2,
                "continuous_mapping_jacobian_range": [
                    1.0 - 2.0 * math.pi * math.pi / 100.0,
                    1.0 + 2.0 * math.pi * math.pi / 100.0,
                ],
            },
            "physics": {
                "equations": "compressible Euler",
                "gamma": 1.4,
                "mach": mach,
                "vortex_strength": 5.0,
                "initial_center": [5.0, 5.0],
                "background_velocity": [1.0, 1.0],
                "period": 10.0,
            },
            "algorithm_common": {
                "reconstruction": "weno_z",
                "reconstruction_variables": "characteristic",
                "riemann": "hllc",
                "time_integrator": "ssprk3",
                "cfl": 0.5,
            },
            "metric_acceptance": {
                "maximum_volume_sum_error": 1.0e-10,
                "maximum_reference_relative_difference": 1.0e-3,
                "maximum_gcl_closure_linf": 1.0e-10,
                "maximum_cross_profile_jacobian_relative_linf": 1.0e-3,
            },
            "metric_profiles": metric_profiles,
            "metric_differences": metric_differences,
            "profiles": profiles,
            "records": records,
        }
        summary_path = CASE_DIR / "case03-summary.json"
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
    print(f"case03 passed: {CASE_DIR / 'case03-summary.json'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"case03 failed: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1)
