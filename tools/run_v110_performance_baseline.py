#!/usr/bin/env python3
"""Measure the three fixed v1.1.0 stage-P workloads without field output."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import statistics
import subprocess
import sys
from pathlib import Path

from process_metrics import run_measured
from run_release_matrix import clean_work_directory, render


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    return parser.parse_args()


def set_key(text: str, key: str, value: str) -> str:
    prefix = f"{key} = "
    lines = text.splitlines()
    matches = [index for index, line in enumerate(lines) if line.startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one configuration key {key}, found {len(matches)}")
    lines[matches[0]] = prefix + value
    return "\n".join(lines) + "\n"


def disable_output(text: str) -> str:
    for key in (
        "output.field.enabled",
        "output.history.enabled",
        "output.statistics.enabled",
        "output.checkpoint.enabled",
    ):
        text = set_key(text, key, "false")
    return text


def run_checked(command: list[str], log: Path) -> None:
    code, _, _, _ = run_measured(command, log)
    if code != 0:
        raise RuntimeError(f"command failed ({code}): {' '.join(command)}; see {log}")


def measure_case(
    executable: Path,
    config_text: str,
    case_root: Path,
    warmups: int,
    repetitions: int,
) -> dict[str, object]:
    samples: list[float] = []
    rss_samples: list[int | None] = []
    for index in range(warmups + repetitions):
        kind = "warmup" if index < warmups else "sample"
        ordinal = index + 1 if kind == "warmup" else index - warmups + 1
        run_root = case_root / f"{kind}-{ordinal}"
        run_root.mkdir(parents=True)
        config = run_root / "case.wcns"
        config.write_text(
            set_key(config_text, "output.directory", (run_root / "output").as_posix()),
            encoding="utf-8",
        )
        code, output, elapsed, peak_rss = run_measured(
            [str(executable), "--config", str(config)], run_root / "run.log"
        )
        if code != 0 and "reason=maximum_steps" not in output:
            raise RuntimeError(f"performance run failed; see {run_root / 'run.log'}")
        if kind == "sample":
            samples.append(elapsed)
            rss_samples.append(peak_rss)
    mean = statistics.fmean(samples)
    deviation = statistics.pstdev(samples)
    return {
        "wall_seconds": samples,
        "peak_process_tree_rss_bytes": rss_samples,
        "median_wall_seconds": statistics.median(samples),
        "mean_wall_seconds": mean,
        "standard_deviation_wall_seconds": deviation,
        "coefficient_of_variation": deviation / mean if mean > 0.0 else math.inf,
    }


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=10)
        return (result.stdout or result.stderr).splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unavailable"


def main() -> int:
    args = parse_args()
    if args.repetitions < 2 or args.warmups < 1:
        raise RuntimeError("stage-P baseline requires at least one warmup and two samples")
    root = args.work_dir.resolve()
    clean_work_directory(root)

    generator = args.generator.resolve()
    executable = args.run.resolve()
    repository = Path(__file__).resolve().parents[1]

    vortex_mesh = root / "vortex-100x100.cgns"
    run_checked(
        [str(generator), "periodic-square", str(vortex_mesh), "100", "100", "10.0"],
        root / "generate-vortex.log",
    )
    vortex = render(
        (repository / "cases/config/isentropic_vortex.wcns.in").read_text(encoding="utf-8"),
        {
            "CASE_NAME": "p-baseline-vortex-100x100",
            "MESH_PATH": vortex_mesh.as_posix(),
            "ALGORITHM_PROFILE": "scmm6_wcns",
            "RECONSTRUCTION": "weno_z",
            "RIEMANN": "hllc",
            "MOLAR_MASS": "0.029",
            "REFERENCE_VELOCITY": "340.0",
            "REFERENCE_TEMPERATURE": "288.15",
            "MIN_CELLS": "8",
            "CFL": "0.1",
            "END_TIME": "100.0",
            "OUTPUT_DIRECTORY": "replaced-per-run",
        },
    )
    vortex = disable_output(set_key(vortex, "run.max_steps", "20"))

    cylinder_mesh = root / "case07-cylinder-48x32.cgns"
    run_checked(
        [
            str(generator),
            "cylinder-o",
            str(cylinder_mesh),
            "48",
            "32",
            "4",
            "1.0",
            "8.0",
            "2.5",
        ],
        root / "generate-cylinder.log",
    )
    cylinder = (
        repository / "cases/manual/case07_2d_cylinder/configs/cylinder_mach5_euler.wcns"
    ).read_text(encoding="utf-8")
    cylinder = set_key(cylinder, "case.name", "p-baseline-case07-cylinder")
    cylinder = set_key(cylinder, "mesh.path", cylinder_mesh.as_posix())
    cylinder = set_key(cylinder, "run.max_steps", "20")
    cylinder = set_key(cylinder, "run.t_end", "100.0")
    cylinder = disable_output(set_key(cylinder, "output.directory", "replaced-per-run"))

    viscous_mesh = root / "viscous-36x48x36.cgns"
    run_checked(
        [
            str(generator),
            "periodic-channel",
            str(viscous_mesh),
            "36",
            "48",
            "36",
            "2",
            "2",
            "2.0",
            "1.0",
            "2.0",
            "1.5",
        ],
        root / "generate-viscous.log",
    )
    viscous = render(
        (repository / "cases/config/viscous_channel.wcns.in").read_text(encoding="utf-8"),
        {
            "CASE_NAME": "p-baseline-viscous-36x48x36",
            "MESH_PATH": viscous_mesh.as_posix(),
            "ALGORITHM_PROFILE": "scmm6_wcns",
            "RIEMANN": "hllc",
            "REFERENCE_VISCOSITY": "0.1",
            "MIN_CELLS": "8",
            "INITIAL_TYPE": "couette",
            "LOWER_VELOCITY": "0.0",
            "UPPER_VELOCITY": "0.0",
            "VELOCITY_CURVATURE": "0.0",
            "LOWER_TEMPERATURE": "1.0",
            "UPPER_TEMPERATURE": "1.0",
            "TEMPERATURE_CURVATURE": "0.0",
            "CFL": "0.05",
            "MAX_STEPS": "3",
            "MIN_STEPS": "3",
            "CHECK_INTERVAL": "1",
            "CONSECUTIVE_CHECKS": "1",
            "L2_ABSOLUTE": "1.0e-30",
            "L2_RELATIVE": "1.0e-30",
            "LINF_ABSOLUTE": "1.0e-30",
            "LINF_RELATIVE": "1.0e-30",
            "OUTPUT_DIRECTORY": "replaced-per-run",
        },
    )
    viscous = set_key(viscous, "run.mode", "unsteady")
    viscous += "run.t_end = 100.0\n"
    viscous = disable_output(viscous)

    cases = []
    for case_id, grid, zones, config_text in (
        ("vortex-2d-100x100", [100, 100, 1], 1, vortex),
        ("case07-cylinder-4zone-48x32", [48, 32, 1], 4, cylinder),
        ("viscous-3d-36x48x36", [36, 48, 36], 4, viscous),
    ):
        result = measure_case(
            executable, config_text, root / case_id, args.warmups, args.repetitions
        )
        cases.append({"id": case_id, "grid": grid, "zones": zones, **result})

    summary = {
        "schema_version": 1,
        "stage": "P",
        "program_version": "1.0.0-baseline-for-1.1.0",
        "machine": {
            "os": platform.platform(),
            "cpu": platform.processor() or os.environ.get("PROCESSOR_IDENTIFIER", "unknown"),
            "logical_cpus": os.cpu_count(),
            "python": platform.python_version(),
            "cmake": command_version(["cmake", "--version"]),
            "build_type": "Release",
        },
        "protocol": {
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "maximum_cv": 0.05,
            "intermediate_output": False,
        },
        "cases": cases,
    }
    path = root / "performance-baseline.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    unstable = [case["id"] for case in cases if case["coefficient_of_variation"] > 0.05]
    print(f"stage-P performance baseline written: {path}")
    if unstable:
        print("warning: CV exceeds 5% for " + ", ".join(unstable), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"stage-P performance baseline failed: {error}", file=sys.stderr)
        raise SystemExit(1)
