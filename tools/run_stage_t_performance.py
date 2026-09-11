#!/usr/bin/env python3
"""Run the reproducible v1.1.0 stage-T performance acceptance matrix."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
from pathlib import Path

from process_metrics import run_measured
from run_release_matrix import clean_work_directory, render
from run_v110_performance_baseline import disable_output, set_key


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--allocation-probe", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--detailed", action="store_true")
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--mpi-run", type=Path)
    parser.add_argument("--scaling-repetitions", type=int, default=5)
    parser.add_argument("--scaling-warmups", type=int, default=1)
    parser.add_argument("--skip-scaling", action="store_true")
    parser.add_argument(
        "--reuse-serial-manifest", type=Path,
        help="reuse serial/allocation results from a prior stage-T manifest",
    )
    return parser.parse_args()


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command, check=True, capture_output=True, text=True, timeout=10
        )
        return (result.stdout or result.stderr).splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unavailable"


def run_checked(command: list[str], log: Path) -> None:
    code, _, _, _ = run_measured(command, log)
    if code != 0:
        raise RuntimeError(f"command failed ({code}): {' '.join(command)}; see {log}")


def statistics_for(samples: list[float]) -> dict[str, float | list[float]]:
    mean = statistics.fmean(samples)
    deviation = statistics.pstdev(samples)
    return {
        "wall_seconds": samples,
        "median_wall_seconds": statistics.median(samples),
        "mean_wall_seconds": mean,
        "min_wall_seconds": min(samples),
        "max_wall_seconds": max(samples),
        "standard_deviation_wall_seconds": deviation,
        "coefficient_of_variation": deviation / mean if mean > 0.0 else math.inf,
    }


def measure_case(
    command_prefix: list[str],
    executable: Path,
    config_text: str,
    case_root: Path,
    warmups: int,
    repetitions: int,
    cells: int,
    steps: int,
    detailed: bool,
) -> dict[str, object]:
    samples: list[float] = []
    rss_samples: list[int | None] = []
    commands: list[list[str]] = []
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
        command = [*command_prefix, str(executable), "--config", str(config)]
        commands.append(command)
        code, output, elapsed, peak_rss = run_measured(
            command,
            run_root / "run.log",
            sample_seconds=0.002 if detailed else 0.01,
        )
        if code != 0 and "reason=maximum_steps" not in output:
            raise RuntimeError(f"performance run failed; see {run_root / 'run.log'}")
        if kind == "sample":
            samples.append(elapsed)
            rss_samples.append(peak_rss)
    result: dict[str, object] = statistics_for(samples)
    median = float(result["median_wall_seconds"])
    result.update({
        "peak_process_tree_rss_bytes": rss_samples,
        "cell_stages_per_second": 3.0 * cells * steps / median,
        "commands": commands,
        "detailed_timing": detailed,
    })
    return result


def serial_cases(repository: Path, root: Path, generator: Path) -> list[dict[str, object]]:
    vortex_mesh = root / "vortex-100x100.cgns"
    run_checked(
        [str(generator), "periodic-square", str(vortex_mesh), "100", "100", "10.0"],
        root / "generate-vortex.log",
    )
    vortex = render(
        (repository / "cases/config/isentropic_vortex.wcns.in").read_text(encoding="utf-8"),
        {
            "CASE_NAME": "t-candidate-vortex-100x100",
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
        [str(generator), "cylinder-o", str(cylinder_mesh), "48", "32", "4", "1.0", "8.0", "2.5"],
        root / "generate-cylinder.log",
    )
    cylinder = (
        repository / "cases/manual/case07_2d_cylinder/configs/cylinder_mach5_euler.wcns"
    ).read_text(encoding="utf-8")
    cylinder = set_key(cylinder, "case.name", "t-candidate-case07-cylinder")
    cylinder = set_key(cylinder, "mesh.path", cylinder_mesh.as_posix())
    cylinder = set_key(cylinder, "run.max_steps", "20")
    cylinder = set_key(cylinder, "run.t_end", "100.0")
    cylinder = disable_output(set_key(cylinder, "output.directory", "replaced-per-run"))

    viscous_mesh = root / "viscous-36x48x36.cgns"
    run_checked(
        [str(generator), "periodic-channel", str(viscous_mesh), "36", "48", "36", "2", "2", "2.0", "1.0", "2.0", "1.5"],
        root / "generate-viscous.log",
    )
    viscous = viscous_config(repository, viscous_mesh, "t-candidate-viscous", 3)
    return [
        {"id": "vortex-2d-100x100", "grid": [100, 100, 1], "zones": 1,
         "cells": 100 * 100, "steps": 20, "config": vortex},
        {"id": "case07-cylinder-4zone-48x32", "grid": [48, 32, 1], "zones": 4,
         "cells": 48 * 32, "steps": 20, "config": cylinder},
        {"id": "viscous-3d-36x48x36", "grid": [36, 48, 36], "zones": 4,
         "cells": 36 * 48 * 36, "steps": 3, "config": viscous},
    ]


def viscous_config(repository: Path, mesh: Path, name: str, steps: int) -> str:
    text = render(
        (repository / "cases/config/viscous_channel.wcns.in").read_text(encoding="utf-8"),
        {
            "CASE_NAME": name,
            "MESH_PATH": mesh.as_posix(),
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
            "MAX_STEPS": str(steps),
            "MIN_STEPS": str(steps),
            "CHECK_INTERVAL": "1",
            "CONSECUTIVE_CHECKS": "1",
            "L2_ABSOLUTE": "1.0e-30",
            "L2_RELATIVE": "1.0e-30",
            "LINF_ABSOLUTE": "1.0e-30",
            "LINF_RELATIVE": "1.0e-30",
            "OUTPUT_DIRECTORY": "replaced-per-run",
        },
    )
    text = set_key(text, "run.mode", "unsteady") + "run.t_end = 100.0\n"
    return disable_output(text)


def allocation_result(probe: Path, root: Path) -> dict[str, object]:
    code, output, elapsed, rss = run_measured([str(probe)], root / "allocation.log")
    if code != 0:
        raise RuntimeError(f"allocation probe failed; see {root / 'allocation.log'}")
    values = {name: int(value) for name, value in re.findall(r"(\w+)=(\d+)", output)}
    baseline = values["baseline_allocations"]
    candidate = values["second_allocations"]
    return {
        **values,
        "reduction": 1.0 - candidate / baseline,
        "wall_seconds": elapsed,
        "peak_process_tree_rss_bytes": rss,
    }


def scaling_matrix(
    repository: Path,
    root: Path,
    generator: Path,
    mpiexec: Path,
    executable: Path,
    warmups: int,
    repetitions: int,
    detailed: bool,
) -> dict[str, object]:
    root.mkdir(parents=True, exist_ok=True)
    ranks = [1, 2, 4, 8]
    strong_grid = (48, 72, 48)
    strong_mesh = root / "strong-48x72x48.cgns"
    run_checked(
        [str(generator), "periodic-channel", str(strong_mesh), *map(str, strong_grid),
         "2", "2", "2.0", "1.0", "2.0", "1.5"],
        root / "generate-strong.log",
    )
    strong_config = viscous_config(repository, strong_mesh, "t-strong", 2)
    strong: list[dict[str, object]] = []
    for rank_count in ranks:
        measured = measure_case(
            [str(mpiexec), "-n", str(rank_count)], executable, strong_config,
            root / "strong" / f"r{rank_count}", warmups, repetitions,
            math.prod(strong_grid), 2, detailed,
        )
        strong.append({
            "ranks": rank_count,
            "cells_per_rank": math.prod(strong_grid) // rank_count,
            **measured,
        })
    t1 = float(strong[0]["median_wall_seconds"])
    for item in strong:
        rank_count = int(item["ranks"])
        item["efficiency"] = t1 / (rank_count * float(item["median_wall_seconds"]))

    weak: list[dict[str, object]] = []
    # Both periodic directions are split into two source zones, so their
    # global counts must remain divisible by two at every rank count.
    local_grid = (24, 24, 28)
    local_cells = math.prod(local_grid)
    for rank_count in ranks:
        grid = (local_grid[0] * rank_count, local_grid[1], local_grid[2])
        mesh = root / f"weak-r{rank_count}.cgns"
        run_checked(
            [str(generator), "periodic-channel", str(mesh), *map(str, grid),
             "2", "2", "2.0", "1.0", "2.0", "1.5"],
            root / f"generate-weak-r{rank_count}.log",
        )
        config = viscous_config(repository, mesh, f"t-weak-r{rank_count}", 2)
        measured = measure_case(
            [str(mpiexec), "-n", str(rank_count)], executable, config,
            root / "weak" / f"r{rank_count}", warmups, repetitions,
            local_cells * rank_count, 2, detailed,
        )
        weak.append({
            "ranks": rank_count,
            "cells_per_rank": local_cells,
            "grid": list(grid),
            **measured,
        })
    weak_t1 = float(weak[0]["median_wall_seconds"])
    for item in weak:
        item["efficiency"] = weak_t1 / float(item["median_wall_seconds"])
    return {"strong": strong, "weak": weak}


def main() -> int:
    args = parse_args()
    if args.repetitions < 5 or args.warmups < 1:
        raise RuntimeError("stage T requires at least one warmup and five samples")
    if not args.skip_scaling and args.scaling_repetitions < 5:
        raise RuntimeError("stage T scaling requires five samples per rank")
    root = args.work_dir.resolve()
    clean_work_directory(root)
    repository = Path(__file__).resolve().parents[1]
    executable = args.run.resolve()
    generator = args.generator.resolve()
    baseline = json.loads(args.baseline.resolve().read_text(encoding="utf-8"))
    baseline_by_id = {item["id"]: item for item in baseline["cases"]}

    serial: list[dict[str, object]] = []
    allocation: dict[str, object]
    if args.reuse_serial_manifest is not None:
        prior = json.loads(
            args.reuse_serial_manifest.resolve().read_text(encoding="utf-8")
        )
        if prior.get("stage") != "T" or not prior.get("serial"):
            raise RuntimeError("reused serial manifest is not a stage-T result")
        serial = prior["serial"]
        allocation = prior["allocation"]
    else:
        for case in serial_cases(repository, root, generator):
            measured = measure_case(
                [], executable, str(case["config"]),
                root / "serial" / str(case["id"]), args.warmups,
                args.repetitions, int(case["cells"]), int(case["steps"]),
                args.detailed,
            )
            baseline_seconds = float(
                baseline_by_id[str(case["id"])]["median_wall_seconds"]
            )
            serial.append({
                key: value for key, value in case.items() if key != "config"
            } | measured | {
                "baseline_median_wall_seconds": baseline_seconds,
                "speedup": baseline_seconds / float(measured["median_wall_seconds"]),
            })
        allocation = allocation_result(args.allocation_probe.resolve(), root)

    scaling: dict[str, object] | None = None
    if not args.skip_scaling:
        if args.mpiexec is None or args.mpi_run is None:
            raise RuntimeError("scaling requires --mpiexec and --mpi-run")
        scaling = scaling_matrix(
            repository, root / "scaling", generator, args.mpiexec.resolve(),
            args.mpi_run.resolve(), args.scaling_warmups,
            args.scaling_repetitions, args.detailed,
        )

    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repository, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
    except subprocess.SubprocessError:
        commit = "unknown"
    summary = {
        "schema_version": 1,
        "stage": "T",
        "commit": commit,
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
            "detailed_timing": args.detailed,
            "intermediate_output": False,
        },
        "allocation": allocation,
        "serial": serial,
        "scaling": scaling,
    }
    output = root / "stage-t-performance.json"
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    failures: list[str] = []
    for item in serial:
        if float(item["coefficient_of_variation"]) > 0.05:
            failures.append(f"{item['id']}: CV exceeds 5%")
        minimum = 1.0 if item["id"] == "case07-cylinder-4zone-48x32" else 1.5
        if float(item["speedup"]) < minimum:
            failures.append(f"{item['id']}: speedup is below {minimum}")
    if float(summary["allocation"]["reduction"]) < 0.90:
        failures.append("allocation reduction is below 90%")
    if scaling is not None:
        strong = scaling["strong"]
        four = next(item for item in strong if item["ranks"] == 4)
        if int(four["cells_per_rank"]) < 15552:
            failures.append("strong scaling has too few cells per rank")
        if float(four["efficiency"]) < 0.70:
            failures.append("four-rank strong-scaling efficiency is below 70%")
        for family in (scaling["strong"], scaling["weak"]):
            for item in family:
                if float(item["coefficient_of_variation"]) > 0.05:
                    failures.append(
                        f"scaling r{item['ranks']}: CV exceeds 5%"
                    )
    print(f"stage-T performance manifest written: {output}")
    if failures:
        for failure in failures:
            print("FAIL: " + failure, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"stage-T performance run failed: {error}", file=sys.stderr)
        raise SystemExit(1)
