#!/usr/bin/env python3
"""Run and compare frozen SCMM6/Roe WENO5 and MDCD-HYBRID case04 jobs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY = CASE_DIR.parents[2]
METHODS = {
    "weno5": {
        "config": "scmm6_roe_weno5_480x120.wcns",
        "case": "case04-scmm6-roe-weno5-480x120",
        "output": "scmm6_roe_weno5_480x120",
    },
    "mdcd_hybrid": {
        "config": "scmm6_roe_mdcd_hybrid_480x120.wcns",
        "case": "case04-scmm6-roe-mdcd-hybrid-480x120",
        "output": "scmm6_roe_mdcd_hybrid_480x120",
    },
}


def default_executable(name: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return REPOSITORY / "build-rc-mpi" / f"{name}{suffix}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=default_executable("wcns_run"))
    parser.add_argument(
        "--generator", type=Path,
        default=default_executable("wcns_generate_release_cgns"),
    )
    parser.add_argument(
        "--validator", type=Path,
        default=default_executable("wcns_validate_release_case"),
    )
    parser.add_argument("--mpi-exec", type=Path, default=Path("mpiexec"))
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument(
        "--method", choices=("both", "weno5", "mdcd_hybrid"), default="both",
    )
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--dry-run-only", action="store_true")
    parser.add_argument(
        "--postprocess-only", action="store_true",
        help="validate and compare already completed result directories",
    )
    parser.add_argument("--skip-analysis", action="store_true")
    return parser.parse_args()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {resolved}")
    return resolved


def clean_comparison_outputs() -> None:
    targets = [
        CASE_DIR / "results/scmm6_roe_weno5_480x120",
        CASE_DIR / "results/scmm6_roe_mdcd_hybrid_480x120",
        CASE_DIR / "logs/comparison",
        CASE_DIR / "validation/comparison",
        CASE_DIR / "comparison",
    ]
    for target in targets:
        resolved = target.resolve()
        if CASE_DIR not in resolved.parents:
            raise RuntimeError(f"refusing unsafe cleanup target: {resolved}")
        if resolved.exists():
            shutil.rmtree(resolved)


def execute_streaming(
    command: list[str], log_path: Path, accepted: tuple[int, ...] = (0,),
) -> tuple[int, float]:
    print("running:", subprocess.list2cmdline(command), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        process = subprocess.Popen(
            command,
            cwd=CASE_DIR,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        returncode = process.wait()
    elapsed = time.perf_counter() - started
    if returncode not in accepted:
        raise RuntimeError(
            f"command failed with exit code {returncode}; see {log_path}")
    return returncode, elapsed


def require_single(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def manifest_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    if values.get("stop_reason") != "physical_time_reached":
        raise RuntimeError(f"unexpected stop reason in {path}: {values.get('stop_reason')}")
    if abs(float(values["time"]) - 0.2) > 1.0e-14:
        raise RuntimeError(f"final time is not 0.2 in {path}")
    return values


def main() -> int:
    args = parse_args()
    if args.ranks <= 0:
        raise ValueError("--ranks must be positive")
    if args.postprocess_only and (args.clean or args.dry_run_only):
        raise ValueError("--postprocess-only cannot be combined with --clean/--dry-run-only")
    if args.clean:
        clean_comparison_outputs()

    validator = require_file(args.validator, "field validator")
    selected = list(METHODS) if args.method == "both" else [args.method]
    elapsed: dict[str, float] = {}
    if not args.postprocess_only:
        run = require_file(args.run, "solver")
        generator = require_file(args.generator, "grid generator")
        mpi_exec = shutil.which(str(args.mpi_exec))
        if mpi_exec is None:
            raise FileNotFoundError(f"MPI launcher was not found: {args.mpi_exec}")
        prefix = [mpi_exec, "-n", str(args.ranks)]

        grid = CASE_DIR / "grids/double_mach_480x120.cgns"
        if not grid.exists():
            execute_streaming(
                [
                    str(generator), "rectangle", "grids/double_mach_480x120.cgns",
                    "480", "120", "4", "4.0", "1.0", "false",
                ],
                CASE_DIR / "logs/comparison/generate-grid.log",
            )
        else:
            print(f"reusing existing mesh: {grid}")

        for method in selected:
            config = METHODS[method]["config"]
            execute_streaming(
                prefix + [str(run), "--config", config, "--dry-run"],
                CASE_DIR / f"logs/comparison/dry-run-{method}-r{args.ranks}.log",
            )
        if args.dry_run_only:
            return 0

        for method in selected:
            specification = METHODS[method]
            output = CASE_DIR / "results" / specification["output"]
            if output.exists():
                raise RuntimeError(f"output exists; archive it or use --clean: {output}")
            _, elapsed[method] = execute_streaming(
                prefix + [
                    str(run), "--config", specification["config"],
                ],
                CASE_DIR / f"logs/comparison/run-{method}-r{args.ranks}.log",
            )

    validation = CASE_DIR / "validation/comparison"
    validation.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, dict[str, Path | dict[str, str]]] = {}
    for method in selected:
        specification = METHODS[method]
        directory = CASE_DIR / "results" / specification["output"]
        final_cgns = require_single(directory, "*.field.*.cgns")
        final_tecplot = require_single(directory, "*.field.*.dat")
        manifest = require_single(directory, "*.manifest.r*.txt")
        history = require_single(directory, "*.history.r*.txt")
        statistics = require_single(directory, "*.statistics.r*.txt")
        values = manifest_values(manifest)
        execute_streaming(
            [str(validator), "finite", str(final_cgns)],
            validation / f"finite-{method}.txt",
        )
        execute_streaming(
            [
                str(validator), "tecplot-consistency", str(final_cgns),
                str(final_tecplot), "1e-12",
            ],
            validation / f"tecplot-consistency-{method}.txt",
        )
        outputs[method] = {
            "cgns": final_cgns,
            "tecplot": final_tecplot,
            "manifest": manifest,
            "history": history,
            "statistics": statistics,
            "manifest_values": values,
        }

    if args.method == "both" and not args.skip_analysis:
        execute_streaming(
            [
                str(validator), "field-error",
                str(outputs["weno5"]["cgns"]),
                str(outputs["mdcd_hybrid"]["cgns"]),
            ],
            validation / "field-error-mdcd-minus-weno5.txt",
        )
        execute_streaming(
            [
                sys.executable, str(CASE_DIR / "analyze_case04.py"),
                str(outputs["weno5"]["tecplot"]),
                str(outputs["mdcd_hybrid"]["tecplot"]),
                "--weno-history", str(outputs["weno5"]["history"]),
                "--mdcd-history", str(outputs["mdcd_hybrid"]["history"]),
                "--weno-statistics", str(outputs["weno5"]["statistics"]),
                "--mdcd-statistics", str(outputs["mdcd_hybrid"]["statistics"]),
                "--output-directory", str(CASE_DIR / "comparison"),
            ],
            validation / "analysis.log",
        )

    manifest_ranks = {
        int(output["manifest_values"]["mpi_ranks"])
        for output in outputs.values()
    }
    if len(manifest_ranks) != 1:
        raise RuntimeError("comparison results use different MPI rank counts")
    summary = {
        "status": "completed" if not args.dry_run_only else "dry_run_completed",
        "ranks": manifest_ranks.pop(),
        "methods": {},
    }
    for method in selected:
        values = outputs[method]["manifest_values"]
        assert isinstance(values, dict)
        summary["methods"][method] = {
            "configuration": METHODS[method]["config"],
            "elapsed_seconds": elapsed.get(method, float(values["wall_time"])),
            "final_step": int(values["step"]),
            "final_time": float(values["time"]),
            "time_step": float(values["time_step"]),
            "wall_time_reported": float(values["wall_time"]),
            "stop_reason": values["stop_reason"],
            "git_commit": values["git_commit"],
        }
    validation.joinpath("run-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
