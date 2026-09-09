#!/usr/bin/env python3
"""Generate, dry-run, briefly advance and validate case05."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY = CASE_DIR.parents[2]
CONFIG = "channel_retau180_feasibility.wcns"


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
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--dry-run-only", action="store_true")
    parser.add_argument(
        "--clean", action="store_true",
        help="delete only case05 grids/logs/results/validation before running",
    )
    return parser.parse_args()


def remove_generated() -> None:
    for relative in ("grids", "logs", "results", "validation"):
        target = (CASE_DIR / relative).resolve()
        if target.parent != CASE_DIR:
            raise RuntimeError(f"refusing unsafe cleanup target: {target}")
        if target.exists():
            shutil.rmtree(target)


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {resolved}")
    return resolved


def execute(
    command: list[str], log_path: Path, acceptable_returncodes: tuple[int, ...] = (0,),
) -> None:
    print("running:", subprocess.list2cmdline(command), flush=True)
    completed = subprocess.run(
        command, cwd=CASE_DIR, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(completed.stdout, encoding="utf-8")
    print(completed.stdout, end="")
    if completed.returncode not in acceptable_returncodes:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; see {log_path}")


def launch_prefix(mpi_exec: str | None, ranks: int) -> list[str]:
    return [] if ranks == 1 else [str(mpi_exec), "-n", str(ranks)]


def main() -> int:
    args = parse_args()
    if args.ranks <= 0:
        raise ValueError("--ranks must be positive")
    if args.clean:
        remove_generated()
    for directory in ("grids", "logs", "validation"):
        (CASE_DIR / directory).mkdir(exist_ok=True)

    generator = require_file(args.generator, "grid generator")
    grid = CASE_DIR / "grids/channel_36x48x36_beta1p75.cgns"
    if not grid.exists():
        execute([
            str(generator), "periodic-channel",
            "grids/channel_36x48x36_beta1p75.cgns",
            "36", "48", "36", "2", "2",
            "6.283185307179586", "2.0", "3.141592653589793",
            "1.75", "-1.0",
        ], CASE_DIR / "logs/generate-grid.log")
    else:
        print(f"reusing existing grid: {grid}")
    if args.generate_only:
        return 0

    run = require_file(args.run, "solver")
    validator = require_file(args.validator, "field validator")
    mpi_exec: str | None = None
    if args.ranks > 1:
        mpi_exec = shutil.which(str(args.mpi_exec))
        if mpi_exec is None:
            raise FileNotFoundError(f"MPI launcher was not found: {args.mpi_exec}")
    prefix = launch_prefix(mpi_exec, args.ranks)
    execute(
        prefix + [str(run), "--config", CONFIG, "--dry-run"],
        CASE_DIR / f"logs/dry-run-r{args.ranks}.log",
    )
    if args.dry_run_only:
        return 0

    output = CASE_DIR / f"results/feasibility-r{args.ranks}"
    configured_output = CASE_DIR / "results/feasibility-r4"
    if args.ranks != 4:
        raise RuntimeError(
            "the frozen feasibility config writes results/feasibility-r4; "
            "use --ranks 4 for recorded acceptance, or copy the config and change "
            "output.directory before testing another rank count")
    if output.exists() or configured_output.exists():
        raise RuntimeError("feasibility output exists; archive it or rerun with --clean")
    execute(
        prefix + [str(run), "--config", CONFIG],
        CASE_DIR / f"logs/run-feasibility-r{args.ranks}.log",
        # wcns_run 用 2 明确表示达到 max_steps 后安全停止，本短测正是该情形。
        acceptable_returncodes=(2,),
    )

    fields = sorted(output.glob("*.field.*.cgns"))
    statistics = sorted(output.glob("*.statistics.r4.txt"))
    manifests = sorted(output.glob("*.manifest.r4.txt"))
    if len(fields) < 2 or len(statistics) != 1 or len(manifests) != 1:
        raise RuntimeError("feasibility run did not commit the expected output set")
    execute(
        [str(validator), "finite", str(fields[-1])],
        CASE_DIR / "validation/final-field-finite.txt",
    )
    execute(
        [sys.executable, str(CASE_DIR / "validate_case05.py"), str(statistics[0])],
        CASE_DIR / "validation/statistics-check.json",
    )
    print("case05 short feasibility validation completed; no turbulence acceptance claimed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
