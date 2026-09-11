#!/usr/bin/env python3
"""Generate and run the classical double-Mach-reflection case.

The script deliberately keeps validation separate from simulation: a dry run proves
that the mesh, configuration and special boundaries assemble; the full command then
advances the production solver to t=0.2 and writes the configured outputs.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY = CASE_DIR.parents[2]


def default_executable(name: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return REPOSITORY / "build-rc-mpi" / f"{name}{suffix}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        type=Path,
        default=default_executable("wcns_run"),
        help="path to the production solver",
    )
    parser.add_argument(
        "--generator",
        type=Path,
        default=default_executable("wcns_generate_release_cgns"),
        help="path to the release CGNS generator",
    )
    parser.add_argument(
        "--mpi-exec",
        type=Path,
        default=Path("mpiexec"),
        help="MPI launcher; use an absolute path if it is not on PATH",
    )
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument(
        "--generate-only",
        action="store_true",
        help="only create the 960x240 CGNS mesh",
    )
    parser.add_argument(
        "--dry-run-only",
        action="store_true",
        help="generate the mesh and execute configuration/assembly validation only",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="remove this case's generated grid, logs and results before running",
    )
    return parser.parse_args()


def remove_generated() -> None:
    for relative in ("grids", "logs", "results"):
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


def execute(command: list[str], log_name: str) -> None:
    log_path = CASE_DIR / "logs" / log_name
    print("running:", subprocess.list2cmdline(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=CASE_DIR,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    print(completed.stdout, end="")
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}; see {log_path}")


def main() -> int:
    args = parse_args()
    if args.ranks <= 0:
        raise ValueError("--ranks must be positive")
    if args.clean:
        remove_generated()

    generator = require_file(args.generator, "grid generator")
    (CASE_DIR / "grids").mkdir(exist_ok=True)
    (CASE_DIR / "logs").mkdir(exist_ok=True)

    grid = CASE_DIR / "grids" / "double_mach_960x240.cgns"
    if not grid.exists():
        execute(
            [
                str(generator),
                "rectangle",
                "grids/double_mach_960x240.cgns",
                "960",
                "240",
                "4",
                "4.0",
                "1.0",
                "false",
            ],
            "generate-grid.log",
        )
    else:
        print(f"reusing existing mesh: {grid}")

    if args.generate_only:
        return 0

    run = require_file(args.run, "solver")
    execute(
        [str(run), "--config", "double_mach_960x240.wcns", "--dry-run"],
        "dry-run.log",
    )
    if args.dry_run_only:
        return 0

    mpi_exec = shutil.which(str(args.mpi_exec))
    if mpi_exec is None:
        raise FileNotFoundError(f"MPI launcher was not found: {args.mpi_exec}")
    execute(
        [
            mpi_exec,
            "-n",
            str(args.ranks),
            str(run),
            "--config",
            "double_mach_960x240.wcns",
        ],
        f"run-r{args.ranks}.log",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
