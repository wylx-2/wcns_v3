#!/usr/bin/env python3
"""Exercise deterministic serial/MPI failure paths through wcns_run."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from run_release_matrix import clean_work_directory, render


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--timeout", type=float, default=60.0)
    return parser.parse_args()


def execute(
    command: list[str],
    log: Path,
    expected: tuple[int, ...],
    timeout: float,
    fragment: str = "",
) -> dict[str, object]:
    started = time.perf_counter()
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        log.write_text(output, encoding="utf-8")
        raise RuntimeError(f"command timed out (possible MPI deadlock): {log}") from error
    elapsed = time.perf_counter() - started
    log.write_text(result.stdout, encoding="utf-8")
    if result.returncode not in expected:
        raise RuntimeError(
            f"expected exit {expected}, got {result.returncode}: {log}"
        )
    if fragment and fragment.lower() not in result.stdout.lower():
        raise RuntimeError(f"expected diagnostic '{fragment}' is absent: {log}")
    return {
        "command": command,
        "return_code": result.returncode,
        "wall_seconds": elapsed,
        "diagnostic": fragment,
        "log": str(log),
        "last_line": result.stdout.strip().splitlines()[-1] if result.stdout.strip() else "",
    }


def digest_tree(directory: Path) -> str:
    result = hashlib.sha256()
    for path in sorted(item for item in directory.rglob("*") if item.is_file()):
        result.update(path.relative_to(directory).as_posix().encode("utf-8"))
        result.update(path.read_bytes())
    return result.hexdigest()


def assert_no_final_solution(directory: Path) -> None:
    if not directory.exists():
        return
    forbidden = [
        path for path in directory.iterdir()
        if ".field." in path.name or ".checkpoint." in path.name
        or ".manifest." in path.name
    ]
    if forbidden:
        raise RuntimeError(f"failed setup committed final-looking files: {forbidden}")


def main() -> int:
    args = parse_args()
    ranks = [int(value) for value in args.ranks.split(",")]
    if any(rank < 1 for rank in ranks) or (
        any(rank > 1 for rank in ranks) and args.mpi_exec is None
    ):
        raise RuntimeError("failure-matrix rank configuration is invalid")
    if args.timeout <= 0.0:
        raise RuntimeError("failure-matrix timeout must be positive")

    run_executable = args.run.resolve()
    generator = args.generator.resolve()
    mpi_executable = args.mpi_exec.resolve() if args.mpi_exec else None
    template = args.template.read_text(encoding="utf-8")
    clean_work_directory(args.work_dir)
    root = args.work_dir.resolve()
    records: list[dict[str, object]] = []

    valid_mesh = root / "valid.cgns"
    small_mesh = root / "small.cgns"
    one_sided_mesh = root / "one-sided.cgns"
    records.append(execute([
        str(generator), "periodic-square", str(valid_mesh), "16", "16", "1.0",
    ], root / "generate-valid.log", (0,), args.timeout))
    records.append(execute([
        str(generator), "rectangle", str(small_mesh), "12", "8", "1",
        "1.0", "1.0", "false",
    ], root / "generate-small.log", (0,), args.timeout))
    records.append(execute([
        str(generator), "invalid-one-sided", str(one_sided_mesh), "16", "8",
    ], root / "generate-one-sided.log", (0,), args.timeout))
    corrupt_mesh = root / "corrupt.cgns"
    corrupt_mesh.write_bytes(valid_mesh.read_bytes()[:64])

    def config_text(
        name: str,
        mesh: Path,
        output: Path,
        *,
        gamma: str = "1.4",
        min_cells: str = "4",
        source_enabled: str = "false",
        source_lines: str = "# no source parameters",
        cfl: str = "0.1",
        max_steps: str = "100",
        end_time: str = "0.01",
        field_quantities: str = "rho,u,v,w,p,T,rho_u,rho_v,rho_w,rho_E",
        history_enabled: str = "false",
        statistics_enabled: str = "false",
        statistic_quantities: str = "total_mass",
        checkpoint_enabled: str = "false",
        restart: Path | None = None,
    ) -> str:
        return render(template, {
            "CASE_NAME": name,
            "MESH_PATH": mesh.as_posix(),
            "GAMMA": gamma,
            "MIN_CELLS": min_cells,
            "SOURCE_ENABLED": source_enabled,
            "SOURCE_LINES": source_lines,
            "CFL": cfl,
            "MAX_STEPS": max_steps,
            "END_TIME": end_time,
            "OUTPUT_DIRECTORY": output.as_posix(),
            "FIELD_QUANTITIES": field_quantities,
            "HISTORY_ENABLED": history_enabled,
            "STATISTICS_ENABLED": statistics_enabled,
            "STATISTIC_QUANTITIES": statistic_quantities,
            "CHECKPOINT_ENABLED": checkpoint_enabled,
            "RESTART_LINE": (
                f"restart.path = {restart.as_posix()}" if restart
                else "# fresh start"
            ),
        })

    def write_config(name: str, text: str) -> Path:
        path = root / f"{name}.wcns"
        path.write_text(text, encoding="utf-8")
        return path

    def command(config: Path, rank: int) -> list[str]:
        result = [str(run_executable), "--config", str(config)]
        if rank > 1:
            assert mpi_executable is not None
            result = [str(mpi_executable), "-n", str(rank), *result]
        return result

    failures: list[tuple[str, Path, Path, str]] = []
    parser_output = root / "parser-output"
    base = config_text("parser", valid_mesh, parser_output)
    failures.extend([
        ("unknown-key", write_config("unknown-key", base + "unknown.option = 1\n"),
         parser_output, "unknown configuration key"),
        ("duplicate-key", write_config("duplicate-key", base + "run.cfl = 0.2\n"),
         parser_output, "duplicate configuration key"),
        ("direct-re", write_config("direct-re", base + "Re = 1000\n"),
         parser_output, "derived values and cannot be configured"),
        ("direct-ma", write_config("direct-ma", base + "Ma = 0.2\n"),
         parser_output, "derived values and cannot be configured"),
    ])
    unknown_field_output = root / "unknown-field-output"
    failures.append((
        "unknown-field",
        write_config("unknown-field", config_text(
            "unknown-field", valid_mesh, unknown_field_output,
            field_quantities="rho,not_registered",
        )),
        unknown_field_output,
        "unknown field quantity",
    ))
    duplicate_field_output = root / "duplicate-field-output"
    failures.append((
        "duplicate-field",
        write_config("duplicate-field", config_text(
            "duplicate-field", valid_mesh, duplicate_field_output,
            field_quantities="rho,rho",
        )),
        duplicate_field_output,
        "selected twice",
    ))
    unknown_statistic_output = root / "unknown-statistic-output"
    failures.append((
        "unknown-statistic",
        write_config("unknown-statistic", config_text(
            "unknown-statistic", valid_mesh, unknown_statistic_output,
            history_enabled="true", statistics_enabled="true",
            statistic_quantities="not_registered",
        )),
        unknown_statistic_output,
        "unknown statistic quantity",
    ))
    infeasible_output = root / "infeasible-output"
    failures.append((
        "infeasible-partition",
        write_config("infeasible-partition", config_text(
            "infeasible-partition", small_mesh, infeasible_output,
            min_cells="8",
        )),
        infeasible_output,
        "maximum feasible",
    ))
    one_sided_output = root / "one-sided-output"
    failures.append((
        "one-sided-connectivity",
        write_config("one-sided-connectivity", config_text(
            "one-sided-connectivity", one_sided_mesh, one_sided_output,
        )),
        one_sided_output,
        "reciprocal",
    ))
    corrupt_output = root / "corrupt-output"
    failures.append((
        "corrupt-cgns",
        write_config("corrupt-cgns", config_text(
            "corrupt-cgns", corrupt_mesh, corrupt_output,
        )),
        corrupt_output,
        "cg_open",
    ))
    output_file = root / "output-is-file"
    output_file.write_text("not a directory\n", encoding="utf-8")
    failures.append((
        "output-not-directory",
        write_config("output-not-directory", config_text(
            "output-not-directory", valid_mesh, output_file,
        )),
        output_file,
        "not a directory",
    ))

    for name, config, output, fragment in failures:
        for rank in ranks:
            # A 12x8 single zone is deliberately infeasible only when rank > 1.
            if name == "infeasible-partition" and rank == 1:
                continue
            log = root / f"{name}-r{rank}.log"
            records.append(execute(
                command(config, rank), log, (1,), args.timeout, fragment
            ))
            if output != output_file:
                assert_no_final_solution(output)
                if output.exists():
                    shutil.rmtree(output)

    # A real successful run creates immutable outputs. Reusing the directory
    # with allow_existing=false must fail without changing any byte.
    conflict_output = root / "conflict-output"
    conflict_config = write_config("output-conflict", config_text(
        "output-conflict", valid_mesh, conflict_output,
        history_enabled="true", statistics_enabled="true",
    ))
    records.append(execute(
        command(conflict_config, 1), root / "output-conflict-first.log",
        (0,), args.timeout, "physical_time_reached",
    ))
    before = digest_tree(conflict_output)
    for rank in ranks:
        records.append(execute(
            command(conflict_config, rank), root / f"output-conflict-r{rank}.log",
            (1,), args.timeout, "output directory already exists",
        ))
        if digest_tree(conflict_output) != before:
            raise RuntimeError("rejected output conflict changed committed files")

    # Create a valid checkpoint, then change a numerical signature input.
    checkpoint_output = root / "checkpoint-source-output"
    checkpoint_config = write_config("checkpoint-source", config_text(
        "checkpoint-source", valid_mesh, checkpoint_output,
        max_steps="1", end_time="10.0", checkpoint_enabled="true",
    ))
    records.append(execute(
        command(checkpoint_config, 1), root / "checkpoint-source.log",
        (2,), args.timeout, "maximum_steps",
    ))
    latest = next(checkpoint_output.glob("*.checkpoint.latest.cgns"), None)
    if latest is None:
        raise RuntimeError("checkpoint source did not write latest")
    for rank in ranks:
        mismatch_output = root / f"restart-mismatch-r{rank}-output"
        mismatch_config = write_config(
            f"restart-mismatch-r{rank}", config_text(
                f"restart-mismatch-r{rank}", valid_mesh, mismatch_output,
                gamma="1.3", restart=latest,
            ))
        records.append(execute(
            command(mismatch_config, rank), root / f"restart-mismatch-r{rank}.log",
            (1,), args.timeout, "signature differs",
        ))
        assert_no_final_solution(mismatch_output)

    # A finite configuration intentionally overflows its first source residual.
    # The run must return the dedicated numerical-failure exit and commit only
    # diagnostic history/manifest, never a final field or checkpoint.
    nonfinite_output = root / "nonfinite-output"
    nonfinite_config = write_config("nonfinite", config_text(
        "nonfinite", valid_mesh, nonfinite_output,
        source_enabled="true",
        source_lines=(
            "source.models = uniform_conservative\n"
            "source.uniform.rho = 1e308\n"
            "source.uniform.momentum_x = 1e308\n"
            "source.uniform.momentum_y = 1e308\n"
            "source.uniform.momentum_z = 0\n"
            "source.uniform.energy = 1e308"
        ),
        cfl="1e308", max_steps="2", end_time="10.0",
        history_enabled="true", statistics_enabled="false",
    ))
    for rank in ranks:
        output = nonfinite_output if rank == 1 else root / f"nonfinite-r{rank}-output"
        config = nonfinite_config
        if rank > 1:
            config = write_config(f"nonfinite-r{rank}", config_text(
                f"nonfinite-r{rank}", valid_mesh, output,
                source_enabled="true",
                source_lines=(
                    "source.models = uniform_conservative\n"
                    "source.uniform.rho = 1e308\n"
                    "source.uniform.momentum_x = 1e308\n"
                    "source.uniform.momentum_y = 1e308\n"
                    "source.uniform.momentum_z = 0\n"
                    "source.uniform.energy = 1e308"
                ),
                cfl="1e308", max_steps="2", end_time="10.0",
                history_enabled="true", statistics_enabled="false",
            ))
        records.append(execute(
            command(config, rank), root / f"nonfinite-r{rank}.log",
            (3,), args.timeout, "numerical_failure",
        ))
        if list(output.glob("*.field.*")) or list(output.glob("*.checkpoint.*")):
            raise RuntimeError("numerical failure committed a final field/checkpoint")
        manifest = next(output.glob("*.manifest.*.txt"), None)
        if manifest is None or "stop_reason=numerical_failure" not in manifest.read_text(
            encoding="utf-8"
        ):
            raise RuntimeError("numerical failure manifest is absent or ambiguous")

    summary = {
        "matrix_version": 1,
        "status": "passed",
        "ranks": ranks,
        "timeout_seconds": args.timeout,
        "failure_cases": [name for name, *_ in failures] + [
            "output-conflict", "restart-signature", "numerical-failure",
            "unknown-dependency-unit", "cyclic-dependency-unit",
        ],
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"failure matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"failure matrix failed: {error}", file=sys.stderr)
        raise SystemExit(1)
