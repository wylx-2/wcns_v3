#!/usr/bin/env python3
"""Run the isolated v1.1.0 release gate for the Case07 Mach 5 case."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import time
from pathlib import Path


MARKER = ".wcns-stage-u-case07"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec", required=True, type=Path)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--expected-version", default="1.1.0")
    parser.add_argument("--expected-commit")
    parser.add_argument("--timeout", type=float, default=900.0)
    return parser.parse_args()


def clean_work_directory(root: Path) -> None:
    root = root.resolve()
    marker = root / MARKER
    if root.exists():
        if not marker.is_file() or marker.read_text(encoding="utf-8") != MARKER + "\n":
            raise RuntimeError(
                f"refusing to clean unmarked Case07 work directory: {root}")
        for child in root.iterdir():
            if child == marker:
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    else:
        root.mkdir(parents=True)
    marker.write_text(MARKER + "\n", encoding="utf-8", newline="\n")


def set_key(text: str, key: str, value: str) -> str:
    expression = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    matches = list(expression.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {key}, found {len(matches)}")
    match = matches[0]
    return text[:match.start()] + f"{key} = {value}" + text[match.end():]


def key_value_file(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in result:
            raise RuntimeError(f"duplicate key in {path}: {key}")
        result[key] = value
    return result


def series_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    header: list[str] | None = None
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# step "):
            header = line[2:].split()
        elif line and not line.startswith("#"):
            rows.append(line.split())
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise RuntimeError(f"invalid fixed-schema series: {path}")
    return header, rows


def maximum_column(header: list[str], rows: list[list[str]], name: str) -> int:
    if name not in header:
        raise RuntimeError(f"history lacks {name}")
    index = header.index(name)
    return max(int(float(row[index])) for row in rows)


def finite_numeric_series(path: Path) -> int:
    count = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "TITLE=", "VARIABLES=", "ZONE ")):
            continue
        values = line.split()
        try:
            numbers = [float(value) for value in values]
        except ValueError as error:
            raise RuntimeError(f"non-numeric Case07 load row: {path}") from error
        if not numbers or not all(math.isfinite(value) for value in numbers):
            raise RuntimeError(f"non-finite Case07 load row: {path}")
        count += 1
    if count < 2:
        raise RuntimeError(f"Case07 load series is empty: {path}")
    return count


def one_file(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {pattern} in {root}, found {len(matches)}")
    return matches[0]


def execute(command: list[str], log: Path, timeout: float) -> dict[str, object]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )
    elapsed = time.perf_counter() - started
    log.write_text(completed.stdout, encoding="utf-8", newline="\n")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Case07 command failed ({completed.returncode}); see {log}")
    return {
        "command": command,
        "return_code": completed.returncode,
        "wall_seconds": elapsed,
        "log": str(log),
        "last_line": completed.stdout.strip().splitlines()[-1],
    }


def main() -> int:
    args = parse_args()
    if args.ranks != 4:
        raise RuntimeError("stage U freezes Case07 at 4 ranks")
    if args.timeout <= 0.0:
        raise RuntimeError("Case07 timeout must be positive")

    run = args.run.resolve()
    validator = args.validator.resolve()
    mpi = args.mpi_exec.resolve()
    source_config = args.config.resolve()
    for path in (run, validator, mpi, source_config):
        if not path.is_file():
            raise RuntimeError(f"required Case07 input does not exist: {path}")

    root = args.work_dir.resolve()
    clean_work_directory(root)
    output = root / "output"
    config_path = root / "case07-mach5-robust-u.wcns"
    text = source_config.read_text(encoding="utf-8")
    mesh_value = re.search(r"^mesh\.path\s*=\s*(.+)$", text, re.MULTILINE)
    if mesh_value is None:
        raise RuntimeError("Case07 source config lacks mesh.path")
    mesh = (source_config.parent / mesh_value.group(1).strip()).resolve()
    if not mesh.is_file():
        raise RuntimeError(f"Case07 mesh does not exist: {mesh}")

    text = set_key(text, "case.name", "stage-u-case07-mach5")
    text = set_key(text, "mesh.path", mesh.as_posix())
    text = set_key(text, "output.directory", output.as_posix())
    text = set_key(text, "output.allow_existing", "false")
    text = set_key(text, "output.field.format", "both")
    config_path.write_text(text, encoding="utf-8", newline="\n")

    command = [str(mpi), "-n", "4", str(run), "--config", str(config_path)]
    records = [execute(command, root / "run.log", args.timeout)]
    final_line = records[0]["last_line"]
    match = re.fullmatch(
        r"WCNS run stopped: reason=physical_time_reached step=(\d+) time=(\S+)",
        str(final_line),
    )
    if match is None or int(match.group(1)) <= 0 or abs(float(match.group(2)) - 8.0) > 1e-12:
        raise RuntimeError(f"Case07 did not stop exactly at t=8: {final_line}")
    step = int(match.group(1))

    manifest_path = one_file(output, "*.manifest.r4.txt")
    manifest = key_value_file(manifest_path)
    if manifest.get("program_version") != args.expected_version:
        raise RuntimeError("Case07 manifest program version mismatch")
    if args.expected_commit is not None and manifest.get("git_commit") != args.expected_commit:
        raise RuntimeError("Case07 manifest source commit mismatch")
    if manifest.get("stop_reason") != "physical_time_reached":
        raise RuntimeError("Case07 manifest stop reason mismatch")

    fields = sorted(output.glob("*.field.*.cgns"))
    final_fields = [path for path in fields if f"step{step:08d}" in path.name]
    if len(final_fields) != 1:
        raise RuntimeError("Case07 final CGNS field is missing or duplicated")
    validate_log = root / "validate-final.log"
    records.append(execute(
        [str(validator), "finite", str(final_fields[0])],
        validate_log,
        min(args.timeout, 60.0),
    ))

    history_path = one_file(output, "*.history.r4.txt")
    header, rows = series_rows(history_path)
    last = rows[-1]
    if int(float(last[header.index("step")])) != step:
        raise RuntimeError("Case07 history does not end at the final step")
    if last[header.index("stop_reason")] != "physical_time_reached":
        raise RuntimeError("Case07 history final stop reason mismatch")
    minimum_names = ("minimum_rho", "minimum_p", "minimum_T", "minimum_e")
    minimum_state = {
        name: min(
            float(row[header.index(name)])
            for row in rows
            if row[header.index(name)].lower() != "nan"
        )
        for name in minimum_names
    }
    if not all(math.isfinite(value) and value > 0.0
               for value in minimum_state.values()):
        raise RuntimeError("Case07 history contains a non-positive minimum state")

    loads_path = one_file(output, "*.loads.r4.txt")
    load_rows = finite_numeric_series(loads_path)
    temporary = sorted(str(path) for path in output.rglob("*.tmp"))
    if temporary:
        raise RuntimeError(f"Case07 left temporary files: {temporary}")
    latest = one_file(output, "*.checkpoint.latest.cgns")
    if latest.stat().st_size <= 0:
        raise RuntimeError("Case07 latest checkpoint is empty")

    summary = {
        "schema_version": 1,
        "stage": "U",
        "case": "case07_2d_cylinder_mach5_robust",
        "status": "passed",
        "ranks": 4,
        "step": step,
        "time": 8.0,
        "wall_seconds": records[0]["wall_seconds"],
        "program_version": manifest["program_version"],
        "source_commit": manifest.get("git_commit", "unknown"),
        "maximum_level1_faces": maximum_column(
            header, rows, "robustness_level1_faces"),
        "maximum_level2_faces": maximum_column(
            header, rows, "robustness_level2_faces"),
        "maximum_level3_faces": maximum_column(
            header, rows, "robustness_level3_faces"),
        "maximum_troubled_cells": maximum_column(
            header, rows, "troubled_cells"),
        "maximum_local_recomputations": maximum_column(
            header, rows, "local_recomputations"),
        "maximum_step_retries": maximum_column(
            header, rows, "step_retries"),
        "minimum_candidate_state": minimum_state,
        "load_rows": load_rows,
        "final_field": str(final_fields[0]),
        "history": str(history_path),
        "loads": str(loads_path),
        "manifest": str(manifest_path),
        "latest_checkpoint": str(latest),
        "records": records,
    }
    summary_path = root / "matrix-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"stage U Case07 matrix passed: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"stage U Case07 matrix failed: {error}")
        raise SystemExit(1)
