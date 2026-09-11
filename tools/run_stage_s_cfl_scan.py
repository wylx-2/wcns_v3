#!/usr/bin/env python3
"""Classify empirical CFL outcomes for several already-rendered cases."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path


def set_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {key}, found {len(matches)}")
    return text[: matches[0].start()] + f"{key} = {value}" + text[matches[0].end() :]


def safe_token(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def parse_spec(value: str) -> tuple[str, Path, list[float]]:
    parts = value.split(":", 2)
    if len(parts) != 3 or not parts[0]:
        raise argparse.ArgumentTypeError("case must be LABEL:CONFIG:CFL,CFL,...")
    cfls = [float(item) for item in parts[2].split(",")]
    if len(cfls) < 2 or any(not math.isfinite(item) or item <= 0.0 for item in cfls):
        raise argparse.ArgumentTypeError("each case requires at least two positive CFLs")
    return parts[0], Path(parts[1]), cfls


def diagnostic_maximum(output: Path, column: str) -> int:
    files = sorted(output.glob("*.history.*.txt"))
    maximum = 0
    for path in files:
        header: list[str] | None = None
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith("# step "):
                header = line[2:].split()
                continue
            if not line or line.startswith("#") or header is None:
                continue
            values = line.split()
            if len(values) == len(header) and column in header:
                maximum = max(maximum, int(float(values[header.index(column)])))
    return maximum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--case", action="append", type=parse_spec, required=True)
    args = parser.parse_args()
    root = args.work_dir.resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    records: list[dict[str, object]] = []
    for label, config_path, cfl_values in args.case:
        base = config_path.read_text(encoding="utf-8")
        successful = 0
        for cfl in cfl_values:
            token = safe_token(cfl)
            output = root / f"{label}-cfl-{token}-output"
            text = set_key(base, "case.name", f"stage-s-{label}-cfl-{token}")
            text = set_key(text, "run.cfl", format(cfl, ".17g"))
            text = set_key(text, "output.directory", output.as_posix())
            text = set_key(text, "output.allow_existing", "false")
            config = root / f"{label}-cfl-{token}.wcns"
            config.write_text(text, encoding="utf-8")
            completed = subprocess.run(
                [str(args.run), "--config", str(config)],
                capture_output=True,
                text=True,
            )
            log_text = completed.stdout + completed.stderr
            log = root / f"{label}-cfl-{token}.log"
            log.write_text(log_text, encoding="utf-8")
            matches = re.findall(r"reason=([a-z_]+) step=(\d+) time=(\S+)", log_text)
            match = matches[-1] if matches else None
            reason = match[0] if match else "startup_failure"
            if completed.returncode == 1 or match is None:
                raise RuntimeError(f"unexpected startup failure for {label} CFL={cfl}: {log}")
            if reason != "numerical_failure":
                successful += 1
            records.append(
                {
                    "case": label,
                    "config_source": str(config_path.resolve()),
                    "cfl": cfl,
                    "exit_code": completed.returncode,
                    "stop_reason": reason,
                    "step": int(match[1]),
                    "time": float(match[2]),
                    "maximum_troubled_cells": diagnostic_maximum(output, "troubled_cells"),
                    "maximum_local_recomputations": diagnostic_maximum(
                        output, "local_recomputations"
                    ),
                    "maximum_step_retries": diagnostic_maximum(output, "step_retries"),
                    "log": str(log),
                }
            )
        if successful == 0:
            raise RuntimeError(f"CFL scan has no successful point for {label}")
    summary = {
        "matrix_version": 1,
        "stage": "S",
        "status": "passed",
        "interpretation": "empirical outcomes only; not universal stability limits",
        "records": records,
    }
    path = root / "matrix-summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"stage S empirical CFL scan passed: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"stage S CFL scan failed: {error}", file=sys.stderr)
        raise SystemExit(1)
