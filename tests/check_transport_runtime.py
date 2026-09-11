#!/usr/bin/env python3
"""Run and verify the production Sutherland transport path."""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mpi-exec", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()
    if args.ranks < 1 or (args.ranks > 1 and args.mpi_exec is None):
        raise RuntimeError("invalid transport runtime rank selection")
    if args.output.exists():
        shutil.rmtree(args.output)
    command = [str(args.run), "--config", str(args.config)]
    if args.ranks > 1:
        command = [str(args.mpi_exec), "-n", str(args.ranks), *command]
    completed = subprocess.run(command, capture_output=True, text=True)
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise RuntimeError(f"transport runtime failed ({completed.returncode}):\n{output}")
    required = (
        "law=sutherland",
        "mu_Tref_over_mu_ref=1.25",
        "Pr=0.70999999999999996",
        "S_over_Tref=0.38313378448724628",
        "T_range=[",
        "mu_range=[",
    )
    missing = [token for token in required if token not in output]
    if missing:
        raise RuntimeError(f"startup transport diagnostics missing {missing}:\n{output}")
    files = sorted(args.output.glob("*.boundary.*.txt"))
    if len(files) != 1:
        raise RuntimeError(f"expected one boundary output, found {len(files)}")
    rows: list[dict[str, str]] = []
    with files[0].open(encoding="utf-8") as stream:
        header: list[str] | None = None
        for line in stream:
            clean = line.strip()
            if clean.startswith("# patch_index "):
                header = clean[2:].split()
                continue
            if not clean or clean.startswith("#"):
                continue
            if header is None:
                raise RuntimeError("boundary output is missing its column header")
            values = clean.split()
            if len(values) != len(header):
                raise RuntimeError("boundary row width differs from its header")
            rows.append(dict(zip(header, values)))
    if not rows:
        raise RuntimeError("boundary output has no face records")
    s = 110.4 / 288.15
    expected = {
        1.0: 1.25,
        2.0: 1.25 * math.pow(2.0, 1.5) * (1.0 + s) / (2.0 + s),
    }
    seen: set[float] = set()
    for row in rows:
        temperature = float(row["T_w"])
        viscosity = float(row["mu_w"])
        target_temperature = min(expected, key=lambda value: abs(value - temperature))
        if abs(temperature - target_temperature) > 1.0e-13:
            raise RuntimeError(f"unexpected wall temperature {temperature}")
        target = expected[target_temperature]
        if abs(viscosity - target) > 1.0e-13 * target:
            raise RuntimeError(
                f"wall viscosity {viscosity} differs from Sutherland value {target}"
            )
        seen.add(target_temperature)
    if seen != set(expected):
        raise RuntimeError(f"both wall temperatures were not observed: {seen}")
    print(
        "production Sutherland transport passed: "
        f"ranks={args.ranks} rows={len(rows)} mu(T=1)={expected[1.0]:.17g} "
        f"mu(T=2)={expected[2.0]:.17g}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"transport runtime check failed: {error}", file=sys.stderr)
        raise SystemExit(1)
