#!/usr/bin/env python3
"""Validate the production x-z plane statistic names and uniform-field values."""

from __future__ import annotations

import math
import sys
from pathlib import Path


def main() -> int:
    path = Path(sys.argv[1])
    lines = path.read_text(encoding="utf-8").splitlines()
    expected = [
        "step",
        "time",
        "xz_mean_u_j0",
        "xz_mass_flow_x_j0",
        "xz_mean_u_j7",
        "xz_mass_flow_x_j7",
    ]
    if not lines or lines[0].removeprefix("# ").split() != expected:
        raise RuntimeError("x-z plane statistic header is incorrect")
    if len(lines) < 3:
        raise RuntimeError("x-z plane statistic output lacks initial/final rows")
    for line in lines[1:]:
        values = [float(value) for value in line.split()]
        if len(values) != len(expected) or not all(math.isfinite(v) for v in values):
            raise RuntimeError("x-z plane statistic row is malformed")
        for value in values[2:]:
            if abs(value - 0.2) > 2.0e-12:
                raise RuntimeError(f"uniform x-z plane statistic is {value}, expected 0.2")
    print(f"validated x-z plane statistics: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
