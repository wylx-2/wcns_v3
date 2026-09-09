#!/usr/bin/env python3
"""Validate the short case05 statistics file without claiming turbulent convergence."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


REQUIRED = {
    "step", "time", "channel_wall_shear_lower", "channel_wall_shear_upper",
    "channel_wall_shear_mean", "channel_friction_velocity", "channel_re_tau",
    "xz_mean_u_j0", "xz_mass_flow_x_j0", "xz_mean_u_j23",
    "xz_mass_flow_x_j23", "xz_mean_u_j47", "xz_mass_flow_x_j47",
}


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_case05.py <statistics.txt>")
    path = Path(sys.argv[1])
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip()]
    if len(lines) < 3:
        raise RuntimeError("statistics output lacks initial/final records")
    header = lines[0][2:] if lines[0].startswith("# ") else lines[0]
    names = header.split()
    missing = sorted(REQUIRED.difference(names))
    if missing:
        raise RuntimeError(f"statistics output lacks columns: {missing}")
    records: list[dict[str, float]] = []
    for line in lines[1:]:
        values = [float(value) for value in line.split()]
        if len(values) != len(names) or not all(math.isfinite(v) for v in values):
            raise RuntimeError("statistics output contains a malformed row")
        records.append(dict(zip(names, values)))
    initial, final = records[0], records[-1]
    if initial["step"] != 0 or final["step"] != 5:
        raise RuntimeError("feasibility run did not complete the configured five steps")
    if not 0.0 < final["time"] < 1.0e-2:
        raise RuntimeError(f"unexpected five-step final time: {final['time']}")
    for record in (initial, final):
        for name in (
            "channel_wall_shear_lower", "channel_wall_shear_upper",
            "channel_wall_shear_mean", "channel_friction_velocity",
            "channel_re_tau",
        ):
            if record[name] <= 0.0:
                raise RuntimeError(f"non-positive wall statistic: {name}")
    if not 150.0 <= initial["channel_re_tau"] <= 210.0:
        raise RuntimeError(
            f"initial measured Re_tau is inconsistent with the target: "
            f"{initial['channel_re_tau']}")
    case_directory = Path(__file__).resolve().parent
    try:
        displayed_path = str(path.resolve().relative_to(case_directory))
    except ValueError:
        displayed_path = str(path)
    summary = {
        "status": "short_feasibility_passed_not_turbulence_accepted",
        "statistics_file": displayed_path.replace("\\", "/"),
        "initial_re_tau": initial["channel_re_tau"],
        "final_re_tau": final["channel_re_tau"],
        "initial_wall_shear": initial["channel_wall_shear_mean"],
        "final_wall_shear": final["channel_wall_shear_mean"],
        "final_step": int(final["step"]),
        "final_time": final["time"],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
