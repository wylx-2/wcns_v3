#!/usr/bin/env python3
"""Create the retained Stage-Q Case07 review figures."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_case07 import join_zones, plot_field, read_tecplot


def read_history(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    names = lines[0].removeprefix("# ").split()
    numeric_names = names[:-1]
    values = np.asarray(
        [[float(value) for value in line.split()[:-1]] for line in lines[1:] if line.strip()],
        dtype=float,
    )
    return numeric_names, values


def plot_robustness_history(path: Path, output: Path) -> None:
    names, values = read_history(path)
    column = {name: values[:, index] for index, name in enumerate(names)}
    total_faces = sum(column[f"robustness_level{level}_faces"] for level in range(4))
    lower_order_faces = sum(column[f"robustness_level{level}_faces"] for level in range(1, 4))
    ratio = np.divide(
        100.0 * lower_order_faces,
        total_faces,
        out=np.zeros_like(lower_order_faces),
        where=total_faces > 0.0,
    )

    figure, axes = plt.subplots(3, 1, figsize=(9.0, 8.0), sharex=True, constrained_layout=True)
    axes[0].plot(column["time"], ratio, color="tab:red", lw=0.9)
    axes[0].set_ylabel("lower-order\nowner faces (%)")
    axes[0].set_ylim(bottom=0.0)

    axes[1].plot(column["time"], column["troubled_cells"], color="tab:purple", lw=0.9)
    axes[1].set_ylabel("troubled cells")
    axes[1].set_ylim(bottom=0.0)

    valid_dt = np.isfinite(column["accepted_dt"]) & (column["accepted_dt"] > 0.0)
    axes[2].semilogy(
        column["time"][valid_dt],
        column["accepted_dt"][valid_dt],
        color="tab:blue",
        lw=0.9,
    )
    axes[2].set_ylabel("accepted dt")
    axes[2].set_xlabel(r"$tU_\infty/D$")

    for axis in axes:
        axis.grid(True, alpha=0.3)
    figure.suptitle("Stage Q Case07 robustness history (4 MPI ranks)")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    case_root = Path(__file__).resolve().parents[1]
    result_root = case_root / "results" / "mach5-robust-v110"
    final_field = next(result_root.glob("*.field.step00012114.time8p000000000eP00.dat"))
    history = next(result_root.glob("*.history.r4.txt"))
    variables, zones = read_tecplot(final_field)
    fields = join_zones(variables, zones)
    figure_root = case_root / "figures"
    figure_root.mkdir(parents=True, exist_ok=True)

    plot_field(
        fields,
        "rho",
        figure_root / "mach5-robust-density.png",
        "Stage Q Mach 5 cylinder: density",
        (-4.0, 3.0, -3.0, 3.0),
    )
    plot_field(
        fields,
        "mach",
        figure_root / "mach5-robust-mach.png",
        "Stage Q Mach 5 cylinder: Mach number",
        (-4.0, 3.0, -3.0, 3.0),
    )
    plot_robustness_history(history, figure_root / "mach5-robust-history.png")

    if not all(math.isfinite(float(np.nanmin(fields[name]))) for name in ("rho", "p", "T", "mach")):
        raise ValueError("non-finite retained Case07 field")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
