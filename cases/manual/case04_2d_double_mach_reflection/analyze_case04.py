#!/usr/bin/env python3
"""Compare final Tecplot fields for the two frozen case04 methods."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


CASE_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("weno5", type=Path)
    parser.add_argument("mdcd_hybrid", type=Path)
    parser.add_argument("--weno-history", type=Path)
    parser.add_argument("--mdcd-history", type=Path)
    parser.add_argument("--weno-statistics", type=Path)
    parser.add_argument("--mdcd-statistics", type=Path)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=CASE_DIR / "comparison",
    )
    return parser.parse_args()


def quoted(line: str) -> list[str]:
    return re.findall(r'"([^"]+)"', line)


def extent(line: str, name: str) -> int:
    match = re.search(rf"(?:^|[, ]){name}=(\d+)", line)
    if match is None:
        raise RuntimeError(f"Tecplot zone is missing {name}")
    return int(match.group(1))


def read_tecplot(path: Path) -> dict[str, np.ndarray]:
    variables: list[str] | None = None
    records: list[list[float]] = []
    with path.open("r", encoding="utf-8") as stream:
        while True:
            line = stream.readline()
            if not line:
                break
            stripped = line.strip()
            if stripped.startswith("VARIABLES="):
                variables = quoted(stripped)
                continue
            if not stripped.startswith("ZONE "):
                continue
            if variables is None:
                raise RuntimeError("ZONE precedes VARIABLES")
            count = extent(stripped, "I") * extent(stripped, "J")
            for _ in range(count):
                values = [float(value) for value in stream.readline().split()]
                if len(values) != len(variables):
                    raise RuntimeError("Tecplot row has the wrong column count")
                records.append(values)
    if variables is None or not records:
        raise RuntimeError(f"empty Tecplot field: {path}")
    array = np.asarray(records, dtype=np.float64)
    return {name: array[:, column] for column, name in enumerate(variables)}


def coordinate_axis(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cluster nominally identical centres emitted by independent CGNS zones."""
    tolerance = 1.0e-12 * max(1.0, float(np.max(np.abs(values))))
    centres: list[float] = []
    counts: list[int] = []
    for value in np.sort(values):
        scalar = float(value)
        if not centres or abs(scalar - centres[-1]) > tolerance:
            centres.append(scalar)
            counts.append(1)
        else:
            counts[-1] += 1
            centres[-1] += (scalar - centres[-1]) / counts[-1]
    axis = np.asarray(centres)
    insertion = np.searchsorted(axis, values)
    upper = np.minimum(insertion, axis.size - 1)
    lower = np.maximum(insertion - 1, 0)
    indices = np.where(
        np.abs(values - axis[lower]) <= np.abs(values - axis[upper]),
        lower,
        upper,
    )
    if np.max(np.abs(values - axis[indices])) > tolerance:
        raise RuntimeError("Tecplot coordinate clustering exceeded tolerance")
    return axis, indices


def structured(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    x_values, i = coordinate_axis(data["X"])
    y_values, j = coordinate_axis(data["Y"])
    expected = x_values.size * y_values.size
    if expected != data["X"].size:
        raise RuntimeError("Tecplot cells do not form one rectangular structured grid")
    occupancy = np.zeros((y_values.size, x_values.size), dtype=np.int8)
    np.add.at(occupancy, (j, i), 1)
    if not np.all(occupancy == 1):
        raise RuntimeError("Tecplot cells contain duplicate or missing coordinates")
    result: dict[str, np.ndarray] = {}
    for name, values in data.items():
        if name in {"X", "Y"}:
            continue
        field = np.full((y_values.size, x_values.size), np.nan)
        field[j, i] = values
        if not np.all(np.isfinite(field)):
            raise RuntimeError(f"field {name} is incomplete or non-finite")
        result[name] = field
    return x_values, y_values, result


def field_statistics(
    reference: np.ndarray,
    value: np.ndarray,
) -> dict[str, float]:
    difference = value - reference
    reference_l2 = float(np.sqrt(np.mean(reference * reference)))
    return {
        "weno5_minimum": float(reference.min()),
        "weno5_maximum": float(reference.max()),
        "weno5_mean": float(reference.mean()),
        "mdcd_minimum": float(value.min()),
        "mdcd_maximum": float(value.max()),
        "mdcd_mean": float(value.mean()),
        "difference_l1": float(np.mean(np.abs(difference))),
        "difference_l2": float(np.sqrt(np.mean(difference * difference))),
        "difference_linf": float(np.max(np.abs(difference))),
        "difference_relative_l2": float(
            np.sqrt(np.mean(difference * difference)) / max(reference_l2, np.finfo(float).tiny)
        ),
    }


def read_table(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().lstrip("# ").split()
    values = np.loadtxt(path, comments="#", dtype=str)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    if values.shape[1] != len(header):
        raise RuntimeError(f"table column count does not match header: {path}")
    result: dict[str, np.ndarray] = {}
    for column, name in enumerate(header):
        if name != "stop_reason":
            result[name] = values[:, column].astype(np.float64)
    return result


def history_figure(
    weno: dict[str, np.ndarray],
    mdcd: dict[str, np.ndarray],
    output: Path,
) -> dict[str, dict[str, float | int]]:
    figure, axis = plt.subplots(figsize=(10, 5), constrained_layout=True)
    tiny = np.finfo(float).tiny
    axis.semilogy(weno["time"], np.maximum(weno["total_l2"], tiny), label="WENO5")
    axis.semilogy(
        mdcd["time"],
        np.maximum(mdcd["total_l2"], tiny),
        label="MDCD-HYBRID",
    )
    axis.set_xlabel("time")
    axis.set_ylabel("total residual L2")
    axis.set_title("Non-steady residual histories")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.savefig(output, dpi=180)
    plt.close(figure)

    def summarize(table: dict[str, np.ndarray]) -> dict[str, float | int]:
        reconstruction = table["reconstruction_fallbacks"]
        riemann = table["riemann_fallbacks"]
        return {
            "history_rows": int(table["step"].size),
            "minimum_dt_excluding_initial": float(table["dt"][1:].min()),
            "maximum_dt": float(table["dt"].max()),
            "final_total_l2": float(table["total_l2"][-1]),
            "maximum_reconstruction_fallbacks_per_residual": int(reconstruction.max()),
            "rows_with_reconstruction_fallbacks": int(np.count_nonzero(reconstruction)),
            "sum_of_reported_reconstruction_fallbacks": int(reconstruction.sum()),
            "maximum_riemann_fallbacks_per_residual": int(riemann.max()),
        }

    return {"weno5": summarize(weno), "mdcd_hybrid": summarize(mdcd)}


def statistics_figure(
    weno: dict[str, np.ndarray],
    mdcd: dict[str, np.ndarray],
    output: Path,
) -> dict[str, dict[str, float]]:
    quantities = (
        "total_mass",
        "total_momentum_x",
        "total_momentum_y",
        "total_energy",
    )
    figure, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    summary: dict[str, dict[str, float]] = {}
    for axis, quantity in zip(axes.flat, quantities):
        axis.plot(weno["time"], weno[quantity], label="WENO5")
        axis.plot(mdcd["time"], mdcd[quantity], label="MDCD-HYBRID")
        axis.set_title(quantity)
        axis.set_xlabel("time")
        axis.grid(alpha=0.25)
        axis.legend()
        weno_final = float(weno[quantity][-1])
        mdcd_final = float(mdcd[quantity][-1])
        summary[quantity] = {
            "weno5_final": weno_final,
            "mdcd_final": mdcd_final,
            "mdcd_minus_weno5": mdcd_final - weno_final,
            "relative_difference_to_weno5": (
                (mdcd_final - weno_final) / max(abs(weno_final), np.finfo(float).tiny)
            ),
        }
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return summary


def density_gradient(
    density: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    dy, dx = np.gradient(density, y, x)
    return np.sqrt(dx * dx + dy * dy)


def density_figure(
    x: np.ndarray,
    y: np.ndarray,
    weno: np.ndarray,
    mdcd: np.ndarray,
    output: Path,
) -> None:
    difference = np.abs(mdcd - weno)
    figure, axes = plt.subplots(3, 1, figsize=(14, 8.5), constrained_layout=True)
    minimum = min(float(weno.min()), float(mdcd.min()))
    maximum = max(float(weno.max()), float(mdcd.max()))
    for axis, field, title in (
        (axes[0], weno, "SCMM6 + Roe + WENO5 (WENO-JS)"),
        (axes[1], mdcd, "SCMM6 + Roe + MDCD-HYBRID"),
    ):
        image = axis.imshow(
            field,
            origin="lower",
            extent=(x[0], x[-1], y[0], y[-1]),
            aspect="equal",
            cmap="turbo",
            vmin=minimum,
            vmax=maximum,
            interpolation="nearest",
        )
        axis.set_title(title)
        axis.set_ylabel("y")
        figure.colorbar(image, ax=axis, label="density")
    image = axes[2].imshow(
        difference,
        origin="lower",
        extent=(x[0], x[-1], y[0], y[-1]),
        aspect="equal",
        cmap="magma",
        interpolation="nearest",
    )
    axes[2].set_title("Absolute density difference |MDCD-HYBRID - WENO5|")
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("y")
    figure.colorbar(image, ax=axes[2], label="absolute difference")
    figure.savefig(output, dpi=220)
    plt.close(figure)


def schlieren_figure(
    x: np.ndarray,
    y: np.ndarray,
    weno: np.ndarray,
    mdcd: np.ndarray,
    output: Path,
) -> tuple[float, float]:
    gradients = [density_gradient(weno, x, y), density_gradient(mdcd, x, y)]
    scale = max(float(np.percentile(field, 99.8)) for field in gradients)
    figure, axes = plt.subplots(2, 1, figsize=(14, 5.8), constrained_layout=True)
    for axis, gradient, title in zip(
        axes,
        gradients,
        ("WENO5 density-gradient indicator", "MDCD-HYBRID density-gradient indicator"),
    ):
        normalized = np.log1p(gradient) / np.log1p(scale)
        image = axis.imshow(
            normalized,
            origin="lower",
            extent=(x[0], x[-1], y[0], y[-1]),
            aspect="equal",
            cmap="gray_r",
            vmin=0.0,
            vmax=1.0,
            interpolation="nearest",
        )
        axis.set_title(title)
        axis.set_ylabel("y")
        figure.colorbar(image, ax=axis, label="normalized log(1+|grad rho|)")
    axes[-1].set_xlabel("x")
    figure.savefig(output, dpi=220)
    plt.close(figure)
    return float(gradients[0].max()), float(gradients[1].max())


def section_figure(
    x: np.ndarray,
    y: np.ndarray,
    weno: np.ndarray,
    mdcd: np.ndarray,
    output: Path,
    csv_path: Path,
) -> list[dict[str, float]]:
    targets = (0.2, 0.5, 0.8)
    figure, axes = plt.subplots(3, 1, figsize=(12, 9), constrained_layout=True)
    rows: list[list[float | str]] = []
    summaries: list[dict[str, float]] = []
    for axis, target in zip(axes, targets):
        index = int(np.argmin(np.abs(y - target)))
        actual_y = float(y[index])
        difference = mdcd[index] - weno[index]
        axis.plot(x, weno[index], label="WENO5", linewidth=1.0)
        axis.plot(x, mdcd[index], label="MDCD-HYBRID", linewidth=1.0, alpha=0.85)
        axis.set_title(f"Density at y={actual_y:.9f}")
        axis.set_ylabel("density")
        axis.grid(alpha=0.25)
        axis.legend()
        summaries.append(
            {
                "target_y": target,
                "sample_y": actual_y,
                "difference_l2": float(np.sqrt(np.mean(difference * difference))),
                "difference_linf": float(np.max(np.abs(difference))),
            }
        )
        rows.extend(
            [
                actual_y,
                float(x_value),
                float(weno_value),
                float(mdcd_value),
                float(mdcd_value - weno_value),
            ]
            for x_value, weno_value, mdcd_value in zip(x, weno[index], mdcd[index])
        )
    axes[-1].set_xlabel("x")
    figure.savefig(output, dpi=180)
    plt.close(figure)
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("y", "x", "rho_weno5", "rho_mdcd_hybrid", "difference"))
        writer.writerows(rows)
    return summaries


def main() -> int:
    args = parse_args()
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    weno_x, weno_y, weno = structured(read_tecplot(args.weno5.resolve()))
    mdcd_x, mdcd_y, mdcd = structured(read_tecplot(args.mdcd_hybrid.resolve()))
    if not np.array_equal(weno_x, mdcd_x) or not np.array_equal(weno_y, mdcd_y):
        raise RuntimeError("the two Tecplot files use different cell centres")
    common = sorted(set(weno) & set(mdcd))
    if set(weno) != set(mdcd):
        raise RuntimeError("the two Tecplot files contain different fields")

    statistics = {name: field_statistics(weno[name], mdcd[name]) for name in common}
    density_figure(
        weno_x,
        weno_y,
        weno["rho"],
        mdcd["rho"],
        output / "density-comparison.png",
    )
    gradient_maxima = schlieren_figure(
        weno_x,
        weno_y,
        weno["rho"],
        mdcd["rho"],
        output / "density-gradient-comparison.png",
    )
    sections = section_figure(
        weno_x,
        weno_y,
        weno["rho"],
        mdcd["rho"],
        output / "density-sections.png",
        output / "density-sections.csv",
    )
    optional_paths = (
        args.weno_history,
        args.mdcd_history,
        args.weno_statistics,
        args.mdcd_statistics,
    )
    if any(path is not None for path in optional_paths) and not all(
        path is not None for path in optional_paths
    ):
        raise RuntimeError("all four history/statistics paths must be provided together")
    history_summary: dict[str, dict[str, float | int]] = {}
    integral_summary: dict[str, dict[str, float]] = {}
    if all(path is not None for path in optional_paths):
        history_summary = history_figure(
            read_table(args.weno_history.resolve()),
            read_table(args.mdcd_history.resolve()),
            output / "residual-history-comparison.png",
        )
        integral_summary = statistics_figure(
            read_table(args.weno_statistics.resolve()),
            read_table(args.mdcd_statistics.resolve()),
            output / "integral-statistics-comparison.png",
        )
    summary = {
        "status": "comparison_completed",
        "grid_cells": [int(weno_x.size), int(weno_y.size)],
        "weno5_definition": "weno_js",
        "density_gradient_maximum": {
            "weno5": gradient_maxima[0],
            "mdcd_hybrid": gradient_maxima[1],
        },
        "fields": statistics,
        "density_sections": sections,
        "history": history_summary,
        "integral_statistics": integral_summary,
    }
    (output / "comparison-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
