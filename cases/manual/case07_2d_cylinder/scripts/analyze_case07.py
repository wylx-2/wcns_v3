#!/usr/bin/env python3
"""Post-process Case07 Tecplot cell-centred fields.

The aerodynamic coefficients are diagnostic estimates: pressure and tangential
velocity are extrapolated from the first off-wall cell of the generated circular
O-grid.  They are suitable for qualitative case checks, not a grid-converged
surface-force validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class CaseInfo:
    directory: str
    reynolds: float | None
    mach: float
    viscous: bool


CASES = {
    "re20": CaseInfo("re20", 20.0, 0.2, True),
    "re40": CaseInfo("re40", 40.0, 0.2, True),
    "re100": CaseInfo("re100", 100.0, 0.2, True),
    "re200": CaseInfo("re200", 200.0, 0.2, True),
    "mach5-euler": CaseInfo("mach5-euler", None, 5.0, False),
}

VARIABLE_RE = re.compile(r'"([^"]+)"')
ZONE_RE = re.compile(r'ZONE T="([^"]+)", I=(\d+), J=(\d+)')
TIME_RE = re.compile(r"time([0-9pPeEM]+)\.dat$")


def filename_time(path: Path) -> float:
    match = TIME_RE.search(path.name)
    if not match:
        raise ValueError(f"cannot extract time from {path.name}")
    return float(match.group(1).replace("p", ".").replace("P", "+").replace("M", "-"))


def read_tecplot(path: Path) -> tuple[list[str], list[tuple[str, np.ndarray]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3 or not lines[1].startswith("VARIABLES="):
        raise ValueError(f"unexpected Tecplot header in {path}")
    variables = VARIABLE_RE.findall(lines[1])
    zones: list[tuple[str, np.ndarray]] = []
    line = 2
    while line < len(lines):
        match = ZONE_RE.match(lines[line])
        if not match:
            raise ValueError(f"unexpected zone line {line + 1} in {path}")
        name, ni_text, nj_text = match.groups()
        ni, nj = int(ni_text), int(nj_text)
        line += 1
        count = ni * nj
        values = np.asarray(
            [[float(value) for value in row.split()] for row in lines[line : line + count]],
            dtype=float,
        )
        if values.shape != (count, len(variables)):
            raise ValueError(f"truncated zone {name} in {path}")
        zones.append((name, values.reshape(nj, ni, len(variables))))
        line += count
    return variables, zones


def join_zones(
    variables: list[str], zones: list[tuple[str, np.ndarray]]
) -> dict[str, np.ndarray]:
    ordered = sorted(zones, key=lambda item: int(re.search(r"(\d+)$", item[0]).group(1)))
    joined = np.concatenate([values for _, values in ordered], axis=1)
    return {name: joined[:, :, index] for index, name in enumerate(variables)}


def polar_vorticity(fields: dict[str, np.ndarray]) -> np.ndarray:
    x, y, u, v = fields["X"], fields["Y"], fields["u"], fields["v"]
    radial = np.mean(np.hypot(x, y), axis=1)
    ntheta = x.shape[1]
    theta = -2.0 * math.pi * (np.arange(ntheta) + 0.5) / ntheta
    cosine, sine = np.cos(theta)[None, :], np.sin(theta)[None, :]
    radial_velocity = u * cosine + v * sine
    tangential_velocity = -u * sine + v * cosine
    radial_term = np.gradient(
        radial[:, None] * tangential_velocity, radial, axis=0, edge_order=2
    ) / radial[:, None]
    delta_theta = -2.0 * math.pi / ntheta
    angular_term = (
        np.roll(radial_velocity, -1, axis=1)
        - np.roll(radial_velocity, 1, axis=1)
    ) / (2.0 * delta_theta * radial[:, None])
    return radial_term - angular_term


def surface_diagnostics(
    fields: dict[str, np.ndarray], info: CaseInfo
) -> dict[str, float | list[float]]:
    ntheta = fields["X"].shape[1]
    theta = -2.0 * math.pi * (np.arange(ntheta) + 0.5) / ntheta
    delta_theta = 2.0 * math.pi / ntheta
    radius = 0.5
    diameter = 1.0
    dynamic_pressure = 0.5
    x, y = fields["X"][0], fields["Y"][0]
    p = fields["p"][0]
    u, v = fields["u"][0], fields["v"][0]
    cell_radius = np.hypot(x, y)
    cosine, sine = np.cos(theta), np.sin(theta)
    pressure_inf = 1.0 / (1.4 * info.mach * info.mach)
    pressure_force_x = -(p - pressure_inf) * cosine * radius * delta_theta
    pressure_force_y = -(p - pressure_inf) * sine * radius * delta_theta
    viscous_force_x = np.zeros_like(p)
    viscous_force_y = np.zeros_like(p)
    if info.viscous and info.reynolds is not None:
        tangent_x, tangent_y = -sine, cosine
        tangent_velocity = u * tangent_x + v * tangent_y
        viscosity_ratio = fields.get("viscosity", np.ones_like(fields["p"]))[0]
        shear = viscosity_ratio * tangent_velocity / (
            info.reynolds * (cell_radius - radius)
        )
        viscous_force_x = shear * tangent_x * radius * delta_theta
        viscous_force_y = shear * tangent_y * radius * delta_theta
    scale = dynamic_pressure * diameter
    return {
        "cd_pressure": float(np.sum(pressure_force_x) / scale),
        "cl_pressure": float(np.sum(pressure_force_y) / scale),
        "cd_viscous_estimate": float(np.sum(viscous_force_x) / scale),
        "cl_viscous_estimate": float(np.sum(viscous_force_y) / scale),
        "cd_estimate": float(np.sum(pressure_force_x + viscous_force_x) / scale),
        "cl_estimate": float(np.sum(pressure_force_y + viscous_force_y) / scale),
        "theta": theta.tolist(),
        "x": x.tolist(),
        "y": y.tolist(),
        "cp_first_cell": ((p - pressure_inf) / dynamic_pressure).tolist(),
    }


def flow_diagnostics(fields: dict[str, np.ndarray]) -> dict[str, float]:
    x, y, u, v = fields["X"], fields["Y"], fields["u"], fields["v"]
    wake = (x > 0.5) & (np.abs(y) < 1.0)
    reversed_wake = wake & (u < 0.0)
    vorticity = polar_vorticity(fields)
    result = {
        "rho_min": float(np.min(fields["rho"])),
        "rho_max": float(np.max(fields["rho"])),
        "p_min": float(np.min(fields["p"])),
        "p_max": float(np.max(fields["p"])),
        "mach_min": float(np.min(fields["mach"])),
        "mach_max": float(np.max(fields["mach"])),
        "wake_u_min": float(np.min(u[wake])),
        "wake_abs_v_max": float(np.max(np.abs(v[wake]))),
        "vorticity_abs_max": float(np.max(np.abs(vorticity))),
        "reverse_flow_cell_count": int(np.count_nonzero(reversed_wake)),
        "reverse_flow_x_max": (
            float(np.max(x[reversed_wake])) if np.any(reversed_wake) else float("nan")
        ),
    }
    return result


def bow_shock_diagnostics(fields: dict[str, np.ndarray]) -> dict[str, float]:
    x, y, rho = fields["X"], fields["Y"], fields["rho"]
    upstream_columns = np.where(np.mean(x, axis=0) < 0.0)[0]
    column = min(upstream_columns, key=lambda i: abs(float(y[0, i])))
    order = np.argsort(x[:, column])
    line_x = x[order, column]
    line_rho = rho[order, column]
    valid = line_x < -0.55
    gradient = np.abs(np.gradient(line_rho[valid], line_x[valid]))
    shock_x = float(line_x[valid][int(np.argmax(gradient))])
    return {
        "bow_shock_x_estimate": shock_x,
        "bow_shock_standoff_estimate": -0.5 - shock_x,
        "upstream_line_rho_max": float(np.max(line_rho[valid])),
    }


def plot_field(
    fields: dict[str, np.ndarray], quantity: str,
    output: Path, title: str, limits: tuple[float, float, float, float]
) -> None:
    figure, axis = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    values = fields[quantity]
    levels = np.linspace(float(np.nanpercentile(values, 1)), float(np.nanpercentile(values, 99)), 41)
    closed_x = np.concatenate([fields["X"], fields["X"][:, :1]], axis=1)
    closed_y = np.concatenate([fields["Y"], fields["Y"][:, :1]], axis=1)
    closed_values = np.concatenate([values, values[:, :1]], axis=1)
    contour = axis.contourf(
        closed_x, closed_y, closed_values,
        levels=levels, extend="both", cmap="turbo"
    )
    circle = plt.Circle((0.0, 0.0), 0.5, color="white", ec="black", lw=1.0, zorder=10)
    axis.add_patch(circle)
    axis.set_aspect("equal")
    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_xlabel("x/D")
    axis.set_ylabel("y/D")
    axis.set_title(title)
    figure.colorbar(contour, ax=axis, label=quantity)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_derived_vorticity(
    fields: dict[str, np.ndarray], output: Path, title: str
) -> None:
    vorticity = polar_vorticity(fields)
    limit = float(np.nanpercentile(np.abs(vorticity), 98))
    figure, axis = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    contour = axis.tricontourf(
        fields["X"].ravel(), fields["Y"].ravel(), vorticity.ravel(),
        levels=np.linspace(-limit, limit, 41), cmap="RdBu_r", extend="both"
    )
    axis.add_patch(plt.Circle((0.0, 0.0), 0.5, color="white", ec="black", zorder=10))
    axis.set_aspect("equal")
    axis.set_xlim(-2.0, 5.0)
    axis.set_ylim(-2.0, 2.0)
    axis.set_xlabel("x/D")
    axis.set_ylabel("y/D")
    axis.set_title(title)
    figure.colorbar(contour, ax=axis, label=r"$\omega_z D/U_\infty$")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def shedding_frequency(rows: list[dict[str, float]]) -> dict[str, float]:
    if len(rows) < 20:
        return {"strouhal_estimate": float("nan"), "cl_rms_late": float("nan")}
    start = len(rows) // 2
    times = np.asarray([row["time"] for row in rows[start:]])
    def estimate(name: str) -> tuple[float, float]:
        signal = np.asarray([row[name] for row in rows[start:]])
        centered = signal - np.mean(signal)
        shifted_time = times - times[0]
        best_frequency = float("nan")
        best_error = float("inf")
        # A sinusoid plus constant/linear drift avoids the coarse FFT-bin limit
        # of the deliberately short qualitative runs.
        for frequency in np.linspace(0.05, 0.5, 4501):
            phase = 2.0 * math.pi * frequency * shifted_time
            design = np.column_stack((
                np.sin(phase), np.cos(phase),
                np.ones_like(phase), shifted_time,
            ))
            coefficients, *_ = np.linalg.lstsq(design, signal, rcond=None)
            error = float(np.mean((design @ coefficients - signal) ** 2))
            if error < best_error:
                best_error = error
                best_frequency = float(frequency)
        return best_frequency, float(np.sqrt(np.mean(centered * centered)))

    lift_frequency, lift_rms = estimate("cl_estimate")
    probe_frequency, probe_rms = estimate("wake_probe_v")
    return {
        "strouhal_from_lift_estimate": lift_frequency,
        "strouhal_from_probe_estimate": probe_frequency,
        "cl_rms_late": lift_rms,
        "wake_probe_v_rms_late": probe_rms,
    }


def read_solver_history(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    names = lines[0].removeprefix("# ").split()
    rows = [line.split() for line in lines[1:] if line.strip()]
    last = dict(zip(names, rows[-1]))
    return {
        "step": int(last["step"]),
        "solver_time": float(last["time"]),
        "last_dt": float(last["dt"]),
        "wall_time_seconds": float(last["wall_time"]),
        "final_total_l2": float(last["total_l2"]),
        "reconstruction_fallbacks": int(last["reconstruction_fallbacks"]),
        "riemann_fallbacks": int(last["riemann_fallbacks"]),
        "stop_reason": last["stop_reason"],
    }


def analyze_case(root: Path, key: str, info: CaseInfo) -> dict[str, object]:
    result_directory = root / "results" / info.directory
    files = sorted(result_directory.glob("*.field.*.dat"), key=filename_time)
    if not files:
        raise FileNotFoundError(f"no field files in {result_directory}")
    history_rows: list[dict[str, float]] = []
    final_variables: list[str] = []
    final_zones: list[tuple[str, np.ndarray]] = []
    final_fields: dict[str, np.ndarray] = {}
    final_surface: dict[str, float | list[float]] = {}
    for path in files:
        variables, zones = read_tecplot(path)
        fields = join_zones(variables, zones)
        surface = surface_diagnostics(fields, info)
        probe_distance = (fields["X"] - 2.0) ** 2 + (fields["Y"] - 0.5) ** 2
        probe = np.unravel_index(int(np.argmin(probe_distance)), probe_distance.shape)
        history_rows.append({
            "time": filename_time(path),
            "cd_pressure": float(surface["cd_pressure"]),
            "cd_viscous_estimate": float(surface["cd_viscous_estimate"]),
            "cd_estimate": float(surface["cd_estimate"]),
            "cl_estimate": float(surface["cl_estimate"]),
            "wake_probe_v": float(fields["v"][probe]),
        })
        final_variables, final_zones, final_fields, final_surface = (
            variables, zones, fields, surface
        )

    with (root / "results" / f"{key}.surface-history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(history_rows[0]))
        writer.writeheader()
        writer.writerows(history_rows)
    with (root / "results" / f"{key}.surface-final.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["theta", "x", "y", "cp_first_cell"])
        writer.writerows(zip(
            final_surface["theta"], final_surface["x"],
            final_surface["y"], final_surface["cp_first_cell"]
        ))

    final = {
        "last_field": str(files[-1].relative_to(root)),
        "time": filename_time(files[-1]),
        **flow_diagnostics(final_fields),
        **{name: value for name, value in final_surface.items() if not isinstance(value, list)},
    }
    solver_histories = list(result_directory.glob("*.history.r*.txt"))
    if len(solver_histories) != 1:
        raise ValueError(f"expected one completed solver history in {result_directory}")
    final.update(read_solver_history(solver_histories[0]))
    if key in {"re100", "re200"}:
        final.update(shedding_frequency(history_rows))
    if key == "mach5-euler":
        final.update(bow_shock_diagnostics(final_fields))

    if key == "mach5-euler":
        plot_field(
            final_fields, "rho", root / "figures" / "mach5-density.png",
            "Mach 5 cylinder: density", (-4.0, 3.0, -3.0, 3.0)
        )
        plot_field(
            final_fields, "mach", root / "figures" / "mach5-mach.png",
            "Mach 5 cylinder: Mach number", (-4.0, 3.0, -3.0, 3.0)
        )
    else:
        plot_derived_vorticity(
            final_fields, root / "figures" / f"{key}-vorticity.png",
            f"Cylinder {key}: vorticity"
        )
    return final


def plot_histories(root: Path) -> None:
    figure, axis = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    for key, info in CASES.items():
        path = next((root / "results" / info.directory).glob("*.history.r*.txt"))
        lines = path.read_text(encoding="utf-8").splitlines()
        names = lines[0].removeprefix("# ").split()
        table = np.asarray([
            [float(value) if value not in {"running", "physical_time_reached"} else np.nan
             for value in line.split()]
            for line in lines[1:] if line.strip()
        ])
        axis.semilogy(table[:, names.index("time")], table[:, names.index("total_l2")], label=key)
    axis.set_xlabel(r"$tU_\infty/D$")
    axis.set_ylabel("global residual L2")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend()
    figure.savefig(root / "figures" / "residual-history.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=False, constrained_layout=True)
    for key in ("re20", "re40", "re100", "re200"):
        with (root / "results" / f"{key}.surface-history.csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        time = np.asarray([float(row["time"]) for row in rows])
        drag = np.asarray([float(row["cd_estimate"]) for row in rows])
        axes[0].plot(time, drag, marker="o", ms=2.0, label=key)
        if key in {"re100", "re200"}:
            lift = np.asarray([float(row["cl_estimate"]) for row in rows])
            probe = np.asarray([float(row["wake_probe_v"]) for row in rows])
            axes[1].plot(time, lift, label=rf"{key} estimated $C_L$")
            axes[1].plot(time, probe, label=rf"{key} probe $v(2D,0.5D)$")
    axes[0].set_xlabel(r"$tU_\infty/D$")
    axes[0].set_ylabel("estimated drag coefficient")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].set_xlabel(r"$tU_\infty/D$")
    axes[1].set_ylabel("oscillatory signals")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    figure.savefig(root / "figures" / "force-and-probe-history.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    for key in CASES:
        with (root / "results" / f"{key}.surface-final.csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        angle = np.mod(np.degrees([float(row["theta"]) for row in rows]), 360.0)
        cp = np.asarray([float(row["cp_first_cell"]) for row in rows])
        order = np.argsort(angle)
        axis.plot(angle[order], cp[order], marker="o", ms=2.5, label=key)
    axis.set_xlabel("polar angle from downstream axis (degree)")
    axis.set_ylabel(r"first-cell $C_p$")
    axis.grid(True, alpha=0.3)
    axis.legend(ncol=2)
    figure.savefig(root / "figures" / "surface-cp.png", dpi=180)
    plt.close(figure)


def plot_grid_from_initial(root: Path) -> None:
    sources = (
        ("re20", "32x20 low-speed grid"),
        ("mach5-euler", "48x32 Mach 5 grid"),
    )
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 5.0), constrained_layout=True)
    for axis, (key, title) in zip(axes, sources):
        info = CASES[key]
        path = min(
            (root / "results" / info.directory).glob("*.field.*.dat"),
            key=filename_time,
        )
        variables, zones = read_tecplot(path)
        fields = join_zones(variables, zones)
        x = np.concatenate([fields["X"], fields["X"][:, :1]], axis=1)
        y = np.concatenate([fields["Y"], fields["Y"][:, :1]], axis=1)
        for radial_index in range(x.shape[0]):
            axis.plot(x[radial_index], y[radial_index], color="0.45", lw=0.35)
        for angular_index in range(0, x.shape[1] - 1, max(1, (x.shape[1] - 1) // 16)):
            axis.plot(x[:, angular_index], y[:, angular_index], color="0.45", lw=0.35)
        axis.add_patch(plt.Circle((0.0, 0.0), 0.5, color="white", ec="black"))
        axis.set_aspect("equal")
        axis.set_xlim(-2.0, 2.0)
        axis.set_ylim(-2.0, 2.0)
        axis.set_xlabel("x/D")
        axis.set_ylabel("y/D")
        axis.set_title(title)
    figure.savefig(root / "figures" / "grid-near-cylinder.png", dpi=180)
    plt.close(figure)


def strict_json_value(value: object) -> object:
    """Replace non-finite floats so the report is valid RFC 8259 JSON."""
    if isinstance(value, dict):
        return {key: strict_json_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [strict_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case-root", type=Path,
        default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    root = args.case_root.resolve()
    summary: dict[str, object] = {}
    for key, info in CASES.items():
        summary[key] = analyze_case(root, key, info)
    plot_histories(root)
    plot_grid_from_initial(root)
    serializable_summary = strict_json_value(summary)
    output = root / "results" / "analysis-summary.json"
    output.write_text(
        json.dumps(serializable_summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(json.dumps(serializable_summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
