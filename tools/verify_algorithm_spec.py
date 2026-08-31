#!/usr/bin/env python3
"""Exact, dependency-free checks for the frozen Stage-G WCNS formulas."""

from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
from typing import Iterable, Sequence


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def exact_moments(
    name: str,
    nodes: Sequence[F],
    weights: Sequence[F],
    target: F,
    max_degree: int,
    derivative: int = 0,
) -> None:
    require(len(nodes) == len(weights), f"{name}: node/weight size mismatch")
    for degree in range(max_degree + 1):
        actual = sum((w * x**degree for x, w in zip(nodes, weights)), F(0))
        if derivative == 0:
            expected = target**degree
        elif derivative == 1:
            expected = F(0) if degree == 0 else degree * target ** (degree - 1)
        else:
            raise ValueError("only interpolation and first derivatives are checked")
        require(actual == expected, f"{name}: degree {degree}: {actual} != {expected}")


def solve_full_column_rank(matrix: Sequence[Sequence[F]], rhs: Sequence[F]) -> list[F]:
    """Solve a possibly overdetermined, consistent exact system."""
    rows = len(matrix)
    cols = len(matrix[0])
    augmented = [list(matrix[r]) + [rhs[r]] for r in range(rows)]
    pivot_row = 0
    pivot_columns: list[int] = []

    for column in range(cols):
        candidate = next(
            (row for row in range(pivot_row, rows) if augmented[row][column] != 0),
            None,
        )
        if candidate is None:
            continue
        augmented[pivot_row], augmented[candidate] = (
            augmented[candidate],
            augmented[pivot_row],
        )
        pivot = augmented[pivot_row][column]
        augmented[pivot_row] = [value / pivot for value in augmented[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = augmented[row][column]
            if factor != 0:
                augmented[row] = [
                    value - factor * base
                    for value, base in zip(augmented[row], augmented[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1

    for row in range(pivot_row, rows):
        require(
            any(augmented[row][column] != 0 for column in range(cols))
            or augmented[row][-1] == 0,
            "inconsistent exact linear system",
        )
    require(len(pivot_columns) == cols, "linear system does not have a unique solution")

    result = [F(0)] * cols
    for row, column in enumerate(pivot_columns):
        result[column] = augmented[row][-1]
    return result


def check_interpolation_and_derivatives() -> None:
    half_nodes_6 = [F(-5, 2), F(-3, 2), F(-1, 2), F(1, 2), F(3, 2), F(5, 2)]
    i6 = [F(3, 256), F(-25, 256), F(150, 256), F(150, 256), F(-25, 256), F(3, 256)]
    exact_moments("I6 interior", half_nodes_6, i6, F(0), 5)

    centers_5 = [F(1, 2), F(3, 2), F(5, 2), F(7, 2), F(9, 2)]
    for name, target, numerators in (
        ("I6 low face 1/2", F(0), [315, -420, 378, -180, 35]),
        ("I6 low face 3/2", F(1), [35, 140, -70, 28, -5]),
        ("I6 low face 5/2", F(2), [-5, 60, 90, -20, 3]),
    ):
        exact_moments(name, centers_5, [F(value, 128) for value in numerators], target, 4)

    vertices_6 = [F(index) for index in range(6)]
    exact_moments("I6 vertex-to-center interior", [F(index) for index in range(-2, 4)], i6, F(1, 2), 5)
    exact_moments(
        "I6 vertex-to-center first",
        vertices_6,
        [F(value, 256) for value in [63, 315, -210, 126, -45, 7]],
        F(1, 2),
        5,
    )
    exact_moments(
        "I6 vertex-to-center second",
        vertices_6,
        [F(value, 256) for value in [-7, 105, 210, -70, 21, -3]],
        F(3, 2),
        5,
    )

    i4_nodes = [F(-3, 2), F(-1, 2), F(1, 2), F(3, 2)]
    exact_moments("I4 PH interior", i4_nodes, [F(-1, 16), F(9, 16), F(9, 16), F(-1, 16)], F(0), 3)
    centers_4 = [F(1, 2), F(3, 2), F(5, 2), F(7, 2)]
    exact_moments("I4 PH boundary", centers_4, [F(35, 16), F(-35, 16), F(21, 16), F(-5, 16)], F(0), 3)
    exact_moments("I4 PH first internal", centers_4, [F(5, 16), F(15, 16), F(-5, 16), F(1, 16)], F(1), 3)

    d6 = [F(-9, 1920), F(125, 1920), F(-2250, 1920), F(2250, 1920), F(-125, 1920), F(9, 1920)]
    exact_moments("D6 interior", half_nodes_6, d6, F(0), 6, derivative=1)
    face_nodes_5 = [F(index) for index in range(5)]
    exact_moments("D6 low first", face_nodes_5, [F(value, 24) for value in [-22, 17, 9, -5, 1]], F(1, 2), 4, derivative=1)
    exact_moments("D6 low second", face_nodes_5[:4], [F(value, 24) for value in [1, -27, 27, -1]], F(3, 2), 4, derivative=1)

    d4 = [F(1, 24), F(-9, 8), F(9, 8), F(-1, 24)]
    exact_moments("D4 PH interior", i4_nodes, d4, F(0), 4, derivative=1)
    exact_moments("D2 PH boundary", [F(0), F(1)], [F(-1), F(1)], F(1, 2), 2, derivative=1)


def wall_derivative_weights(interior_count: int) -> list[F]:
    nodes = [F(0)] + [F(2 * index + 1, 2) for index in range(interior_count)]
    matrix = [[node**degree for node in nodes] for degree in range(interior_count + 1)]
    rhs = [F(0), F(1)] + [F(0)] * (interior_count - 1)
    return solve_full_column_rank(matrix, rhs)


def check_wall_derivatives() -> None:
    expected_ph = [F(-352, 105), F(35, 8), F(-35, 24), F(21, 40), F(-5, 56)]
    expected_scmm = [
        F(-13016, 3465), F(693, 128), F(-385, 128), F(693, 320),
        F(-495, 448), F(385, 1152), F(-63, 1408),
    ]
    for name, count, expected in (
        ("PH wall derivative", 4, expected_ph),
        ("SCMM6 wall derivative", 6, expected_scmm),
    ):
        weights = wall_derivative_weights(count)
        require(weights == expected, f"{name}: generated coefficients changed")
        nodes = [F(0)] + [F(2 * index + 1, 2) for index in range(count)]
        exact_moments(name, nodes, weights, F(0), count, derivative=1)


def ph_derivative_matrix(cell_count: int) -> list[list[F]]:
    matrix = [[F(0) for _ in range(cell_count + 1)] for _ in range(cell_count)]
    matrix[0][0:2] = [F(-1), F(1)]
    matrix[-1][-2:] = [F(-1), F(1)]
    for cell in range(1, cell_count - 1):
        for face, weight in zip(
            [cell - 1, cell, cell + 1, cell + 2],
            [F(1, 24), F(-9, 8), F(9, 8), F(-1, 24)],
        ):
            matrix[cell][face] = weight
    return matrix


def scmm_derivative_matrix(cell_count: int) -> list[list[F]]:
    matrix = [[F(0) for _ in range(cell_count + 1)] for _ in range(cell_count)]
    for face, value in enumerate([-22, 17, 9, -5, 1]):
        matrix[0][face] = F(value, 24)
    for face, value in enumerate([1, -27, 27, -1]):
        matrix[1][face] = F(value, 24)
    for cell in range(2, cell_count - 2):
        for face, value in zip(range(cell - 2, cell + 4), [-9, 125, -2250, 2250, -125, 9]):
            matrix[cell][face] = F(value, 1920)
    for face, value in zip(range(cell_count - 3, cell_count + 1), [1, -27, 27, -1]):
        matrix[-2][face] = F(value, 24)
    for face, value in zip(range(cell_count - 4, cell_count + 1), [-1, 5, -9, -17, 22]):
        matrix[-1][face] = F(value, 24)
    return matrix


def check_conservation_weights() -> None:
    for name, builder, minimum in (
        ("PH", ph_derivative_matrix, 4),
        ("SCMM6", scmm_derivative_matrix, 5),
    ):
        for cell_count in range(minimum, 13):
            derivative = builder(cell_count)
            transpose = [
                [derivative[cell][face] for cell in range(cell_count)]
                for face in range(cell_count + 1)
            ]
            boundary = [F(-1)] + [F(0)] * (cell_count - 1) + [F(1)]
            weights = solve_full_column_rank(transpose, boundary)
            require(all(weight > 0 for weight in weights), f"{name} N={cell_count}: non-positive integration weight")
            for face in range(cell_count + 1):
                actual = sum(derivative[cell][face] * weights[cell] for cell in range(cell_count))
                require(actual == boundary[face], f"{name} N={cell_count}: D^T w mismatch at face {face}")


Vector = tuple[F, F, F]


def dot(left: Vector, right: Vector) -> F:
    return sum((a * b for a, b in zip(left, right)), F(0))


def cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def check_affine_metrics() -> None:
    r_xi: Vector = (F(2), F(1), F(0))
    r_eta: Vector = (F(-1), F(3), F(1))
    r_zeta: Vector = (F(1), F(0), F(2))
    s_xi = cross(r_eta, r_zeta)
    s_eta = cross(r_zeta, r_xi)
    s_zeta = cross(r_xi, r_eta)
    jacobian = dot(r_xi, s_xi)
    require(jacobian > 0, "affine test mapping must have positive orientation")
    for row, tangent in enumerate((r_xi, r_eta, r_zeta)):
        for column, cofactor in enumerate((s_xi, s_eta, s_zeta)):
            expected = jacobian if row == column else F(0)
            require(dot(tangent, cofactor) == expected, "affine cofactor/Jacobian identity failed")


def check_nondimensionalization() -> None:
    gamma = F(7, 5)
    gas_constant = F(287)
    temperature_ref = F(300)
    velocity_ref = F(340)
    density_ref = F(6, 5)
    length_ref = F(3, 2)
    viscosity_ref = F(9, 500000)
    mach_squared = velocity_ref**2 / (gamma * gas_constant * temperature_ref)
    reynolds = density_ref * velocity_ref * length_ref / viscosity_ref
    require(reynolds > 0 and mach_squared > 0, "reference Re/Ma must be positive")

    rho = F(5, 4)
    temperature = F(7, 6)
    pressure_dimensional = (rho * density_ref) * gas_constant * (temperature * temperature_ref)
    pressure = pressure_dimensional / (density_ref * velocity_ref**2)
    require(pressure == rho * temperature / (gamma * mach_squared), "dimensionless ideal-gas pressure mismatch")

    internal_energy_dimensional = gas_constant * (temperature * temperature_ref) / (gamma - 1)
    internal_energy = internal_energy_dimensional / velocity_ref**2
    require(internal_energy == temperature / (gamma * (gamma - 1) * mach_squared), "dimensionless internal energy mismatch")

    mu_at_reference_temperature = F(11, 10) * viscosity_ref
    sutherland_ratio_at_one = mu_at_reference_temperature / viscosity_ref
    require(sutherland_ratio_at_one == F(11, 10), "Sutherland reference viscosity scale mismatch")

    velocity = (F(2, 5), F(-1, 4), F(1, 10))
    acceleration_dimensional = (F(3), F(-2), F(1, 2))
    acceleration = tuple(
        length_ref * component / velocity_ref**2
        for component in acceleration_dimensional
    )
    density_dimensional = rho * density_ref
    momentum_source_dimensional = tuple(
        density_dimensional * component for component in acceleration_dimensional
    )
    momentum_source = tuple(
        length_ref * component / (density_ref * velocity_ref**2)
        for component in momentum_source_dimensional
    )
    require(
        momentum_source == tuple(rho * component for component in acceleration),
        "dimensionless body-force momentum source mismatch",
    )
    velocity_dimensional = tuple(component * velocity_ref for component in velocity)
    energy_source_dimensional = density_dimensional * sum(
        (u * a for u, a in zip(velocity_dimensional, acceleration_dimensional)),
        F(0),
    )
    energy_source = length_ref * energy_source_dimensional / (density_ref * velocity_ref**3)
    require(
        energy_source
        == rho * sum((u * a for u, a in zip(velocity, acceleration)), F(0)),
        "dimensionless body-force energy source mismatch",
    )


def check_document_contracts() -> None:
    document = Path(__file__).resolve().parents[1] / "算法补充.md"
    text = document.read_text(encoding="utf-8")
    required_fragments: Iterable[str] = (
        "algorithm_profile = phenglei_wcns | scmm6_wcns",
        "禁止把前一套的度量、插值、差分或边界闭合与后一套交叉使用",
        "strong_boundary_face_state = true | false",
        "PhysicalGhostStateOperator",
        "InviscidBoundaryFaceState",
        "ViscousBoundaryTrace",
        "FaceFluxHalo",
        "GeometryOperand",
        "不设计、不实现通量分裂接口",
        "enable_source_terms = false | true",
        "SourceTermOperator",
        "局部源项不需要连接通信",
    )
    for fragment in required_fragments:
        require(fragment in text, f"frozen document contract is missing: {fragment}")


def main() -> None:
    checks = (
        ("interpolation and derivative moments", check_interpolation_and_derivatives),
        ("wall-normal derivative generation", check_wall_derivatives),
        ("positive conservation weights", check_conservation_weights),
        ("affine metric cofactors", check_affine_metrics),
        ("nondimensional identities", check_nondimensionalization),
        ("document contracts", check_document_contracts),
    )
    for label, check in checks:
        check()
        print(f"PASS: {label}")
    print(f"PASS: all {len(checks)} Stage-G specification checks")


if __name__ == "__main__":
    main()
