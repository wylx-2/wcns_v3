#pragma once

#include <wcns/solver/physical_boundary.hpp>
#include <wcns/solver/viscous_gradient.hpp>

#include <vector>

namespace wcns {

[[nodiscard]] Real wall_dirichlet_computational_derivative(
    Real wall_value,
    const std::vector<Real>& inward_center_values,
    const AlgorithmProfile& profile);

[[nodiscard]] Real interpolate_internal_pressure_trace(
    const StructuredBlock& block,
    const AlgorithmProfile& profile,
    Axis axis,
    Index3 face);

[[nodiscard]] ViscousFaceTrace apply_viscous_boundary_trace(
    const StructuredBlock& block,
    const MetricField& metric,
    const BoundaryPatch& patch,
    Index3 face,
    const ViscousFaceTrace& raw_trace,
    Real internal_pressure_trace,
    Normal3 outward_unit_normal,
    const BoundaryData& data,
    const AlgorithmProfile& profile,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors);

} // namespace wcns
