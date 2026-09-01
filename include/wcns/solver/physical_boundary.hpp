#pragma once

#include <wcns/mesh/structured_block.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>

namespace wcns {

struct BoundaryData {
    std::optional<TemperaturePrimitiveState> target_state;
    std::array<Real, 3> wall_velocity {{0.0, 0.0, 0.0}};
    std::optional<Real> wall_temperature;

    void validate(BoundaryType type, int dimension) const;
};

using BoundaryDataMap = std::unordered_map<std::string, BoundaryData>;

struct PhysicalGhostFillResult {
    std::uint64_t version = 0;
    std::size_t state_count = 0;
};

class PhysicalGhostStateOperator {
public:
    [[nodiscard]] static PhysicalGhostFillResult fill(
        StructuredBlock& block,
        const BoundaryDataMap& boundary_data,
        const GasModel& gas,
        const ReferenceScales& reference,
        const NumericalFloors& floors,
        std::uint64_t version);
};

struct InviscidBoundaryOptions {
    bool strong_boundary_face_state = true;

    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

[[nodiscard]] PressurePrimitiveState apply_inviscid_boundary_face_state(
    const BoundaryPatch& patch,
    const PressurePrimitiveState& interior_trace,
    const PressurePrimitiveState& reconstructed_exterior_trace,
    Normal3 outward_unit_normal,
    const BoundaryData& data,
    const InviscidBoundaryOptions& options,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension);

void update_temperature_primitive_interior(
    StructuredBlock& block,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors);

void update_temperature_primitive_cell(
    StructuredBlock& block,
    Index3 index,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors);

} // namespace wcns
