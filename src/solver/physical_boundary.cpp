#include <wcns/solver/physical_boundary.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace wcns {
namespace {

bool finite(Real value)
{
    return std::isfinite(value);
}

void require_unit_normal(Normal3 normal)
{
    const Real norm = std::sqrt(
        normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!finite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsConfigurationError("boundary face normal must be a unit vector");
    }
}

TemperaturePrimitiveState load_temperature(const Field<Real>& field, Index3 index)
{
    TemperaturePrimitiveState result {};
    for (int component = 0; component < fluid_components; ++component) {
        result[static_cast<std::size_t>(component)]
            = field(index.i, index.j, index.k, component);
    }
    return result;
}

void store_temperature(Field<Real>& field, Index3 index,
    const TemperaturePrimitiveState& state)
{
    for (int component = 0; component < fluid_components; ++component) {
        field(index.i, index.j, index.k, component)
            = state[static_cast<std::size_t>(component)];
    }
}

Index3 mirror_index(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    const auto axis = static_cast<std::size_t>(location.axis);
    face[axis] = location.side == Side::Lower ? layer - 1 : extent[axis] - layer;
    return face;
}

Index3 ghost_index(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    const auto axis = static_cast<std::size_t>(location.axis);
    face[axis] = location.side == Side::Lower ? -layer : extent[axis] - 1 + layer;
    return face;
}

TemperaturePrimitiveState reflected_velocity(
    TemperaturePrimitiveState state,
    Normal3 normal,
    const std::array<Real, 3>& wall_velocity,
    bool no_slip)
{
    if (no_slip) {
        state[temperature_velocity_x] = 2.0 * wall_velocity[0]
            - state[temperature_velocity_x];
        state[temperature_velocity_y] = 2.0 * wall_velocity[1]
            - state[temperature_velocity_y];
        state[temperature_velocity_z] = 2.0 * wall_velocity[2]
            - state[temperature_velocity_z];
        return state;
    }
    const Real relative_x = state[temperature_velocity_x] - wall_velocity[0];
    const Real relative_y = state[temperature_velocity_y] - wall_velocity[1];
    const Real relative_z = state[temperature_velocity_z] - wall_velocity[2];
    const Real normal_velocity
        = relative_x * normal.x + relative_y * normal.y + relative_z * normal.z;
    state[temperature_velocity_x] -= 2.0 * normal_velocity * normal.x;
    state[temperature_velocity_y] -= 2.0 * normal_velocity * normal.y;
    state[temperature_velocity_z] -= 2.0 * normal_velocity * normal.z;
    return state;
}

Normal3 patch_outward_normal(
    const StructuredBlock& block, const BoundaryPatch& patch, Index3 face)
{
    const FaceMetric* metric = nullptr;
    switch (patch.face.axis) {
    case Axis::I: metric = &block.face_metrics.i_faces; break;
    case Axis::J: metric = &block.face_metrics.j_faces; break;
    case Axis::K: metric = &block.face_metrics.k_faces; break;
    }
    const Real sign = patch.face.side == Side::Lower ? -1.0 : 1.0;
    Normal3 normal {
        sign * metric->normal_x(face.i, face.j, face.k),
        sign * metric->normal_y(face.i, face.j, face.k),
        sign * metric->normal_z(face.i, face.j, face.k),
    };
    const Real norm = std::sqrt(
        normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!finite(norm) || norm <= 0.0) {
        throw PhysicsConfigurationError("physical boundary has an invalid face normal");
    }
    normal.x /= norm;
    normal.y /= norm;
    normal.z /= norm;
    return normal;
}

TemperaturePrimitiveState wall_ghost(
    const TemperaturePrimitiveState& interior,
    Normal3 normal,
    const BoundaryData& data,
    BoundaryType type,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    const auto interior_pressure
        = pressure_primitive(interior, gas, reference, floors, dimension);
    const bool no_slip = type == BoundaryType::NoSlipAdiabaticWall
        || type == BoundaryType::NoSlipIsothermalWall;
    auto ghost = reflected_velocity(interior, normal, data.wall_velocity, no_slip);
    if (type == BoundaryType::NoSlipIsothermalWall) {
        ghost[temperature_value] = 2.0 * *data.wall_temperature
            - interior[temperature_value];
    } else {
        ghost[temperature_value] = interior[temperature_value];
    }
    ghost[temperature_density] = gas.gamma() * reference.mach() * reference.mach()
        * interior_pressure[4] / ghost[temperature_value];
    if (dimension == 2) {
        ghost[temperature_velocity_z] = 0.0;
    }
    static_cast<void>(pressure_primitive(ghost, gas, reference, floors, dimension));
    return ghost;
}

TemperaturePrimitiveState make_ghost(
    const TemperaturePrimitiveState& interior,
    Normal3 normal,
    const BoundaryData& data,
    BoundaryType type,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    switch (type) {
    case BoundaryType::Farfield:
    case BoundaryType::Inflow:
        return *data.target_state;
    case BoundaryType::Outflow:
        return interior;
    case BoundaryType::SlipWall:
    case BoundaryType::Symmetry:
    case BoundaryType::NoSlipAdiabaticWall:
    case BoundaryType::NoSlipIsothermalWall:
        return wall_ghost(
            interior, normal, data, type, gas, reference, floors, dimension);
    case BoundaryType::Periodic:
        throw PhysicsConfigurationError("periodic boundary must use a connectivity");
    case BoundaryType::Undefined:
        throw PhysicsConfigurationError("undefined physical boundary type");
    }
    throw PhysicsConfigurationError("unsupported physical boundary type");
}

PressurePrimitiveState reflected_face_trace(
    PressurePrimitiveState state,
    Normal3 normal,
    const std::array<Real, 3>& wall_velocity)
{
    const Real relative_x = state[1] - wall_velocity[0];
    const Real relative_y = state[2] - wall_velocity[1];
    const Real relative_z = state[3] - wall_velocity[2];
    const Real normal_velocity
        = relative_x * normal.x + relative_y * normal.y + relative_z * normal.z;
    state[1] -= 2.0 * normal_velocity * normal.x;
    state[2] -= 2.0 * normal_velocity * normal.y;
    state[3] -= 2.0 * normal_velocity * normal.z;
    return state;
}

} // namespace

void BoundaryData::validate(BoundaryType type, int dimension) const
{
    if (dimension != 2 && dimension != 3) {
        throw PhysicsConfigurationError("boundary data dimension must be two or three");
    }
    for (const Real value : wall_velocity) {
        if (!finite(value)) {
            throw PhysicsConfigurationError("wall velocity must be finite");
        }
    }
    if (dimension == 2 && wall_velocity[2] != 0.0) {
        throw PhysicsConfigurationError("two-dimensional wall z velocity must be zero");
    }
    const bool needs_target = type == BoundaryType::Farfield
        || type == BoundaryType::Inflow;
    if (needs_target != target_state.has_value()) {
        throw PhysicsConfigurationError(
            needs_target ? "boundary target state is required"
                         : "boundary target state is not valid for this boundary type");
    }
    const bool needs_temperature = type == BoundaryType::NoSlipIsothermalWall;
    if (needs_temperature != wall_temperature.has_value()) {
        throw PhysicsConfigurationError(
            needs_temperature ? "isothermal wall temperature is required"
                              : "wall temperature is only valid for an isothermal wall");
    }
    if (wall_temperature && (!finite(*wall_temperature) || *wall_temperature <= 0.0)) {
        throw PhysicsConfigurationError("wall temperature must be positive and finite");
    }
}

PhysicalGhostFillResult PhysicalGhostStateOperator::fill(
    StructuredBlock& block,
    const BoundaryDataMap& boundary_data,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    std::uint64_t version)
{
    if (version == 0) {
        throw PhysicsConfigurationError("physical ghost version must be non-zero");
    }
    if (block.ghost_width() < 3) {
        throw PhysicsConfigurationError("WCNS physical ghost fill requires three layers");
    }
    PhysicalGhostFillResult result {version, 0};
    for (const auto& patch : block.boundaries) {
        const auto data_iterator = boundary_data.find(patch.name);
        if (data_iterator == boundary_data.end()) {
            throw PhysicsConfigurationError("physical boundary data is missing for patch " + patch.name);
        }
        const auto& data = data_iterator->second;
        data.validate(patch.type, block.cell_dimension());
        const auto counts = patch.boundary_face_range.counts();
        for (int ok = 0; ok < counts.nk; ++ok) {
            for (int oj = 0; oj < counts.nj; ++oj) {
                for (int oi = 0; oi < counts.ni; ++oi) {
                    const auto face = patch.boundary_face_range.at({oi, oj, ok});
                    const auto normal = patch_outward_normal(block, patch, face);
                    for (int layer = 1; layer <= 3; ++layer) {
                        const auto interior_index = mirror_index(
                            face, patch.face, layer, block.cell_extent());
                        const auto ghost = make_ghost(
                            load_temperature(
                                block.flow.temperature_primitive, interior_index),
                            normal, data, patch.type, gas, reference, floors,
                            block.cell_dimension());
                        const auto destination = ghost_index(
                            face, patch.face, layer, block.cell_extent());
                        const auto pressure = pressure_primitive(
                            ghost, gas, reference, floors, block.cell_dimension());
                        const auto conservative = thermodynamic_conservative(
                            ghost, gas, reference, floors, block.cell_dimension());
                        store_temperature(
                            block.flow.temperature_primitive, destination, ghost);
                        store_state(block.flow.primitive, destination, pressure);
                        store_state(block.flow.conservative, destination, conservative);
                        ++result.state_count;
                    }
                }
            }
        }
    }
    block.flow.physical_ghost_version = version;
    return result;
}

std::string InviscidBoundaryOptions::summary() const
{
    return std::string("strong_boundary_face_state=")
        + (strong_boundary_face_state ? "true" : "false");
}

std::string InviscidBoundaryOptions::restart_signature() const
{
    return "inviscid_boundary_v1;" + summary();
}

PressurePrimitiveState apply_inviscid_boundary_face_state(
    const BoundaryPatch& patch,
    const PressurePrimitiveState& interior_trace,
    const PressurePrimitiveState& reconstructed_exterior_trace,
    Normal3 outward_unit_normal,
    const BoundaryData& data,
    const InviscidBoundaryOptions& options,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    require_unit_normal(outward_unit_normal);
    data.validate(patch.type, dimension);
    static_cast<void>(temperature_primitive(
        interior_trace, gas, reference, floors, dimension));
    static_cast<void>(temperature_primitive(
        reconstructed_exterior_trace, gas, reference, floors, dimension));
    if (!options.strong_boundary_face_state) {
        return reconstructed_exterior_trace;
    }
    switch (patch.type) {
    case BoundaryType::SlipWall:
    case BoundaryType::Symmetry:
    case BoundaryType::NoSlipAdiabaticWall:
    case BoundaryType::NoSlipIsothermalWall:
        return reflected_face_trace(
            interior_trace, outward_unit_normal, data.wall_velocity);
    case BoundaryType::Farfield:
    case BoundaryType::Inflow:
        return pressure_primitive(
            *data.target_state, gas, reference, floors, dimension);
    case BoundaryType::Outflow:
        return interior_trace;
    case BoundaryType::Periodic:
        throw PhysicsConfigurationError("periodic boundary must use a connectivity");
    case BoundaryType::Undefined:
        throw PhysicsConfigurationError("undefined physical boundary type");
    }
    throw PhysicsConfigurationError("unsupported inviscid boundary type");
}

void update_temperature_primitive_interior(
    StructuredBlock& block,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors)
{
    const auto extent = block.cell_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                update_temperature_primitive_cell(
                    block, {i, j, k}, gas, reference, floors);
            }
        }
    }
}

void update_temperature_primitive_cell(
    StructuredBlock& block,
    Index3 index,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors)
{
    const auto conservative = load_conservative(block.flow.conservative, index);
    const auto temperature = temperature_primitive_from_conservative(
        conservative, gas, reference, floors, block.cell_dimension());
    const auto pressure = pressure_primitive(
        temperature, gas, reference, floors, block.cell_dimension());
    store_temperature(block.flow.temperature_primitive, index, temperature);
    store_state(block.flow.primitive, index, pressure);
}

} // namespace wcns
