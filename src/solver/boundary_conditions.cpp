#include <wcns/solver/boundary_conditions.hpp>

#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

const FaceMetric& face_metric(const StructuredBlock& block, Axis axis)
{
    switch (axis) {
    case Axis::I:
        return block.face_metrics.i_faces;
    case Axis::J:
        return block.face_metrics.j_faces;
    case Axis::K:
        return block.face_metrics.k_faces;
    }
    throw std::invalid_argument("invalid boundary axis");
}

Normal3 outward_normal(
    const StructuredBlock& block,
    const BoundaryPatch& patch,
    Index3 face)
{
    const auto& metric = face_metric(block, patch.face.axis);
    const Real sign = patch.face.side == Side::Lower ? -1.0 : 1.0;
    Normal3 normal {
        sign * metric.normal_x(face.i, face.j, face.k),
        sign * metric.normal_y(face.i, face.j, face.k),
        sign * metric.normal_z(face.i, face.j, face.k),
    };
    const Real magnitude = std::sqrt(
        normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(magnitude) || magnitude <= 0.0) {
        throw PhysicsError("boundary face has an invalid normal");
    }
    normal.x /= magnitude;
    normal.y /= magnitude;
    normal.z /= magnitude;
    return normal;
}

PrimitiveState reflected(PrimitiveState state, Normal3 normal)
{
    const Real normal_velocity = state[velocity_x] * normal.x
        + state[velocity_y] * normal.y + state[velocity_z] * normal.z;
    state[velocity_x] -= 2.0 * normal_velocity * normal.x;
    state[velocity_y] -= 2.0 * normal_velocity * normal.y;
    state[velocity_z] -= 2.0 * normal_velocity * normal.z;
    return state;
}

void store_primitive_and_conservative(
    StructuredBlock& block,
    Index3 index,
    const PrimitiveState& primitive,
    const IdealGas& gas)
{
    store_state(block.flow.primitive, index, primitive);
    store_state(block.flow.conservative, index, to_conservative(primitive, gas));
}

Index3 ghost_index(
    Index3 face,
    Axis axis,
    Side side,
    int layer,
    const Extent3& extent)
{
    const auto a = static_cast<std::size_t>(axis);
    face[a] = side == Side::Lower ? -layer : extent[a] - 1 + layer;
    return face;
}

Index3 adjacent_index(Index3 face, Axis axis, Side side, const Extent3& extent)
{
    const auto a = static_cast<std::size_t>(axis);
    face[a] = side == Side::Lower ? 0 : extent[a] - 1;
    return face;
}

Index3 mirror_index(
    Index3 face,
    Axis axis,
    Side side,
    int layer,
    const Extent3& extent)
{
    const auto a = static_cast<std::size_t>(axis);
    face[a] = side == Side::Lower ? layer - 1 : extent[a] - layer;
    return face;
}

PrimitiveState boundary_state(
    StructuredBlock& block,
    const BoundaryPatch& patch,
    Index3 face,
    int layer,
    const PrimitiveState& prescribed)
{
    const auto extent = block.cell_extent();
    switch (patch.type) {
    case BoundaryType::Farfield:
    case BoundaryType::Inflow:
        return prescribed;
    case BoundaryType::Outflow:
        return load_primitive(
            block.flow.primitive,
            adjacent_index(face, patch.face.axis, patch.face.side, extent));
    case BoundaryType::SlipWall:
    case BoundaryType::Symmetry:
        return reflected(
            load_primitive(
                block.flow.primitive,
                mirror_index(face, patch.face.axis, patch.face.side, layer, extent)),
            outward_normal(block, patch, face));
    case BoundaryType::NoSlipAdiabaticWall:
    case BoundaryType::NoSlipIsothermalWall:
        throw PhysicsError("no-slip boundary conditions require a viscous solver");
    case BoundaryType::Periodic:
        throw PhysicsError("periodic boundaries must be represented as connectivities");
    case BoundaryType::Undefined:
        throw PhysicsError("undefined physical boundary type");
    }
    throw PhysicsError("unsupported physical boundary type");
}

} // namespace

void update_primitive_cell(
    StructuredBlock& block,
    Index3 index,
    const IdealGas& gas)
{
    store_state(
        block.flow.primitive,
        index,
        to_primitive(load_conservative(block.flow.conservative, index), gas));
}

void update_primitive_interior(StructuredBlock& block, const IdealGas& gas)
{
    const auto extent = block.cell_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                update_primitive_cell(block, {i, j, k}, gas);
            }
        }
    }
}

void fill_physical_boundaries(
    StructuredBlock& block,
    const PrimitiveState& prescribed,
    const IdealGas& gas)
{
    static_cast<void>(to_conservative(prescribed, gas));
    const int ghost_width = block.ghost_width();
    for (const auto& patch : block.boundaries) {
        if (patch.face.axis == Axis::K && block.cell_dimension() == 2) {
            throw PhysicsError("a 2D block cannot have a K-normal physical boundary");
        }
        const auto counts = patch.cell_face_range.counts();
        for (int ok = 0; ok < counts.nk; ++ok) {
            for (int oj = 0; oj < counts.nj; ++oj) {
                for (int oi = 0; oi < counts.ni; ++oi) {
                    const auto face = patch.cell_face_range.at({oi, oj, ok});
                    for (int layer = 1; layer <= ghost_width; ++layer) {
                        const auto ghost = ghost_index(
                            face,
                            patch.face.axis,
                            patch.face.side,
                            layer,
                            block.cell_extent());
                        store_primitive_and_conservative(
                            block,
                            ghost,
                            boundary_state(block, patch, face, layer, prescribed),
                            gas);
                    }
                }
            }
        }
    }
}

} // namespace wcns
