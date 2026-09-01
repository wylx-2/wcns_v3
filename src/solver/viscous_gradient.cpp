#include <wcns/solver/viscous_gradient.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace wcns {
namespace {

int operand_component(ViscousPrimitive variable, int direction)
{
    if (direction < 0 || direction >= 3) {
        throw std::out_of_range("Cartesian gradient direction is outside [0,3)");
    }
    return static_cast<int>(variable) * 3 + direction;
}

bool contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) return false;
    }
    return true;
}

bool connection_covers(
    const StructuredBlock& block, Axis axis, Side side, Index3 face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis
            && connection.receiver_face.side == side
            && contains(connection.shared_face_range.untyped(), face)) {
            return true;
        }
    }
    return false;
}

const FaceAreaVectors& face_metrics(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(
    const MetricField& metric, Axis axis, Index3 face)
{
    const auto& values = face_metrics(metric, axis);
    return {{
        values.x(face.i, face.j, face.k),
        values.y(face.i, face.j, face.k),
        values.z(face.i, face.j, face.k),
    }};
}

template<std::size_t N>
Real interpolate_component(
    const Field<Real>& field,
    Axis axis,
    Index3 face,
    int component,
    const std::array<int, N>& offsets,
    const std::array<Real, N>& coefficients)
{
    Real result = 0.0;
    for (std::size_t term = 0; term < N; ++term) {
        auto center = face;
        center[static_cast<std::size_t>(axis)] += offsets[term];
        result += coefficients[term]
            * field(center.i, center.j, center.k, component);
    }
    if (!std::isfinite(result)) {
        throw PhysicsError("viscous face interpolation produced a non-finite value");
    }
    return result;
}

Real centered_derivative(
    const Field<Real>& field,
    Axis axis,
    Index3 cell,
    int component,
    AlgorithmProfileKind profile)
{
    const auto value = [&](int offset) {
        auto face = cell;
        face[static_cast<std::size_t>(axis)] += offset;
        const Real result = field(face.i, face.j, face.k, component);
        if (!std::isfinite(result)) {
            throw PhysicsError("gradient operand face halo is incomplete");
        }
        return result;
    };
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return (value(-1) - 27.0 * value(0) + 27.0 * value(1) - value(2))
            / 24.0;
    }
    return (-9.0 * value(-2) + 125.0 * value(-1) - 2250.0 * value(0)
        + 2250.0 * value(1) - 125.0 * value(2) + 9.0 * value(3))
        / 1920.0;
}

} // namespace

GradientOperandFaceField::GradientOperandFaceField(
    Extent3 cells,
    int dimension,
    AlgorithmProfileKind profile,
    std::uint64_t version)
    : profile_(profile)
    , version_(version)
    , halo_layers_(profile == AlgorithmProfileKind::PhengleiWcns ? 1 : 2)
    , dimension_(dimension)
    , i_({cells.ni + 1, cells.nj, cells.nk}, gradient_operand_components,
          halo_layers_, std::numeric_limits<Real>::quiet_NaN())
    , j_({cells.ni, cells.nj + 1, cells.nk}, gradient_operand_components,
          halo_layers_, std::numeric_limits<Real>::quiet_NaN())
    , k_({cells.ni, cells.nj, cells.nk + 1}, gradient_operand_components,
          halo_layers_, std::numeric_limits<Real>::quiet_NaN())
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("gradient operand dimension must be two or three");
    }
}

Field<Real>& GradientOperandFaceField::field(Axis axis)
{
    if (axis == Axis::I) return i_;
    if (axis == Axis::J) return j_;
    if (dimension_ != 3) throw std::out_of_range("two-dimensional gradient has no K faces");
    return k_;
}

const Field<Real>& GradientOperandFaceField::field(Axis axis) const
{
    return const_cast<GradientOperandFaceField*>(this)->field(axis);
}

PrimitiveGradientField::PrimitiveGradientField(
    Extent3 cells,
    int dimension,
    AlgorithmProfileKind profile,
    std::uint64_t version)
    : profile_(profile)
    , version_(version)
    , halo_layers_(profile == AlgorithmProfileKind::PhengleiWcns ? 2 : 3)
    , dimension_(dimension)
    , values_(cells, gradient_operand_components, halo_layers_,
          std::numeric_limits<Real>::quiet_NaN())
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("primitive gradient dimension must be two or three");
    }
}

Real& PrimitiveGradientField::operator()(
    Index3 cell, ViscousPrimitive variable, int direction)
{
    return values_(cell.i, cell.j, cell.k, operand_component(variable, direction));
}

Real PrimitiveGradientField::operator()(
    Index3 cell, ViscousPrimitive variable, int direction) const
{
    return values_(cell.i, cell.j, cell.k, operand_component(variable, direction));
}

TemperaturePrimitiveState interpolate_temperature_face(
    const StructuredBlock& block,
    const AlgorithmProfile& profile,
    Axis axis,
    Index3 face)
{
    const auto& field = block.flow.temperature_primitive;
    TemperaturePrimitiveState result {};
    if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
        constexpr std::array<int, 4> offsets {{-2, -1, 0, 1}};
        constexpr std::array<Real, 4> coefficients {{
            -1.0 / 16.0, 9.0 / 16.0, 9.0 / 16.0, -1.0 / 16.0}};
        for (int component = 0; component < fluid_components; ++component) {
            result[static_cast<std::size_t>(component)] = interpolate_component(
                field, axis, face, component, offsets, coefficients);
        }
    } else {
        constexpr std::array<int, 6> offsets {{-3, -2, -1, 0, 1, 2}};
        constexpr std::array<Real, 6> coefficients {{
            3.0 / 256.0, -25.0 / 256.0, 150.0 / 256.0,
            150.0 / 256.0, -25.0 / 256.0, 3.0 / 256.0}};
        for (int component = 0; component < fluid_components; ++component) {
            result[static_cast<std::size_t>(component)] = interpolate_component(
                field, axis, face, component, offsets, coefficients);
        }
    }
    if (block.cell_dimension() == 2) result[temperature_velocity_z] = 0.0;
    return result;
}

GradientOperandFaceField compute_gradient_face_operands(
    const StructuredBlock& block,
    const MetricField& metric,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    if (metric.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()) {
        throw ProfileError("gradient operands use incompatible metric/profile metadata");
    }
    GradientOperandFaceField result(
        block.cell_extent(), block.cell_dimension(), profile.kind(), version);
    const auto cells = block.cell_extent();
    const auto compute_axis = [&](Axis axis) {
        auto faces = cells;
        ++faces[static_cast<std::size_t>(axis)];
        auto& output = result.field(axis);
        for (int k = 0; k < faces.nk; ++k) {
            for (int j = 0; j < faces.nj; ++j) {
                for (int i = 0; i < faces.ni; ++i) {
                    const Index3 face {i, j, k};
                    const auto state = interpolate_temperature_face(
                        block, profile, axis, face);
                    const auto area = area_vector(metric, axis, face);
                    const std::array<Real, viscous_primitive_components> q {{
                        state[temperature_velocity_x],
                        state[temperature_velocity_y],
                        state[temperature_velocity_z],
                        state[temperature_value],
                    }};
                    for (int variable = 0; variable < viscous_primitive_components;
                         ++variable) {
                        for (int direction = 0; direction < 3; ++direction) {
                            output(i, j, k, variable * 3 + direction)
                                = q[static_cast<std::size_t>(variable)]
                                * area[static_cast<std::size_t>(direction)];
                        }
                    }
                }
            }
        }
    };
    compute_axis(Axis::I);
    compute_axis(Axis::J);
    if (block.cell_dimension() == 3) compute_axis(Axis::K);
    return result;
}

PrimitiveGradientField compute_primitive_gradients(
    const StructuredBlock& block,
    const MetricField& metric,
    const GradientOperandFaceField& operands,
    const AlgorithmProfile& profile)
{
    if (metric.profile() != profile.kind() || operands.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()
        || operands.dimension() != block.cell_dimension()) {
        throw ProfileError("primitive gradient inputs use incompatible metadata");
    }
    PrimitiveGradientField result(
        block.cell_extent(), block.cell_dimension(), profile.kind(), operands.version());
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Index3 cell {i, j, k};
                const Real jacobian = metric.jacobian()(i, j, k);
                if (!std::isfinite(jacobian) || jacobian <= 0.0) {
                    throw PhysicsError("primitive gradient has an invalid Jacobian");
                }
                for (int variable = 0; variable < viscous_primitive_components;
                     ++variable) {
                    for (int direction = 0; direction < 3; ++direction) {
                        Real divergence = 0.0;
                        for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                            const auto axis = static_cast<Axis>(logical);
                            const int count = cells[static_cast<std::size_t>(axis)];
                            const int normal = cell[static_cast<std::size_t>(axis)];
                            Index3 lower = cell;
                            lower[static_cast<std::size_t>(axis)] = 0;
                            Index3 upper = cell;
                            upper[static_cast<std::size_t>(axis)] = count;
                            const bool connected_lower
                                = connection_covers(block, axis, Side::Lower, lower);
                            const bool connected_upper
                                = connection_covers(block, axis, Side::Upper, upper);
                            const int width = profile.kind()
                                    == AlgorithmProfileKind::PhengleiWcns
                                ? 1 : 2;
                            const auto& values = operands.field(axis);
                            const int component = variable * 3 + direction;
                            if ((connected_lower && normal < width)
                                || (connected_upper && normal >= count - width)) {
                                divergence += centered_derivative(
                                    values, axis, cell, component, profile.kind());
                            } else {
                                const auto operators = LineOperators::build(profile, count);
                                const auto& row = operators.derivative_rows()[
                                    static_cast<std::size_t>(normal)];
                                for (const auto [face_index, coefficient] : row) {
                                    auto face = cell;
                                    face[static_cast<std::size_t>(axis)] = face_index;
                                    divergence += coefficient * values(
                                        face.i, face.j, face.k, component);
                                }
                            }
                        }
                        result(cell, static_cast<ViscousPrimitive>(variable), direction)
                            = divergence / jacobian;
                    }
                }
                if (block.cell_dimension() == 2) {
                    for (int variable = 0; variable < viscous_primitive_components;
                         ++variable) {
                        result(cell, static_cast<ViscousPrimitive>(variable), 2) = 0.0;
                    }
                    for (int direction = 0; direction < 3; ++direction) {
                        result(cell, ViscousPrimitive::VelocityZ, direction) = 0.0;
                    }
                }
            }
        }
    }
    return result;
}

} // namespace wcns
