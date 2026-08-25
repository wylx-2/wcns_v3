#include <wcns/solver/spatial_operator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
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
    throw std::invalid_argument("invalid flux axis");
}

Index3 shifted(Index3 index, Axis axis, int offset)
{
    index[static_cast<std::size_t>(axis)] += offset;
    return index;
}

bool interior(Index3 index, const Extent3& extent)
{
    return index.i >= 0 && index.i < extent.ni && index.j >= 0
        && index.j < extent.nj && index.k >= 0 && index.k < extent.nk;
}

bool valid_state(const PrimitiveState& state, const IdealGas& gas)
{
    try {
        static_cast<void>(to_conservative(state, gas));
        return true;
    } catch (const PhysicsError&) {
        return false;
    }
}

EulerFaceStates limited_face_states(
    const StructuredBlock& block,
    Axis axis,
    Index3 face,
    const IdealGas& gas,
    const WcnsParameters& parameters)
{
    auto states = reconstruct_euler_face(block.flow.primitive, axis, face, parameters);
    if (!valid_state(states.left, gas) || !valid_state(states.right, gas)) {
        states.left = load_primitive(block.flow.primitive, shifted(face, axis, -1));
        states.right = load_primitive(block.flow.primitive, face);
    }
    if (!valid_state(states.left, gas) || !valid_state(states.right, gas)) {
        throw PhysicsError("neither WCNS nor first-order face states are physical");
    }
    return states;
}

void accumulate_face(
    StructuredBlock& block,
    Axis axis,
    Index3 face,
    const IdealGas& gas,
    const WcnsParameters& parameters)
{
    const auto& metric = face_metric(block, axis);
    const Normal3 normal {
        metric.normal_x(face.i, face.j, face.k),
        metric.normal_y(face.i, face.j, face.k),
        metric.normal_z(face.i, face.j, face.k),
    };
    const Real area = metric.area(face.i, face.j, face.k);
    if (!std::isfinite(area) || area <= 0.0) {
        throw PhysicsError("flux face has a non-positive area");
    }
    const auto states = limited_face_states(block, axis, face, gas, parameters);
    const auto flux = rusanov_flux(states.left, states.right, normal, gas);
    const auto left = shifted(face, axis, -1);
    const auto right = face;
    const auto extent = block.cell_extent();

    if (interior(left, extent)) {
        const Real inverse_volume = 1.0 / block.cell_metrics.volume(left.i, left.j, left.k);
        for (int component = 0; component < euler_components; ++component) {
            block.flow.residual(left.i, left.j, left.k, component)
                -= area * inverse_volume * flux[static_cast<std::size_t>(component)];
        }
    }
    if (interior(right, extent)) {
        const Real inverse_volume = 1.0 / block.cell_metrics.volume(right.i, right.j, right.k);
        for (int component = 0; component < euler_components; ++component) {
            block.flow.residual(right.i, right.j, right.k, component)
                += area * inverse_volume * flux[static_cast<std::size_t>(component)];
        }
    }
}

Real face_spectral_radius(
    const StructuredBlock& block,
    Axis axis,
    Index3 face,
    const PrimitiveState& primitive,
    const IdealGas& gas)
{
    const auto& metric = face_metric(block, axis);
    const Real nx = metric.normal_x(face.i, face.j, face.k);
    const Real ny = metric.normal_y(face.i, face.j, face.k);
    const Real nz = metric.normal_z(face.i, face.j, face.k);
    const Real normal_velocity = primitive[velocity_x] * nx
        + primitive[velocity_y] * ny + primitive[velocity_z] * nz;
    return metric.area(face.i, face.j, face.k)
        * (std::abs(normal_velocity) + sound_speed(primitive, gas));
}

} // namespace

void SpatialParameters::validate() const
{
    gas.validate();
    wcns.validate();
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw std::invalid_argument("CFL number must be positive and finite");
    }
}

void compute_euler_residual(
    StructuredBlock& block,
    const IdealGas& gas,
    const WcnsParameters& parameters)
{
    gas.validate();
    parameters.validate();
    if (block.ghost_width() < 3) {
        throw PhysicsError("WCNS residual requires at least three ghost layers");
    }
    block.flow.residual.fill(0.0);
    const auto extent = block.cell_extent();

    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i <= extent.ni; ++i) {
                accumulate_face(block, Axis::I, {i, j, k}, gas, parameters);
            }
        }
    }
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j <= extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                accumulate_face(block, Axis::J, {i, j, k}, gas, parameters);
            }
        }
    }
    if (block.cell_dimension() == 3) {
        for (int k = 0; k <= extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    accumulate_face(block, Axis::K, {i, j, k}, gas, parameters);
                }
            }
        }
    }
}

Real stable_time_step(const StructuredBlock& block, Real cfl, const IdealGas& gas)
{
    gas.validate();
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw std::invalid_argument("CFL number must be positive and finite");
    }
    const auto extent = block.cell_extent();
    Real result = std::numeric_limits<Real>::infinity();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const auto primitive = load_primitive(block.flow.primitive, cell);
                Real denominator = 0.0;
                denominator += face_spectral_radius(block, Axis::I, {i, j, k}, primitive, gas);
                denominator += face_spectral_radius(block, Axis::I, {i + 1, j, k}, primitive, gas);
                denominator += face_spectral_radius(block, Axis::J, {i, j, k}, primitive, gas);
                denominator += face_spectral_radius(block, Axis::J, {i, j + 1, k}, primitive, gas);
                if (block.cell_dimension() == 3) {
                    denominator += face_spectral_radius(block, Axis::K, {i, j, k}, primitive, gas);
                    denominator += face_spectral_radius(block, Axis::K, {i, j, k + 1}, primitive, gas);
                }
                const Real volume = block.cell_metrics.volume(i, j, k);
                if (!std::isfinite(volume) || volume <= 0.0 || !std::isfinite(denominator)
                    || denominator <= 0.0) {
                    throw PhysicsError("cannot compute a positive stable time step");
                }
                result = std::min(result, cfl * volume / denominator);
            }
        }
    }
    return result;
}

Real residual_l2(const StructuredBlock& block)
{
    const auto extent = block.cell_extent();
    Real sum = 0.0;
    std::size_t count = 0;
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < euler_components; ++component) {
                    const Real value = block.flow.residual(i, j, k, component);
                    sum += value * value;
                    ++count;
                }
            }
        }
    }
    return std::sqrt(sum / static_cast<Real>(count));
}

} // namespace wcns
