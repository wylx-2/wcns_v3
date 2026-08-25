#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

Real square(Real value)
{
    return value * value;
}

Index3 shifted(Index3 index, Axis axis, int offset)
{
    index[static_cast<std::size_t>(axis)] += offset;
    return index;
}

} // namespace

void WcnsParameters::validate() const
{
    if (!std::isfinite(epsilon) || epsilon <= 0.0 || nonlinear_power <= 0) {
        throw std::invalid_argument("WCNS parameters require positive epsilon and power");
    }
}

Real wcns5_left_interpolation(
    const std::array<Real, 5>& q,
    const WcnsParameters& parameters)
{
    parameters.validate();
    for (const auto value : q) {
        if (!std::isfinite(value)) {
            throw PhysicsError("WCNS stencil contains a non-finite value");
        }
    }

    const std::array<Real, 3> candidates {{
        (3.0 * q[0] - 10.0 * q[1] + 15.0 * q[2]) / 8.0,
        (-q[1] + 6.0 * q[2] + 3.0 * q[3]) / 8.0,
        (3.0 * q[2] + 6.0 * q[3] - q[4]) / 8.0,
    }};
    const std::array<Real, 3> beta {{
        (13.0 / 12.0) * square(q[0] - 2.0 * q[1] + q[2])
            + 0.25 * square(q[0] - 4.0 * q[1] + 3.0 * q[2]),
        (13.0 / 12.0) * square(q[1] - 2.0 * q[2] + q[3])
            + 0.25 * square(q[1] - q[3]),
        (13.0 / 12.0) * square(q[2] - 2.0 * q[3] + q[4])
            + 0.25 * square(3.0 * q[2] - 4.0 * q[3] + q[4]),
    }};
    constexpr std::array<Real, 3> linear_weights {{1.0 / 16.0, 10.0 / 16.0, 5.0 / 16.0}};
    std::array<Real, 3> alpha {};
    Real alpha_sum = 0.0;
    for (std::size_t candidate = 0; candidate < alpha.size(); ++candidate) {
        alpha[candidate] = linear_weights[candidate]
            / std::pow(parameters.epsilon + beta[candidate], parameters.nonlinear_power);
        alpha_sum += alpha[candidate];
    }
    if (!std::isfinite(alpha_sum) || alpha_sum <= 0.0) {
        throw PhysicsError("WCNS nonlinear weights are invalid");
    }
    Real result = 0.0;
    for (std::size_t candidate = 0; candidate < alpha.size(); ++candidate) {
        result += (alpha[candidate] / alpha_sum) * candidates[candidate];
    }
    return result;
}

ScalarFaceStates wcns5_reconstruct(
    const std::array<Real, 6>& stencil,
    const WcnsParameters& parameters)
{
    const std::array<Real, 5> left {{
        stencil[0], stencil[1], stencil[2], stencil[3], stencil[4]}};
    const std::array<Real, 5> right {{
        stencil[5], stencil[4], stencil[3], stencil[2], stencil[1]}};
    return {
        wcns5_left_interpolation(left, parameters),
        wcns5_left_interpolation(right, parameters),
    };
}

EulerFaceStates reconstruct_euler_face(
    const Field<Real>& primitive,
    Axis axis,
    Index3 face,
    const WcnsParameters& parameters)
{
    if (primitive.components() != euler_components || primitive.ghost_width() < 3) {
        throw std::invalid_argument("WCNS Euler reconstruction requires five components and three ghost layers");
    }
    EulerFaceStates result;
    const auto left_center = shifted(face, axis, -1);
    for (int component = 0; component < euler_components; ++component) {
        std::array<Real, 6> stencil {};
        for (int point = 0; point < 6; ++point) {
            const auto index = shifted(left_center, axis, point - 2);
            stencil[static_cast<std::size_t>(point)]
                = primitive(index.i, index.j, index.k, component);
        }
        const auto states = wcns5_reconstruct(stencil, parameters);
        result.left[static_cast<std::size_t>(component)] = states.left;
        result.right[static_cast<std::size_t>(component)] = states.right;
    }
    return result;
}

} // namespace wcns

