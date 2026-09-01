#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

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

Real wcns5_left_scaled(
    const std::array<Real, 5>& q,
    Real scale,
    const WcnsParameters& parameters)
{
    parameters.validate();
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("WCNS reconstruction scale must be positive and finite");
    }
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
    constexpr std::array<Real, 3> weights {{1.0 / 16.0, 10.0 / 16.0, 5.0 / 16.0}};
    std::array<Real, 3> alpha {};
    Real sum = 0.0;
    for (std::size_t candidate = 0; candidate < alpha.size(); ++candidate) {
        const Real normalized = beta[candidate] / (scale * scale);
        alpha[candidate] = weights[candidate]
            / std::pow(parameters.epsilon + normalized, parameters.nonlinear_power);
        sum += alpha[candidate];
    }
    if (!std::isfinite(sum) || sum <= 0.0) {
        throw PhysicsError("WCNS nonlinear weights are invalid");
    }
    Real result = 0.0;
    for (std::size_t candidate = 0; candidate < alpha.size(); ++candidate) {
        result += alpha[candidate] * candidates[candidate] / sum;
    }
    return result;
}

std::array<Real, 6> component_stencil(
    const Field<Real>& field, Axis axis, Index3 face, int component)
{
    std::array<Real, 6> stencil {};
    const auto left = shifted(face, axis, -1);
    for (int point = 0; point < 6; ++point) {
        const auto index = shifted(left, axis, point - 2);
        stencil[static_cast<std::size_t>(point)]
            = field(index.i, index.j, index.k, component);
    }
    return stencil;
}

PressurePrimitiveState convert_reconstructed(
    const std::array<Real, euler_components>& state,
    ReconstructionVariables variables,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    if (variables == ReconstructionVariables::Primitive) {
        static_cast<void>(temperature_primitive(
            state, gas, reference, floors, dimension));
        return state;
    }
    return pressure_primitive(
        temperature_primitive_from_conservative(
            state, gas, reference, floors, dimension),
        gas, reference, floors, dimension);
}

bool try_convert(
    const std::array<Real, euler_components>& left,
    const std::array<Real, euler_components>& right,
    const ReconstructionConfig& config,
    const GasModel& gas,
    const ReferenceScales& reference,
    int dimension,
    EulerFaceStates& output)
{
    try {
        output.left = convert_reconstructed(
            left, config.variables, gas, reference, config.floors, dimension);
        output.right = convert_reconstructed(
            right, config.variables, gas, reference, config.floors, dimension);
        return true;
    } catch (const PhysicsConfigurationError&) {
        return false;
    }
}

} // namespace

void WcnsParameters::validate() const
{
    if (!std::isfinite(epsilon) || epsilon <= 0.0 || nonlinear_power <= 0) {
        throw std::invalid_argument("WCNS parameters require positive epsilon and power");
    }
}

void ReconstructionConfig::validate() const
{
    nonlinear.validate();
    scaling.validate();
    floors.validate();
    if (nonlinear.epsilon != floors.reconstruction_epsilon
        || scaling.epsilon != floors.reconstruction_epsilon
        || scaling.scale_floor != floors.reconstruction_scale) {
        throw std::invalid_argument(
            "reconstruction epsilon and scale floor must come from NumericalFloors");
    }
    switch (kind) {
    case ReconstructionKind::Linear5:
    case ReconstructionKind::WcnsJs: break;
    default: throw std::invalid_argument("unknown reconstruction kind");
    }
    switch (variables) {
    case ReconstructionVariables::Conservative:
    case ReconstructionVariables::Primitive: break;
    default: throw std::invalid_argument("unknown reconstruction variable set");
    }
}

std::string ReconstructionConfig::summary() const
{
    validate();
    const char* kind_name = kind == ReconstructionKind::Linear5
        ? "linear5" : "wcns_js";
    const char* variable_name = variables == ReconstructionVariables::Conservative
        ? "conservative" : "primitive";
    return std::string("reconstruction=") + kind_name
        + ";variables=" + variable_name;
}

std::string ReconstructionConfig::restart_signature() const
{
    return "reconstruction_v1;" + summary();
}

Real wcns5_left_interpolation(
    const std::array<Real, 5>& q,
    const WcnsParameters& parameters)
{
    return wcns5_left_scaled(q, 1.0, parameters);
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

ScalarFaceStates linear5_reconstruct(const std::array<Real, 6>& q)
{
    for (const auto value : q) {
        if (!std::isfinite(value)) {
            throw PhysicsError("linear reconstruction stencil contains a non-finite value");
        }
    }
    return {
        (3.0 * q[0] - 20.0 * q[1] + 90.0 * q[2]
            + 60.0 * q[3] - 5.0 * q[4]) / 128.0,
        (-5.0 * q[1] + 60.0 * q[2] + 90.0 * q[3]
            - 20.0 * q[4] + 3.0 * q[5]) / 128.0,
    };
}

ScalarFaceStates wcns5_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters)
{
    const std::array<Real, 5> left {{
        stencil[0], stencil[1], stencil[2], stencil[3], stencil[4]}};
    const std::array<Real, 5> right {{
        stencil[5], stencil[4], stencil[3], stencil[2], stencil[1]}};
    return {
        wcns5_left_scaled(left, scale, parameters),
        wcns5_left_scaled(right, scale, parameters),
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

EulerFaceStates reconstruct_thermodynamic_face(
    const Field<Real>& conservative,
    const Field<Real>& pressure_primitive_field,
    Axis axis,
    Index3 face,
    const ReconstructionConfig& config,
    const GasModel& gas,
    const ReferenceScales& reference,
    ReconstructionDiagnostics& diagnostics,
    int dimension)
{
    config.validate();
    const auto& source = config.variables == ReconstructionVariables::Conservative
        ? conservative : pressure_primitive_field;
    if (source.components() != euler_components || source.ghost_width() < 3) {
        throw std::invalid_argument("thermodynamic reconstruction requires five components and three ghost layers");
    }
    std::array<Real, euler_components> left {};
    std::array<Real, euler_components> right {};
    const auto reconstruct = [&](bool nonlinear) {
        for (int component = 0; component < euler_components; ++component) {
            const auto stencil = component_stencil(source, axis, face, component);
            const auto states = nonlinear
                ? wcns5_reconstruct_scaled(
                    stencil,
                    std::max(
                        config.scaling.component[static_cast<std::size_t>(component)],
                        config.scaling.scale_floor),
                    config.nonlinear)
                : linear5_reconstruct(stencil);
            left[static_cast<std::size_t>(component)] = states.left;
            right[static_cast<std::size_t>(component)] = states.right;
        }
    };

    EulerFaceStates result;
    if (config.kind == ReconstructionKind::WcnsJs) {
        reconstruct(true);
        ++diagnostics.nonlinear_faces;
        if (try_convert(left, right, config, gas, reference, dimension, result)) {
            return result;
        }
        ++diagnostics.linear_fallbacks;
    }
    reconstruct(false);
    ++diagnostics.linear_faces;
    if (try_convert(left, right, config, gas, reference, dimension, result)) {
        return result;
    }

    const auto left_index = shifted(face, axis, -1);
    const auto right_index = face;
    const auto load_first_order = [&](Index3 index) {
        const auto state = config.variables == ReconstructionVariables::Conservative
            ? load_conservative(source, index) : load_primitive(source, index);
        return convert_reconstructed(
            state, config.variables, gas, reference, config.floors, dimension);
    };
    result.left = load_first_order(left_index);
    result.right = load_first_order(right_index);
    ++diagnostics.first_order_fallbacks;
    return result;
}

} // namespace wcns
