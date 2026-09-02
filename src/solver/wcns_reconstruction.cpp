#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace wcns {
namespace {

bool valid_algorithm_name(std::string_view name)
{
    if (name.empty()) return false;
    for (const char character : name) {
        const bool lower = character >= 'a' && character <= 'z';
        const bool digit = character >= '0' && character <= '9';
        if (!lower && !digit && character != '_') return false;
    }
    return true;
}

std::array<Real, 6> checked_stencil(ScalarStencilView stencil)
{
    if (stencil.size() != 6) {
        throw std::invalid_argument("Stage-L reconstruction requires exactly six stencil values");
    }
    std::array<Real, 6> result {};
    std::copy(stencil.begin(), stencil.end(), result.begin());
    for (const auto value : result) {
        if (!std::isfinite(value)) {
            throw PhysicsError("reconstruction stencil contains a non-finite value");
        }
    }
    return result;
}

std::array<Real, 6> orient_stencil(
    const std::array<Real, 6>& stencil, TraceSide side);
Real weno_z_left_scaled(
    const std::array<Real, 6>& q,
    Real scale,
    const WcnsParameters& parameters);
Real mdcd_linear_left(
    const std::array<Real, 6>& q,
    const WcnsParameters& parameters);
Real mdcd_hybrid_left_scaled(
    const std::array<Real, 6>& q,
    Real scale,
    const WcnsParameters& parameters);

class Linear5Scheme final : public IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override
    {
        return "linear5";
    }

    [[nodiscard]] StencilRequirement stencil_requirement() const noexcept override
    {
        return {};
    }

    [[nodiscard]] Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const override
    {
        if (!std::isfinite(context.scale) || context.scale <= 0.0) {
            throw std::invalid_argument("reconstruction scale must be positive and finite");
        }
        const auto values = checked_stencil(stencil);
        const auto result = linear5_reconstruct(values);
        return side == TraceSide::Left ? result.left : result.right;
    }
};

class WenoJsScheme final : public IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override
    {
        return "weno_js";
    }

    [[nodiscard]] StencilRequirement stencil_requirement() const noexcept override
    {
        return {};
    }

    [[nodiscard]] Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const override
    {
        const auto values = checked_stencil(stencil);
        const auto result = wcns5_reconstruct_scaled(
            values, context.scale, context.parameters);
        return side == TraceSide::Left ? result.left : result.right;
    }
};

class WenoZScheme final : public IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "weno_z"; }
    [[nodiscard]] StencilRequirement stencil_requirement() const noexcept override { return {}; }
    [[nodiscard]] Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const override
    {
        const auto values = orient_stencil(checked_stencil(stencil), side);
        return weno_z_left_scaled(values, context.scale, context.parameters);
    }
};

class MdcdLinearScheme final : public IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "mdcd_linear"; }
    [[nodiscard]] StencilRequirement stencil_requirement() const noexcept override { return {}; }
    [[nodiscard]] Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const override
    {
        if (!std::isfinite(context.scale) || context.scale <= 0.0) {
            throw std::invalid_argument("reconstruction scale must be positive and finite");
        }
        const auto values = orient_stencil(checked_stencil(stencil), side);
        return mdcd_linear_left(values, context.parameters);
    }
};

class MdcdHybridScheme final : public IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "mdcd_hybrid"; }
    [[nodiscard]] StencilRequirement stencil_requirement() const noexcept override { return {}; }
    [[nodiscard]] Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const override
    {
        const auto values = orient_stencil(checked_stencil(stencil), side);
        return mdcd_hybrid_left_scaled(values, context.scale, context.parameters);
    }
};

Real square(Real value)
{
    return value * value;
}

std::array<Real, 6> orient_stencil(
    const std::array<Real, 6>& stencil, TraceSide side)
{
    if (side == TraceSide::Left) return stencil;
    auto result = stencil;
    std::reverse(result.begin(), result.end());
    return result;
}

std::array<Real, 3> weno_candidates(const std::array<Real, 6>& q)
{
    return {{
        (3.0 * q[0] - 10.0 * q[1] + 15.0 * q[2]) / 8.0,
        (-q[1] + 6.0 * q[2] + 3.0 * q[3]) / 8.0,
        (3.0 * q[2] + 6.0 * q[3] - q[4]) / 8.0,
    }};
}

std::array<Real, 3> weno_smoothness(const std::array<Real, 6>& q, Real scale)
{
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("reconstruction scale must be positive and finite");
    }
    std::array<Real, 6> normalized {};
    for (std::size_t index = 0; index < q.size(); ++index) {
        normalized[index] = q[index] / scale;
    }
    return {{
        (13.0 / 12.0) * square(normalized[0] - 2.0 * normalized[1] + normalized[2])
            + 0.25 * square(normalized[0] - 4.0 * normalized[1] + 3.0 * normalized[2]),
        (13.0 / 12.0) * square(normalized[1] - 2.0 * normalized[2] + normalized[3])
            + 0.25 * square(normalized[1] - normalized[3]),
        (13.0 / 12.0) * square(normalized[2] - 2.0 * normalized[3] + normalized[4])
            + 0.25 * square(3.0 * normalized[2] - 4.0 * normalized[3] + normalized[4]),
    }};
}

Real checked_weighted_sum(
    const std::array<Real, 4>& alpha,
    const std::array<Real, 4>& candidates,
    std::size_t count,
    const char* label)
{
    Real sum = 0.0;
    Real result = 0.0;
    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(alpha[index]) || alpha[index] < 0.0) {
            throw PhysicsError(std::string(label) + " contains an invalid nonlinear weight");
        }
        sum += alpha[index];
        result += alpha[index] * candidates[index];
    }
    if (!std::isfinite(sum) || sum <= 0.0 || !std::isfinite(result)) {
        throw PhysicsError(std::string(label) + " has an invalid weight sum");
    }
    return result / sum;
}

Real weno_z_left_scaled(
    const std::array<Real, 6>& q,
    Real scale,
    const WcnsParameters& parameters)
{
    parameters.validate();
    const auto candidates3 = weno_candidates(q);
    const auto beta = weno_smoothness(q, scale);
    const Real tau5 = std::abs(beta[0] - beta[2]);
    constexpr std::array<Real, 3> optimal {{1.0 / 16.0, 10.0 / 16.0, 5.0 / 16.0}};
    std::array<Real, 4> alpha {};
    std::array<Real, 4> candidates {};
    for (std::size_t index = 0; index < 3; ++index) {
        alpha[index] = optimal[index]
            * (1.0 + std::pow(
                tau5 / (parameters.epsilon + beta[index]),
                parameters.weno_z_power));
        candidates[index] = candidates3[index];
    }
    return checked_weighted_sum(alpha, candidates, 3, "WENO-Z");
}

Real mdcd_linear_left(
    const std::array<Real, 6>& q,
    const WcnsParameters& parameters)
{
    parameters.validate();
    const Real dispersion = parameters.mdcd_dispersion;
    const Real dissipation = parameters.mdcd_dissipation;
    const std::array<Real, 6> coefficients {{
        3.0 * (dispersion + dissipation) / 8.0,
        (-18.0 * dispersion - 30.0 * dissipation - 1.0) / 16.0,
        (12.0 * dispersion + 60.0 * dissipation + 9.0) / 16.0,
        (12.0 * dispersion - 60.0 * dissipation + 9.0) / 16.0,
        (-18.0 * dispersion + 30.0 * dissipation - 1.0) / 16.0,
        3.0 * (dispersion - dissipation) / 8.0,
    }};
    Real result = 0.0;
    for (std::size_t index = 0; index < q.size(); ++index) {
        result += coefficients[index] * q[index];
    }
    if (!std::isfinite(result)) throw PhysicsError("MDCD-LINEAR result is non-finite");
    return result;
}

Real mdcd_sensor(const std::array<Real, 6>& q, Real scale, const WcnsParameters& parameters)
{
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("MDCD sensor scale must be positive and finite");
    }
    std::array<Real, 6> f {};
    for (std::size_t index = 0; index < q.size(); ++index) f[index] = q[index] / scale;
    const Real a1 = std::abs(f[2] - f[1]) + std::abs(f[2] - 2.0 * f[1] + f[0]);
    const Real b1 = std::abs(f[2] - f[3]) + std::abs(f[2] - 2.0 * f[3] + f[4]);
    const Real a2 = std::abs(f[3] - f[2]) + std::abs(f[3] - 2.0 * f[2] + f[1]);
    const Real b2 = std::abs(f[3] - f[4]) + std::abs(f[3] - 2.0 * f[4] + f[5]);
    const auto psi = [&](Real a, Real b) {
        return (2.0 * a * b + parameters.mdcd_sensor_epsilon)
            / (a * a + b * b + parameters.mdcd_sensor_epsilon);
    };
    const Real result = std::min(psi(a1, b1), psi(a2, b2));
    if (!std::isfinite(result)) throw PhysicsError("MDCD sensor is non-finite");
    return result;
}

Real mdcd_hybrid_left_scaled(
    const std::array<Real, 6>& q,
    Real scale,
    const WcnsParameters& parameters)
{
    parameters.validate();
    if (mdcd_sensor(q, scale, parameters) > parameters.mdcd_sensor_threshold) {
        return mdcd_linear_left(q, parameters);
    }
    const auto candidates3 = weno_candidates(q);
    const std::array<Real, 4> candidates {{
        candidates3[0], candidates3[1], candidates3[2],
        (15.0 * q[3] - 10.0 * q[4] + 3.0 * q[5]) / 8.0,
    }};
    const auto beta3 = weno_smoothness(q, scale);
    std::array<Real, 4> beta {{beta3[0], beta3[1], beta3[2],
        mdcd_six_point_smoothness(q, scale)}};
    for (auto& value : beta) {
        if (!std::isfinite(value) || value < -1.0e-13) {
            throw PhysicsError("MDCD smoothness indicator is invalid");
        }
        value = std::max(0.0, value);
    }
    const Real dispersion = parameters.mdcd_dispersion;
    const Real dissipation = parameters.mdcd_dissipation;
    const std::array<Real, 4> optimal {{
        1.5 * (dispersion + dissipation),
        0.5 - 1.5 * (dispersion - 3.0 * dissipation),
        0.5 - 1.5 * (dispersion + 3.0 * dissipation),
        1.5 * (dispersion - dissipation),
    }};
    const Real tau6 = std::abs(beta[3] - (beta[0] + 4.0 * beta[1] + beta[2]) / 6.0);
    std::array<Real, 4> alpha {};
    for (std::size_t index = 0; index < alpha.size(); ++index) {
        alpha[index] = optimal[index] * square(
            parameters.mdcd_weight_bias
            + tau6 / (beta[index] + parameters.mdcd_weight_epsilon));
    }
    return checked_weighted_sum(alpha, candidates, 4, "MDCD-HYBRID");
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
    if (!std::isfinite(epsilon) || epsilon <= 0.0 || nonlinear_power <= 0
        || weno_z_power <= 0 || !std::isfinite(mdcd_dispersion)
        || !std::isfinite(mdcd_dissipation)
        || !std::isfinite(mdcd_sensor_epsilon)
        || !std::isfinite(mdcd_sensor_threshold)
        || !std::isfinite(mdcd_weight_bias)
        || !std::isfinite(mdcd_weight_epsilon)) {
        throw std::invalid_argument("WCNS parameters require positive epsilon and power");
    }
    if (mdcd_dispersion <= 0.0 || mdcd_dissipation < 0.0
        || mdcd_dissipation >= mdcd_dispersion
        || 3.0 * mdcd_dispersion + 9.0 * mdcd_dissipation >= 1.0
        || mdcd_sensor_epsilon <= 0.0 || mdcd_sensor_threshold <= 0.0
        || mdcd_sensor_threshold >= 1.0 || mdcd_weight_bias <= 0.0
        || mdcd_weight_epsilon <= 0.0) {
        throw std::invalid_argument("MDCD parameters are outside the frozen Stage-L domain");
    }
}

void ReconstructionRegistry::register_scheme(std::string name, Factory factory)
{
    if (!valid_algorithm_name(name) || !factory) {
        throw std::invalid_argument("reconstruction registration has an invalid name or factory");
    }
    auto probe = factory();
    if (!probe || probe->name() != name) {
        throw std::invalid_argument("reconstruction factory name does not match its registry key");
    }
    const auto requirement = probe->stencil_requirement();
    if (requirement.first_offset != -2 || requirement.point_count != 6
        || requirement.halo_width != 3) {
        throw std::invalid_argument("Stage-L reconstruction must declare the frozen six-point stencil");
    }
    if (!factories_.emplace(std::move(name), std::move(factory)).second) {
        throw std::invalid_argument("duplicate reconstruction registration");
    }
}

bool ReconstructionRegistry::contains(std::string_view name) const noexcept
{
    return factories_.find(std::string(name)) != factories_.end();
}

std::unique_ptr<IReconstructionScheme> ReconstructionRegistry::create(
    std::string_view name) const
{
    const auto iterator = factories_.find(std::string(name));
    if (iterator == factories_.end()) {
        throw std::invalid_argument("unknown reconstruction scheme: " + std::string(name));
    }
    auto result = iterator->second();
    if (!result || result->name() != iterator->first) {
        throw std::logic_error("registered reconstruction factory returned an invalid strategy");
    }
    return result;
}

std::vector<std::string> ReconstructionRegistry::names() const
{
    std::vector<std::string> result;
    result.reserve(factories_.size());
    for (const auto& [name, factory] : factories_) {
        static_cast<void>(factory);
        result.push_back(name);
    }
    std::sort(result.begin(), result.end());
    return result;
}

ReconstructionRegistry ReconstructionRegistry::with_builtins()
{
    ReconstructionRegistry result;
    result.register_scheme("linear5", [] { return std::make_unique<Linear5Scheme>(); });
    result.register_scheme("weno_js", [] { return std::make_unique<WenoJsScheme>(); });
    result.register_scheme("weno_z", [] { return std::make_unique<WenoZScheme>(); });
    result.register_scheme("mdcd_linear", [] { return std::make_unique<MdcdLinearScheme>(); });
    result.register_scheme("mdcd_hybrid", [] { return std::make_unique<MdcdHybridScheme>(); });
    return result;
}

std::string_view reconstruction_name(ReconstructionKind kind)
{
    switch (kind) {
    case ReconstructionKind::Linear5: return "linear5";
    case ReconstructionKind::WenoJs: return "weno_js";
    case ReconstructionKind::WenoZ: return "weno_z";
    case ReconstructionKind::MdcdLinear: return "mdcd_linear";
    case ReconstructionKind::MdcdHybrid: return "mdcd_hybrid";
    }
    throw std::invalid_argument("unknown reconstruction kind");
}

void ReconstructionConfig::validate(const ReconstructionRegistry& registry) const
{
    validate();
    if (!registry.contains(scheme)) {
        throw std::invalid_argument("unknown reconstruction scheme: " + scheme);
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
    if (!valid_algorithm_name(scheme)) {
        throw std::invalid_argument("invalid reconstruction scheme name: " + scheme);
    }
    switch (variables) {
    case ReconstructionVariables::Conservative:
    case ReconstructionVariables::Primitive:
    case ReconstructionVariables::Characteristic: break;
    default: throw std::invalid_argument("unknown reconstruction variable set");
    }
}

std::string ReconstructionConfig::summary() const
{
    validate();
    const char* variable_name = "characteristic";
    if (variables == ReconstructionVariables::Conservative) {
        variable_name = "conservative";
    } else if (variables == ReconstructionVariables::Primitive) {
        variable_name = "primitive";
    }
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "reconstruction=" << scheme
           << ";variables=" << variable_name
           << ";epsilon=" << nonlinear.epsilon
           << ";js_power=" << nonlinear.nonlinear_power
           << ";z_power=" << nonlinear.weno_z_power
           << ";mdcd_dispersion=" << nonlinear.mdcd_dispersion
           << ";mdcd_dissipation=" << nonlinear.mdcd_dissipation
           << ";mdcd_sensor_epsilon=" << nonlinear.mdcd_sensor_epsilon
           << ";mdcd_sensor_threshold=" << nonlinear.mdcd_sensor_threshold
           << ";mdcd_weight_bias=" << nonlinear.mdcd_weight_bias
           << ";mdcd_weight_epsilon=" << nonlinear.mdcd_weight_epsilon;
    return output.str();
}

std::string ReconstructionConfig::restart_signature() const
{
    return "reconstruction_v2;" + summary();
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

ScalarFaceStates weno_z_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters)
{
    const auto values = checked_stencil(ScalarStencilView(stencil));
    return {
        weno_z_left_scaled(values, scale, parameters),
        weno_z_left_scaled(orient_stencil(values, TraceSide::Right), scale, parameters),
    };
}

ScalarFaceStates mdcd_linear_reconstruct(
    const std::array<Real, 6>& stencil,
    const WcnsParameters& parameters)
{
    const auto values = checked_stencil(ScalarStencilView(stencil));
    return {
        mdcd_linear_left(values, parameters),
        mdcd_linear_left(orient_stencil(values, TraceSide::Right), parameters),
    };
}

ScalarFaceStates mdcd_hybrid_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters)
{
    const auto values = checked_stencil(ScalarStencilView(stencil));
    return {
        mdcd_hybrid_left_scaled(values, scale, parameters),
        mdcd_hybrid_left_scaled(
            orient_stencil(values, TraceSide::Right), scale, parameters),
    };
}

Real mdcd_six_point_smoothness(
    const std::array<Real, 6>& stencil,
    Real scale)
{
    const auto values = checked_stencil(ScalarStencilView(stencil));
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("MDCD smoothness scale must be positive and finite");
    }
    std::array<Real, 6> f {};
    Real maximum = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
        f[index] = values[index] / scale;
        maximum = std::max(maximum, std::abs(f[index]));
    }
    const Real numerator =
        271779.0 * square(f[0])
        + f[0] * (-2380800.0 * f[1] + 4086352.0 * f[2]
            - 3462252.0 * f[3] + 1458762.0 * f[4] - 245620.0 * f[5])
        + f[1] * (5653317.0 * f[1] - 20427884.0 * f[2]
            + 17905032.0 * f[3] - 7727988.0 * f[4] + 1325006.0 * f[5])
        + f[2] * (19510972.0 * f[2] - 35817664.0 * f[3]
            + 15929912.0 * f[4] - 2792660.0 * f[5])
        + f[3] * (17195652.0 * f[3] - 15880404.0 * f[4]
            + 2863984.0 * f[5])
        + f[4] * (3824847.0 * f[4] - 1429976.0 * f[5])
        + 139633.0 * square(f[5]);
    const Real result = numerator / 120960.0;
    const Real tolerance = 1.0e-11 * std::max(1.0, maximum * maximum);
    if (!std::isfinite(result) || result < -tolerance) {
        throw PhysicsError("MDCD six-point smoothness is negative or non-finite");
    }
    return std::max(0.0, result);
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
    const auto registry = ReconstructionRegistry::with_builtins();
    config.validate(registry);
    if (config.variables == ReconstructionVariables::Characteristic) {
        throw std::invalid_argument("characteristic reconstruction requires the Stage-L L3 path");
    }
    const auto& source = config.variables == ReconstructionVariables::Conservative
        ? conservative : pressure_primitive_field;
    if (source.components() != euler_components || source.ghost_width() < 3) {
        throw std::invalid_argument("thermodynamic reconstruction requires five components and three ghost layers");
    }
    std::array<Real, euler_components> left {};
    std::array<Real, euler_components> right {};
    const auto reconstruct = [&](const IReconstructionScheme& scheme) {
        for (int component = 0; component < euler_components; ++component) {
            const auto stencil = component_stencil(source, axis, face, component);
            ReconstructionContext context;
            context.scale = std::max(
                config.scaling.component[static_cast<std::size_t>(component)],
                config.scaling.scale_floor);
            context.parameters = config.nonlinear;
            left[static_cast<std::size_t>(component)] = scheme.reconstruct_scalar(
                stencil, TraceSide::Left, context);
            right[static_cast<std::size_t>(component)] = scheme.reconstruct_scalar(
                stencil, TraceSide::Right, context);
        }
    };

    EulerFaceStates result;
    if (config.scheme != "linear5") {
        const auto selected = registry.create(config.scheme);
        reconstruct(*selected);
        ++diagnostics.nonlinear_faces;
        if (try_convert(left, right, config, gas, reference, dimension, result)) {
            return result;
        }
        ++diagnostics.linear_fallbacks;
    }
    const auto linear = registry.create("linear5");
    reconstruct(*linear);
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
