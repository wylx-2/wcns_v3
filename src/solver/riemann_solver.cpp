#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
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

Normal3 checked_unit_normal(Normal3 normal)
{
    const Real norm = std::sqrt(
        normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsError("Riemann solver requires a unit normal");
    }
    return normal;
}

Real normal_velocity(const PressurePrimitiveState& state, Normal3 normal)
{
    return state[1] * normal.x + state[2] * normal.y + state[3] * normal.z;
}

bool finite_state(const ConservativeState& state)
{
    return std::all_of(state.begin(), state.end(), [](Real value) {
        return std::isfinite(value);
    });
}

struct RoeAverage {
    Normal3 velocity {};
    Real normal_velocity = 0.0;
    Real enthalpy = 0.0;
    Real sound_speed = 0.0;
};

RoeAverage roe_average(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 normal,
    const IdealGas& gas)
{
    const auto left_conservative = to_conservative(left, gas);
    const auto right_conservative = to_conservative(right, gas);
    const Real root_left = std::sqrt(left[0]);
    const Real root_right = std::sqrt(right[0]);
    const Real denominator = root_left + root_right;
    if (!std::isfinite(denominator) || denominator <= 0.0) {
        throw PhysicsError("Roe average has an invalid density denominator");
    }
    const auto average = [&](Real left_value, Real right_value) {
        return (root_left * left_value + root_right * right_value) / denominator;
    };
    RoeAverage result;
    result.velocity = {
        average(left[1], right[1]),
        average(left[2], right[2]),
        average(left[3], right[3]),
    };
    const Real left_enthalpy = (left_conservative[4] + left[4]) / left[0];
    const Real right_enthalpy = (right_conservative[4] + right[4]) / right[0];
    result.enthalpy = average(left_enthalpy, right_enthalpy);
    result.normal_velocity = result.velocity.x * normal.x
        + result.velocity.y * normal.y + result.velocity.z * normal.z;
    const Real speed_squared = result.velocity.x * result.velocity.x
        + result.velocity.y * result.velocity.y
        + result.velocity.z * result.velocity.z;
    const Real sound_squared
        = (gas.gamma - 1.0) * (result.enthalpy - 0.5 * speed_squared);
    if (!std::isfinite(sound_squared) || sound_squared <= 0.0) {
        throw PhysicsError("Roe average has an invalid sound speed");
    }
    result.sound_speed = std::sqrt(sound_squared);
    return result;
}

RiemannResult rusanov_result(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 normal,
    const GasModel& gas,
    const NumericalFloors& floors,
    std::string requested,
    RiemannFallbackReason reason)
{
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const Real spectral_radius = std::max(
        std::abs(normal_velocity(left, normal)) + sound_speed(left, ideal),
        std::abs(normal_velocity(right, normal)) + sound_speed(right, ideal));
    return {
        rusanov_flux(left, right, normal, ideal),
        spectral_radius,
        std::move(requested),
        "rusanov",
        reason,
    };
}

class RusanovStrategy final : public IRiemannSolver {
public:
    [[nodiscard]] std::string_view name() const noexcept override
    {
        return "rusanov";
    }

    [[nodiscard]] RiemannResult solve(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        return rusanov_result(
            left, right, normal, gas, floors,
            "rusanov", RiemannFallbackReason::None);
    }
};

class HllcStrategy final : public IRiemannSolver {
public:
    explicit HllcStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "hllc"; }

    [[nodiscard]] RiemannResult solve(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto left_conservative = to_conservative(left, ideal);
        const auto right_conservative = to_conservative(right, ideal);
        const auto left_flux = euler_flux(left, normal, ideal);
        const auto right_flux = euler_flux(right, normal, ideal);
        const Real un_left = normal_velocity(left, normal);
        const Real un_right = normal_velocity(right, normal);
        const Real sound_left = sound_speed(left, ideal);
        const Real sound_right = sound_speed(right, ideal);
        RoeAverage average;
        try {
            average = roe_average(left, right, normal, ideal);
        } catch (const PhysicsError&) {
            return rusanov_result(
                left, right, normal, gas, floors,
                "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const Real speed_left = std::min(
            un_left - sound_left,
            average.normal_velocity - average.sound_speed);
        const Real speed_right = std::max(
            un_right + sound_right,
            average.normal_velocity + average.sound_speed);
        const Real spectral_radius
            = std::max(std::abs(speed_left), std::abs(speed_right));
        if (!std::isfinite(speed_left) || !std::isfinite(speed_right)
            || speed_left >= speed_right) {
            return rusanov_result(
                left, right, normal, gas, floors,
                "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        if (speed_left >= 0.0) {
            return {left_flux, spectral_radius, "hllc", "hllc",
                RiemannFallbackReason::None};
        }
        if (speed_right <= 0.0) {
            return {right_flux, spectral_radius, "hllc", "hllc",
                RiemannFallbackReason::None};
        }
        const Real left_term = left[0] * (speed_left - un_left);
        const Real right_term = right[0] * (speed_right - un_right);
        const Real denominator = left_term - right_term;
        const Real denominator_scale
            = std::max({1.0, std::abs(left_term), std::abs(right_term)});
        if (!std::isfinite(denominator)
            || std::abs(denominator)
                <= parameters_.denominator_tolerance * denominator_scale) {
            return rusanov_result(
                left, right, normal, gas, floors,
                "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const Real speed_middle = (right[4] - left[4]
            + left_term * un_left - right_term * un_right) / denominator;
        if (!std::isfinite(speed_middle)
            || speed_middle <= speed_left || speed_middle >= speed_right) {
            return rusanov_result(
                left, right, normal, gas, floors,
                "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const auto star_flux = [&](const PressurePrimitiveState& state,
                                   const ConservativeState& conservative,
                                   const ConservativeState& physical_flux,
                                   Real wave_speed,
                                   Real un) {
            const Real star_denominator = wave_speed - speed_middle;
            const Real scale = std::max({1.0, std::abs(wave_speed), std::abs(speed_middle)});
            if (!std::isfinite(star_denominator)
                || std::abs(star_denominator)
                    <= parameters_.denominator_tolerance * scale) {
                throw PhysicsError("HLLC star-state denominator is invalid");
            }
            const Real star_density
                = state[0] * (wave_speed - un) / star_denominator;
            const Real star_pressure
                = state[4] + state[0] * (wave_speed - un) * (speed_middle - un);
            const Normal3 star_velocity {
                state[1] + (speed_middle - un) * normal.x,
                state[2] + (speed_middle - un) * normal.y,
                state[3] + (speed_middle - un) * normal.z,
            };
            ConservativeState star {{
                star_density,
                star_density * star_velocity.x,
                star_density * star_velocity.y,
                star_density * star_velocity.z,
                ((wave_speed - un) * conservative[4] - state[4] * un
                    + star_pressure * speed_middle) / star_denominator,
            }};
            if (!finite_state(star)) {
                throw PhysicsError("HLLC star state is non-finite");
            }
            static_cast<void>(to_primitive(star, ideal));
            ConservativeState flux = physical_flux;
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index] += wave_speed * (star[index] - conservative[index]);
            }
            if (!finite_state(flux)) throw PhysicsError("HLLC star flux is non-finite");
            return flux;
        };
        try {
            const auto flux = speed_middle >= 0.0
                ? star_flux(left, left_conservative, left_flux, speed_left, un_left)
                : star_flux(right, right_conservative, right_flux, speed_right, un_right);
            return {flux, spectral_radius, "hllc", "hllc",
                RiemannFallbackReason::None};
        } catch (const PhysicsError&) {
            return rusanov_result(
                left, right, normal, gas, floors,
                "hllc", RiemannFallbackReason::InvalidIntermediateState);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
};

class RoeStrategy final : public IRiemannSolver {
public:
    explicit RoeStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "roe"; }

    [[nodiscard]] RiemannResult solve(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto fallback = [&](RiemannFallbackReason reason) {
            auto result = HllcStrategy(parameters_).solve(
                left, right, normal, gas, floors);
            result.requested_solver = "roe";
            result.fallback_reason = reason;
            return result;
        };
        try {
            const auto left_conservative = to_conservative(left, ideal);
            const auto right_conservative = to_conservative(right, ideal);
            const auto left_flux = euler_flux(left, normal, ideal);
            const auto right_flux = euler_flux(right, normal, ideal);
            const auto average = roe_average(left, right, normal, ideal);
            const auto basis = make_roe_characteristic_basis(
                left, right, normal, gas, floors, 3);
            ConservativeState jump {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                jump[index] = right_conservative[index] - left_conservative[index];
            }
            const auto strengths = project_characteristic(jump, basis);
            std::array<Real, euler_components> eigenvalues {{
                average.normal_velocity - average.sound_speed,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity + average.sound_speed,
            }};
            const Real delta = parameters_.entropy_fix_coefficient
                * std::max({sound_speed(left, ideal), sound_speed(right, ideal),
                    average.sound_speed});
            ConservativeState scaled_strengths {};
            for (int wave = 0; wave < euler_components; ++wave) {
                Real magnitude = std::abs(eigenvalues[static_cast<std::size_t>(wave)]);
                if (magnitude < delta) {
                    magnitude = 0.5 * (magnitude * magnitude / delta + delta);
                }
                scaled_strengths[static_cast<std::size_t>(wave)]
                    = magnitude * strengths[static_cast<std::size_t>(wave)];
            }
            const auto dissipation = restore_characteristic(scaled_strengths, basis);
            ConservativeState flux {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index] = 0.5 * (left_flux[index] + right_flux[index])
                    - 0.5 * dissipation[index];
            }
            if (!finite_state(flux)) return fallback(RiemannFallbackReason::NonFiniteFlux);
            return {
                flux,
                std::abs(average.normal_velocity) + average.sound_speed,
                "roe",
                "roe",
                RiemannFallbackReason::None,
            };
        } catch (const PhysicsError&) {
            return fallback(RiemannFallbackReason::InvalidRoeAverage);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
};

} // namespace

void RiemannSolverRegistry::register_solver(std::string name, Factory factory)
{
    if (!valid_algorithm_name(name) || !factory) {
        throw std::invalid_argument("Riemann registration has an invalid name or factory");
    }
    auto probe = factory();
    if (!probe || probe->name() != name) {
        throw std::invalid_argument("Riemann factory name does not match its registry key");
    }
    if (!factories_.emplace(std::move(name), std::move(factory)).second) {
        throw std::invalid_argument("duplicate Riemann solver registration");
    }
}

bool RiemannSolverRegistry::contains(std::string_view name) const noexcept
{
    return factories_.find(std::string(name)) != factories_.end();
}

std::unique_ptr<IRiemannSolver> RiemannSolverRegistry::create(
    std::string_view name) const
{
    const auto iterator = factories_.find(std::string(name));
    if (iterator == factories_.end()) {
        throw std::invalid_argument("unknown Riemann solver: " + std::string(name));
    }
    auto result = iterator->second();
    if (!result || result->name() != iterator->first) {
        throw std::logic_error("registered Riemann factory returned an invalid strategy");
    }
    return result;
}

std::vector<std::string> RiemannSolverRegistry::names() const
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

void RiemannSolverParameters::validate() const
{
    if (!std::isfinite(entropy_fix_coefficient)
        || entropy_fix_coefficient <= 0.0) {
        throw std::invalid_argument(
            "Roe entropy-fix coefficient must be positive and finite");
    }
    if (!std::isfinite(denominator_tolerance)
        || denominator_tolerance <= 0.0 || denominator_tolerance >= 1.0) {
        throw std::invalid_argument(
            "Riemann denominator tolerance must be finite and lie in (0,1)");
    }
}

void RiemannConfig::validate() const
{
    if (!valid_algorithm_name(scheme)) {
        throw std::invalid_argument("Riemann scheme name is invalid");
    }
    parameters.validate();
}

void RiemannConfig::validate(const RiemannSolverRegistry& registry) const
{
    validate();
    if (!registry.contains(scheme)) {
        throw std::invalid_argument("unknown Riemann solver: " + scheme);
    }
}

std::string RiemannConfig::summary() const
{
    validate();
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "riemann_solver=" << scheme
           << ";roe_entropy_fix=" << parameters.entropy_fix_coefficient
           << ";riemann_denominator_tolerance="
           << parameters.denominator_tolerance;
    return stream.str();
}

std::string RiemannConfig::restart_signature() const
{
    return "riemann_config_v1;" + summary();
}

RiemannSolverRegistry RiemannSolverRegistry::with_builtins(
    const RiemannSolverParameters& parameters)
{
    parameters.validate();
    RiemannSolverRegistry result;
    result.register_solver("rusanov", [] { return std::make_unique<RusanovStrategy>(); });
    result.register_solver("hllc", [parameters] {
        return std::make_unique<HllcStrategy>(parameters);
    });
    result.register_solver("roe", [parameters] {
        return std::make_unique<RoeStrategy>(parameters);
    });
    return result;
}

std::string_view riemann_solver_name(RiemannSolverKind kind)
{
    switch (kind) {
    case RiemannSolverKind::Rusanov: return "rusanov";
    case RiemannSolverKind::Hllc: return "hllc";
    case RiemannSolverKind::Roe: return "roe";
    }
    throw std::invalid_argument("unknown Riemann solver kind");
}

RiemannSolver::RiemannSolver(
    RiemannSolverKind kind,
    RiemannSolverParameters parameters)
    : parameters_(parameters)
{
    parameters_.validate();
    const auto registry = RiemannSolverRegistry::with_builtins(parameters_);
    implementation_ = registry.create(riemann_solver_name(kind));
}

RiemannSolver::RiemannSolver(
    std::string_view name,
    const RiemannSolverRegistry& registry,
    RiemannSolverParameters parameters)
    : implementation_(registry.create(name))
    , parameters_(parameters)
{
    parameters_.validate();
}

std::string_view RiemannSolver::name() const noexcept
{
    return implementation_->name();
}

std::string RiemannSolver::summary() const
{
    RiemannConfig config;
    config.scheme = std::string(name());
    config.parameters = parameters_;
    return config.summary();
}

std::string RiemannSolver::restart_signature() const
{
    return "riemann_v3;" + summary();
}

RiemannResult RiemannSolver::solve(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 unit_normal,
    const GasModel& gas,
    const NumericalFloors& floors) const
{
    static_cast<void>(checked_unit_normal(unit_normal));
    auto result = implementation_->solve(
        left, right, unit_normal, gas, floors);
    if (result.requested_solver.empty()) result.requested_solver = std::string(name());
    if (result.used_solver.empty()) result.used_solver = result.requested_solver;
    if (result.requested_solver != name()) {
        throw std::logic_error("Riemann result requested-solver diagnostic is inconsistent");
    }
    if (!std::isfinite(result.spectral_radius) || result.spectral_radius < 0.0) {
        throw PhysicsError("Riemann solver returned an invalid spectral radius");
    }
    for (const auto value : result.flux_per_unit_area) {
        if (!std::isfinite(value)) {
            throw PhysicsError("Riemann solver returned a non-finite flux");
        }
    }
    return result;
}

ConservativeState RiemannSolver::flux(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 unit_normal,
    const GasModel& gas,
    const NumericalFloors& floors) const
{
    return solve(left, right, unit_normal, gas, floors).flux_per_unit_area;
}

} // namespace wcns
