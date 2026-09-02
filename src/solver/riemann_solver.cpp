#include <wcns/solver/riemann_solver.hpp>

#include <algorithm>
#include <cmath>
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
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const Real spectral_radius = std::max(
            std::abs(normal_velocity(left, normal)) + sound_speed(left, ideal),
            std::abs(normal_velocity(right, normal)) + sound_speed(right, ideal));
        return {
            rusanov_flux(left, right, normal, ideal),
            spectral_radius,
            "rusanov",
            "rusanov",
            RiemannFallbackReason::None,
        };
    }
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

RiemannSolverRegistry RiemannSolverRegistry::with_builtins()
{
    RiemannSolverRegistry result;
    result.register_solver("rusanov", [] { return std::make_unique<RusanovStrategy>(); });
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

RiemannSolver::RiemannSolver(RiemannSolverKind kind)
{
    const auto registry = RiemannSolverRegistry::with_builtins();
    implementation_ = registry.create(riemann_solver_name(kind));
}

RiemannSolver::RiemannSolver(
    std::string_view name,
    const RiemannSolverRegistry& registry)
    : implementation_(registry.create(name))
{
}

std::string_view RiemannSolver::name() const noexcept
{
    return implementation_->name();
}

std::string RiemannSolver::summary() const
{
    return "riemann_solver=" + std::string(name());
}

std::string RiemannSolver::restart_signature() const
{
    return "riemann_v2;" + summary();
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
