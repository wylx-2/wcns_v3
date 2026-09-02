#pragma once

#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class RiemannSolverKind {
    Rusanov,
    Hllc,
    Roe,
};

enum class RiemannFallbackReason {
    None,
    InvalidWaveSpeed,
    InvalidIntermediateState,
    InvalidRoeAverage,
    NonFiniteFlux,
};

struct RiemannResult {
    ConservativeState flux_per_unit_area {};
    Real spectral_radius = 0.0;
    std::string requested_solver;
    std::string used_solver;
    RiemannFallbackReason fallback_reason = RiemannFallbackReason::None;
};

class IRiemannSolver {
public:
    virtual ~IRiemannSolver() = default;
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;
    [[nodiscard]] virtual RiemannResult solve(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const = 0;
};

class RiemannSolverRegistry {
public:
    using Factory = std::function<std::unique_ptr<IRiemannSolver>()>;

    void register_solver(std::string name, Factory factory);
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] std::unique_ptr<IRiemannSolver> create(
        std::string_view name) const;
    [[nodiscard]] std::vector<std::string> names() const;

    [[nodiscard]] static RiemannSolverRegistry with_builtins();

private:
    std::unordered_map<std::string, Factory> factories_;
};

[[nodiscard]] std::string_view riemann_solver_name(RiemannSolverKind kind);

class RiemannSolver {
public:
    explicit RiemannSolver(RiemannSolverKind kind = RiemannSolverKind::Rusanov);
    explicit RiemannSolver(
        std::string_view name,
        const RiemannSolverRegistry& registry);

    [[nodiscard]] std::string_view name() const noexcept;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;

    [[nodiscard]] RiemannResult solve(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const;

    [[nodiscard]] ConservativeState flux(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const;

private:
    std::unique_ptr<IRiemannSolver> implementation_;
};

} // namespace wcns
