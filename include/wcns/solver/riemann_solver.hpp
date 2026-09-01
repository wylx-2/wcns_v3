#pragma once

#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <string>

namespace wcns {

enum class RiemannSolverKind {
    Rusanov,
};

class RiemannSolver {
public:
    explicit RiemannSolver(RiemannSolverKind kind = RiemannSolverKind::Rusanov)
        : kind_(kind)
    {
    }

    [[nodiscard]] RiemannSolverKind kind() const noexcept { return kind_; }
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;

    [[nodiscard]] ConservativeState flux(
        const PressurePrimitiveState& left,
        const PressurePrimitiveState& right,
        Normal3 unit_normal,
        const GasModel& gas,
        const NumericalFloors& floors) const;

private:
    RiemannSolverKind kind_;
};

} // namespace wcns
