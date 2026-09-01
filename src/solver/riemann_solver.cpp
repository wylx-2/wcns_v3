#include <wcns/solver/riemann_solver.hpp>

#include <cmath>
#include <stdexcept>

namespace wcns {

std::string RiemannSolver::summary() const
{
    if (kind_ != RiemannSolverKind::Rusanov) {
        throw std::invalid_argument("unknown Riemann solver kind");
    }
    return "riemann_solver=rusanov";
}

std::string RiemannSolver::restart_signature() const
{
    return "riemann_v1;" + summary();
}

ConservativeState RiemannSolver::flux(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 unit_normal,
    const GasModel& gas,
    const NumericalFloors& floors) const
{
    static_cast<void>(summary());
    floors.validate();
    const Real norm = std::sqrt(
        unit_normal.x * unit_normal.x + unit_normal.y * unit_normal.y
        + unit_normal.z * unit_normal.z);
    if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsError("Riemann solver requires a unit normal");
    }
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    return rusanov_flux(left, right, unit_normal, ideal);
}

} // namespace wcns
