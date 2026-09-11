#pragma once

#include <wcns/physics/thermodynamics.hpp>

#include <string>

namespace wcns {

// Classical Woodward-Colella double-Mach-reflection data.  The fixed states
// are the gamma=1.4, incident-Mach-10 benchmark values.
class DoubleMachReflection {
public:
    explicit DoubleMachReflection(Real shock_foot = 1.0 / 6.0);

    void validate(Real gamma, int dimension) const;
    [[nodiscard]] Real shock_foot() const noexcept { return shock_foot_; }
    [[nodiscard]] Real shock_x(Real y, Real time) const;
    [[nodiscard]] bool is_post_shock(Real x, Real y, Real time) const;
    [[nodiscard]] PressurePrimitiveState upstream_state() const noexcept;
    [[nodiscard]] PressurePrimitiveState post_shock_state() const noexcept;
    [[nodiscard]] PressurePrimitiveState exact_state(Real x, Real y, Real time) const;
    [[nodiscard]] std::string restart_signature() const;

private:
    Real shock_foot_ = 1.0 / 6.0;
};

} // namespace wcns
