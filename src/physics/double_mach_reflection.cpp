#include <wcns/physics/double_mach_reflection.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace wcns {
namespace {

constexpr Real classic_gamma = 1.4;
constexpr Real sqrt_three = 1.732050807568877293527446341505872367;

void require_finite(Real value, const char* label)
{
    if (!std::isfinite(value)) {
        throw PhysicsConfigurationError(
            std::string("double-Mach-reflection ") + label + " must be finite");
    }
}

} // namespace

DoubleMachReflection::DoubleMachReflection(Real shock_foot)
    : shock_foot_(shock_foot)
{
    require_finite(shock_foot_, "shock foot");
}

void DoubleMachReflection::validate(Real gamma, int dimension) const
{
    require_finite(shock_foot_, "shock foot");
    if (!std::isfinite(gamma)
        || std::abs(gamma - classic_gamma)
            > 64.0 * std::numeric_limits<Real>::epsilon()) {
        throw PhysicsConfigurationError(
            "classical double-Mach reflection requires gas.gamma=1.4");
    }
    if (dimension != 2) {
        throw PhysicsConfigurationError(
            "classical double-Mach reflection requires a two-dimensional mesh");
    }
}

Real DoubleMachReflection::shock_x(Real y, Real time) const
{
    require_finite(y, "y coordinate");
    require_finite(time, "time");
    if (time < 0.0) {
        throw PhysicsConfigurationError(
            "double-Mach-reflection time must be non-negative");
    }
    // x = x0 + y/tan(60 deg) + (Mach/sin(60 deg))*t.
    return shock_foot_ + y / sqrt_three + 20.0 * time / sqrt_three;
}

bool DoubleMachReflection::is_post_shock(Real x, Real y, Real time) const
{
    require_finite(x, "x coordinate");
    return x <= shock_x(y, time);
}

PressurePrimitiveState DoubleMachReflection::upstream_state() const noexcept
{
    return {1.4, 0.0, 0.0, 0.0, 1.0};
}

PressurePrimitiveState DoubleMachReflection::post_shock_state() const noexcept
{
    return {8.0, 8.25 * sqrt_three / 2.0, -4.125, 0.0, 116.5};
}

PressurePrimitiveState DoubleMachReflection::exact_state(
    Real x, Real y, Real time) const
{
    return is_post_shock(x, y, time)
        ? post_shock_state() : upstream_state();
}

std::string DoubleMachReflection::restart_signature() const
{
    std::ostringstream result;
    result << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "double_mach_reflection_v1;shock_foot=" << shock_foot_
           << ";gamma=1.4;incident_mach=10;shock_angle_deg=60";
    return result.str();
}

} // namespace wcns
