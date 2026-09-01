#include <wcns/solver/transport_model.hpp>

#include <wcns/solver/euler.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace wcns {
namespace {

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

std::string number(Real value)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
    return stream.str();
}

} // namespace

void TransportConfig::validate() const
{
    if (!positive_finite(prandtl)) {
        throw std::invalid_argument("transport Prandtl number must be positive and finite");
    }
    std::visit([](const auto& law) {
        using Law = std::decay_t<decltype(law)>;
        if constexpr (std::is_same_v<Law, ConstantViscosity>) {
            if (!positive_finite(law.viscosity_ratio)) {
                throw std::invalid_argument(
                    "constant viscosity ratio must be positive and finite");
            }
        } else {
            if (!positive_finite(law.reference_viscosity_ratio)
                || !positive_finite(law.constant_temperature_ratio)) {
                throw std::invalid_argument(
                    "Sutherland viscosity ratios must be positive and finite");
            }
        }
    }, viscosity);
}

std::string TransportConfig::summary() const
{
    validate();
    std::string result = "transport Pr=" + number(prandtl);
    std::visit([&](const auto& law) {
        using Law = std::decay_t<decltype(law)>;
        if constexpr (std::is_same_v<Law, ConstantViscosity>) {
            result += " law=constant mu_const_over_mu_ref="
                + number(law.viscosity_ratio);
        } else {
            result += " law=sutherland mu_Tref_over_mu_ref="
                + number(law.reference_viscosity_ratio)
                + " S_over_Tref=" + number(law.constant_temperature_ratio);
        }
    }, viscosity);
    return result;
}

std::string TransportConfig::restart_signature() const
{
    return "transport_v1;" + summary();
}

TransportModel::TransportModel(TransportConfig config)
    : config_(std::move(config))
{
    config_.validate();
}

Real TransportModel::viscosity(Real temperature) const
{
    if (!positive_finite(temperature)) {
        throw PhysicsError("transport temperature must be positive and finite");
    }
    const Real result = std::visit([temperature](const auto& law) {
        using Law = std::decay_t<decltype(law)>;
        if constexpr (std::is_same_v<Law, ConstantViscosity>) {
            return law.viscosity_ratio;
        } else {
            return law.reference_viscosity_ratio
                * std::pow(temperature, 1.5)
                * (1.0 + law.constant_temperature_ratio)
                / (temperature + law.constant_temperature_ratio);
        }
    }, config_.viscosity);
    if (!positive_finite(result)) {
        throw PhysicsError("transport viscosity is not positive and finite");
    }
    return result;
}

Real TransportModel::thermal_coefficient(
    Real temperature,
    const GasModel& gas,
    const ReferenceScales& reference) const
{
    const Real denominator = (gas.gamma() - 1.0)
        * reference.mach() * reference.mach() * config_.prandtl;
    if (!positive_finite(denominator)) {
        throw PhysicsError("thermal coefficient denominator is invalid");
    }
    return viscosity(temperature) / denominator;
}

} // namespace wcns
