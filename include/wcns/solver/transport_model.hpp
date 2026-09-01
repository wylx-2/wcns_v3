#pragma once

#include <wcns/physics/thermodynamics.hpp>

#include <string>
#include <variant>

namespace wcns {

struct ConstantViscosity {
    Real viscosity_ratio = 1.0;
};

struct SutherlandViscosity {
    Real reference_viscosity_ratio = 1.0;
    Real constant_temperature_ratio = 0.383;
};

using ViscosityLaw = std::variant<ConstantViscosity, SutherlandViscosity>;

struct TransportConfig {
    Real prandtl = 0.72;
    ViscosityLaw viscosity = ConstantViscosity {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

class TransportModel {
public:
    explicit TransportModel(TransportConfig config);

    [[nodiscard]] const TransportConfig& config() const noexcept { return config_; }
    [[nodiscard]] Real viscosity(Real temperature) const;
    [[nodiscard]] Real thermal_coefficient(
        Real temperature,
        const GasModel& gas,
        const ReferenceScales& reference) const;

private:
    TransportConfig config_;
};

} // namespace wcns
