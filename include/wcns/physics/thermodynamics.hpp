#pragma once

#include <wcns/core/types.hpp>

#include <array>
#include <optional>
#include <stdexcept>
#include <string>

namespace wcns {

inline constexpr int fluid_components = 5;
inline constexpr Real universal_gas_constant = 8.314;

class PhysicsConfigurationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct GasModelInput {
    Real gamma = 1.4;
    std::optional<Real> molar_mass;
    std::optional<Real> specific_gas_constant;
};

class GasModel {
public:
    [[nodiscard]] static GasModel from_input(const GasModelInput& input);

    [[nodiscard]] Real gamma() const noexcept { return gamma_; }
    [[nodiscard]] Real molar_mass() const noexcept { return molar_mass_; }
    [[nodiscard]] Real specific_gas_constant() const noexcept
    {
        return specific_gas_constant_;
    }

    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;

private:
    GasModel(Real gamma, Real molar_mass, Real specific_gas_constant)
        : gamma_(gamma)
        , molar_mass_(molar_mass)
        , specific_gas_constant_(specific_gas_constant)
    {
    }

    Real gamma_ = 0.0;
    Real molar_mass_ = 0.0;
    Real specific_gas_constant_ = 0.0;
};

struct ReferenceInput {
    Real velocity = 0.0;
    Real density = 0.0;
    Real temperature = 0.0;
    Real length = 0.0;
    Real viscosity = 0.0;
    std::optional<Real> reynolds;
    std::optional<Real> mach;
};

class ReferenceScales {
public:
    [[nodiscard]] static ReferenceScales derive(
        const ReferenceInput& input,
        const GasModel& gas);

    [[nodiscard]] Real velocity() const noexcept { return velocity_; }
    [[nodiscard]] Real density() const noexcept { return density_; }
    [[nodiscard]] Real temperature() const noexcept { return temperature_; }
    [[nodiscard]] Real length() const noexcept { return length_; }
    [[nodiscard]] Real viscosity() const noexcept { return viscosity_; }
    [[nodiscard]] Real reynolds() const noexcept { return reynolds_; }
    [[nodiscard]] Real mach() const noexcept { return mach_; }
    [[nodiscard]] Real dynamic_pressure() const noexcept { return dynamic_pressure_; }
    [[nodiscard]] Real time() const noexcept { return time_; }

    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;

private:
    ReferenceScales(
        Real velocity,
        Real density,
        Real temperature,
        Real length,
        Real viscosity,
        Real reynolds,
        Real mach,
        Real dynamic_pressure,
        Real time)
        : velocity_(velocity)
        , density_(density)
        , temperature_(temperature)
        , length_(length)
        , viscosity_(viscosity)
        , reynolds_(reynolds)
        , mach_(mach)
        , dynamic_pressure_(dynamic_pressure)
        , time_(time)
    {
    }

    Real velocity_ = 0.0;
    Real density_ = 0.0;
    Real temperature_ = 0.0;
    Real length_ = 0.0;
    Real viscosity_ = 0.0;
    Real reynolds_ = 0.0;
    Real mach_ = 0.0;
    Real dynamic_pressure_ = 0.0;
    Real time_ = 0.0;
};

struct NumericalFloors {
    Real density = 1.0e-12;
    Real pressure = 1.0e-12;
    Real temperature = 1.0e-12;
    Real jacobian_absolute = 0.0;
    Real jacobian_relative = 1.0e-12;
    Real face_area_absolute = 0.0;
    Real face_area_relative = 1.0e-12;
    Real reconstruction_scale = 1.0e-12;
    Real reconstruction_epsilon = 1.0e-6;

    void validate() const;
    [[nodiscard]] Real jacobian_floor(Real reference_volume) const;
    [[nodiscard]] Real face_area_floor(Real reference_area) const;
};

struct ReconstructionScaling {
    std::array<Real, fluid_components> component {{1.0, 1.0, 1.0, 1.0, 1.0}};
    Real scale_floor = 1.0e-12;
    Real epsilon = 1.0e-6;

    void validate() const;
    [[nodiscard]] Real normalized_smoothness(Real beta, int component_index) const;
};

enum TemperaturePrimitiveIndex : int {
    temperature_density = 0,
    temperature_velocity_x = 1,
    temperature_velocity_y = 2,
    temperature_velocity_z = 3,
    temperature_value = 4,
};

using TemperaturePrimitiveState = std::array<Real, fluid_components>;
using PressurePrimitiveState = std::array<Real, fluid_components>;
using ThermodynamicConservativeState = std::array<Real, fluid_components>;

[[nodiscard]] PressurePrimitiveState pressure_primitive(
    const TemperaturePrimitiveState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

[[nodiscard]] TemperaturePrimitiveState temperature_primitive(
    const PressurePrimitiveState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

[[nodiscard]] ThermodynamicConservativeState thermodynamic_conservative(
    const TemperaturePrimitiveState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

[[nodiscard]] TemperaturePrimitiveState temperature_primitive_from_conservative(
    const ThermodynamicConservativeState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

[[nodiscard]] Real thermodynamic_sound_speed(
    const TemperaturePrimitiveState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

[[nodiscard]] Real thermodynamic_total_enthalpy(
    const TemperaturePrimitiveState& state,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors = {},
    int dimension = 3);

} // namespace wcns
