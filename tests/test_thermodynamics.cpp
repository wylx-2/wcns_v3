#include "test_support.hpp"

#include <wcns/physics/thermodynamics.hpp>

#include <cmath>
#include <limits>
#include <string>

void test_thermodynamics()
{
    using namespace wcns;

    const auto gas = GasModel::from_input(GasModelInput {1.4, {}, 287.0});
    WCNS_REQUIRE_NEAR(gas.gamma(), 1.4, 0.0);
    WCNS_REQUIRE_NEAR(gas.specific_gas_constant(), 287.0, 0.0);
    WCNS_REQUIRE_NEAR(gas.molar_mass(), universal_gas_constant / 287.0, 1.0e-16);
    WCNS_REQUIRE(gas.summary().find("specific_gas_constant=287") != std::string::npos);
    WCNS_REQUIRE(gas.restart_signature() == gas.summary());

    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        GasModel::from_input(GasModelInput {1.4, {}, {}}));
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        GasModel::from_input(GasModelInput {1.4, 0.029, 287.0}));
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        GasModel::from_input(GasModelInput {1.0, 0.029, {}}));

    const ReferenceInput input {340.0, 1.2, 300.0, 1.5, 1.8e-5, {}, {}};
    const auto reference = ReferenceScales::derive(input, gas);
    WCNS_REQUIRE_NEAR(reference.reynolds(), 1.2 * 340.0 * 1.5 / 1.8e-5, 1.0e-8);
    WCNS_REQUIRE_NEAR(
        reference.mach(),
        340.0 / std::sqrt(1.4 * 287.0 * 300.0),
        1.0e-15);
    WCNS_REQUIRE_NEAR(reference.dynamic_pressure(), 1.2 * 340.0 * 340.0, 1.0e-12);
    WCNS_REQUIRE_NEAR(reference.time(), 1.5 / 340.0, 1.0e-18);
    WCNS_REQUIRE(reference.summary().find("Re=") != std::string::npos);
    WCNS_REQUIRE(reference.summary().find("Ma=") != std::string::npos);
    WCNS_REQUIRE(reference.restart_signature() == reference.summary());

    auto forbidden_reynolds = input;
    forbidden_reynolds.reynolds = reference.reynolds();
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        ReferenceScales::derive(forbidden_reynolds, gas));
    auto forbidden_mach = input;
    forbidden_mach.mach = reference.mach();
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        ReferenceScales::derive(forbidden_mach, gas));
    auto invalid_reference = input;
    invalid_reference.temperature = 0.0;
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        ReferenceScales::derive(invalid_reference, gas));

    const NumericalFloors floors;
    floors.validate();
    WCNS_REQUIRE_NEAR(floors.jacobian_floor(2.0), 2.0e-12, 0.0);
    WCNS_REQUIRE_NEAR(floors.face_area_floor(3.0), 3.0e-12, 0.0);
    auto invalid_floors = floors;
    invalid_floors.temperature = 0.0;
    WCNS_REQUIRE_THROWS(PhysicsConfigurationError, invalid_floors.validate());
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        floors.jacobian_floor(0.0));

    ReconstructionScaling scaling;
    scaling.component = {{2.0, 4.0, 1.0, 1.0, 8.0}};
    scaling.validate();
    WCNS_REQUIRE_NEAR(scaling.normalized_smoothness(16.0, 1), 1.0, 0.0);
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        scaling.normalized_smoothness(-1.0, 0));
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        scaling.normalized_smoothness(1.0, fluid_components));

    const TemperaturePrimitiveState temperature_state {{1.2, 0.3, -0.2, 0.0, 1.1}};
    const auto pressure_state = pressure_primitive(
        temperature_state, gas, reference, floors, 2);
    WCNS_REQUIRE_NEAR(
        pressure_state[4],
        1.2 * 1.1 / (1.4 * reference.mach() * reference.mach()),
        1.0e-14);
    const auto recovered_temperature = temperature_primitive(
        pressure_state, gas, reference, floors, 2);
    for (int component = 0; component < fluid_components; ++component) {
        WCNS_REQUIRE_NEAR(
            recovered_temperature[static_cast<std::size_t>(component)],
            temperature_state[static_cast<std::size_t>(component)],
            1.0e-14);
    }

    const auto conservative = thermodynamic_conservative(
        temperature_state, gas, reference, floors, 2);
    const auto recovered_conservative = temperature_primitive_from_conservative(
        conservative, gas, reference, floors, 2);
    for (int component = 0; component < fluid_components; ++component) {
        WCNS_REQUIRE_NEAR(
            recovered_conservative[static_cast<std::size_t>(component)],
            temperature_state[static_cast<std::size_t>(component)],
            1.0e-14);
    }
    WCNS_REQUIRE(thermodynamic_sound_speed(
        temperature_state, gas, reference, floors, 2) > 0.0);
    WCNS_REQUIRE(thermodynamic_total_enthalpy(
        temperature_state, gas, reference, floors, 2) > 0.0);

    auto invalid_2d = temperature_state;
    invalid_2d[temperature_velocity_z] = 1.0e-6;
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        thermodynamic_conservative(invalid_2d, gas, reference, floors, 2));
    auto non_finite = temperature_state;
    non_finite[temperature_value] = std::numeric_limits<Real>::quiet_NaN();
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        thermodynamic_conservative(non_finite, gas, reference, floors, 3));
}
