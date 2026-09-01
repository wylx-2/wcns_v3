#include "test_support.hpp"

#include <wcns/solver/transport_model.hpp>
#include <wcns/solver/viscous_flux.hpp>

#include <cmath>
#include <stdexcept>

namespace {

wcns::GasModel test_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales test_reference(const wcns::GasModel& gas)
{
    wcns::ReferenceInput input;
    input.velocity = 100.0;
    input.density = 1.0;
    input.temperature = 300.0;
    input.length = 1.0;
    input.viscosity = 1.0e-5;
    return wcns::ReferenceScales::derive(input, gas);
}

} // namespace

// 验收常黏度、Sutherland 定律及热传导系数的无量纲定义和非法配置。
void test_transport_model()
{
    using namespace wcns;
    const auto gas = test_gas();
    const auto reference = test_reference(gas);

    TransportConfig constant_config;
    constant_config.prandtl = 0.8;
    constant_config.viscosity = ConstantViscosity {2.0};
    const TransportModel constant(constant_config);
    WCNS_REQUIRE_NEAR(constant.viscosity(0.25), 2.0, 0.0);
    WCNS_REQUIRE_NEAR(
        constant.thermal_coefficient(3.0, gas, reference),
        2.0 / ((gas.gamma() - 1.0) * reference.mach()
            * reference.mach() * constant_config.prandtl),
        1.0e-13);
    WCNS_REQUIRE(constant_config.restart_signature().find("transport_v1;") == 0);

    TransportConfig sutherland_config;
    sutherland_config.viscosity = SutherlandViscosity {1.25, 0.4};
    const TransportModel sutherland(sutherland_config);
    WCNS_REQUIRE_NEAR(sutherland.viscosity(1.0), 1.25, 1.0e-15);
    WCNS_REQUIRE_NEAR(
        sutherland.viscosity(4.0),
        1.25 * 8.0 * 1.4 / 4.4,
        1.0e-14);
    WCNS_REQUIRE_THROWS(PhysicsError, sutherland.viscosity(0.0));

    TransportConfig invalid;
    invalid.prandtl = 0.0;
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid.validate());
    invalid.prandtl = 0.72;
    invalid.viscosity = SutherlandViscosity {1.0, -0.1};
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid.validate());
}

// 验收 Stokes 应力、Fourier 能量项、二维 z 约束以及本构中不含 Reynolds 缩放。
void test_viscous_cartesian_flux()
{
    using namespace wcns;
    const auto gas = test_gas();
    const auto reference = test_reference(gas);
    TransportConfig config;
    config.prandtl = 0.8;
    config.viscosity = ConstantViscosity {2.0};
    const TransportModel transport(config);

    ViscousFaceTrace trace;
    trace.state = {{1.0, 1.0, 2.0, 0.0, 2.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityX)] = {{1.0, 2.0, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityY)] = {{3.0, 4.0, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityZ)] = {{0.0, 0.0, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::Temperature)] = {{5.0, 6.0, 0.0}};

    const auto flux = compute_viscous_cartesian_flux(
        trace, transport, gas, reference, NumericalFloors {}, 2);
    const Real isotropic = (2.0 / 3.0) * 2.0 * 5.0;
    const Real tau_xx = 4.0 - isotropic;
    const Real tau_yy = 16.0 - isotropic;
    const Real tau_xy = 10.0;
    const Real chi = transport.thermal_coefficient(2.0, gas, reference);
    WCNS_REQUIRE_NEAR(flux.x[momentum_x], tau_xx, 1.0e-14);
    WCNS_REQUIRE_NEAR(flux.x[momentum_y], tau_xy, 1.0e-14);
    WCNS_REQUIRE_NEAR(flux.y[momentum_x], tau_xy, 1.0e-14);
    WCNS_REQUIRE_NEAR(flux.y[momentum_y], tau_yy, 1.0e-14);
    WCNS_REQUIRE_NEAR(
        flux.x[total_energy], tau_xx + 2.0 * tau_xy + 5.0 * chi, 1.0e-12);
    WCNS_REQUIRE_NEAR(
        flux.y[total_energy], tau_xy + 2.0 * tau_yy + 6.0 * chi, 1.0e-12);
    WCNS_REQUIRE(flux.x[density] == 0.0);
    WCNS_REQUIRE(flux.y[momentum_z] == 0.0);
    for (const Real value : flux.z) WCNS_REQUIRE(value == 0.0);

    trace.gradients[0][2] = 1.0;
    WCNS_REQUIRE_THROWS(
        PhysicsError,
        compute_viscous_cartesian_flux(
            trace, transport, gas, reference, NumericalFloors {}, 2));
}
