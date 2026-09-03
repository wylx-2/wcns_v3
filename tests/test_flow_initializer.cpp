#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/runtime/flow_initializer.hpp>

#include <cmath>

namespace {

wcns::GasModel initializer_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales initializer_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive(
        {340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}},
        gas);
}

wcns::StructuredBlock initializer_block()
{
    wcns::StructuredBlock block(0, "initial", 0, 2, 2, {9, 9, 1}, 3);
    for (int j = 0; j < 9; ++j) {
        for (int i = 0; i < 9; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<wcns::Real>(i) / 8.0;
            block.coordinates.y(i, j, 0) = static_cast<wcns::Real>(j) / 8.0;
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    return block;
}

} // namespace

void test_flow_initializer()
{
    using namespace wcns;
    const auto gas = initializer_gas();
    const auto reference = initializer_reference(gas);
    const NumericalFloors floors;
    const auto profile = ProfileFactory::create(
        AlgorithmProfileKind::PhengleiWcns);

    {
        auto block = initializer_block();
        const auto metrics = initialize_metric_field(block, profile).metric;
        InitialConditionConfig config;
        config.type = "uniform";
        config.parameters = {
            {"rho", 1.1}, {"u", 0.2}, {"v", -0.1},
            {"w", 0.0}, {"temperature", 0.9},
        };
        FlowInitializer::initialize_block(
            block, metrics, config, gas, reference, floors);
        const auto expected = thermodynamic_conservative(
            {1.1, 0.2, -0.1, 0.0, 0.9},
            gas,
            reference,
            floors,
            2);
        const auto actual = load_conservative(
            block.flow.conservative, {3, 4, 0});
        for (std::size_t component = 0; component < expected.size(); ++component) {
            WCNS_REQUIRE_NEAR(actual[component], expected[component], 1.0e-14);
        }
        WCNS_REQUIRE(std::isnan(block.flow.conservative(-1, 0, 0, 0)));
    }
    {
        auto block = initializer_block();
        const auto metrics = initialize_metric_field(block, profile).metric;
        InitialConditionConfig config;
        config.type = "sod_x";
        config.parameters = {{"x0", 0.5}};
        FlowInitializer::initialize_block(
            block, metrics, config, gas, reference, floors);
        const auto left = temperature_primitive_from_conservative(
            load_conservative(block.flow.conservative, {0, 0, 0}),
            gas, reference, floors, 2);
        const auto right = temperature_primitive_from_conservative(
            load_conservative(block.flow.conservative, {7, 0, 0}),
            gas, reference, floors, 2);
        WCNS_REQUIRE_NEAR(left[0], 1.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(right[0], 0.125, 1.0e-14);
        WCNS_REQUIRE_NEAR(
            pressure_primitive(left, gas, reference, floors, 2)[4],
            1.0,
            1.0e-14);
        WCNS_REQUIRE_NEAR(
            pressure_primitive(right, gas, reference, floors, 2)[4],
            0.1,
            1.0e-14);
    }
    {
        InitialConditionConfig config;
        config.type = "quadrant_riemann";
        const auto ne = FlowInitializer::evaluate(
            config, {0.75, 0.75, 0.0},
            gas, reference, floors, 2);
        const auto sw = FlowInitializer::evaluate(
            config, {0.25, 0.25, 0.0},
            gas, reference, floors, 2);
        WCNS_REQUIRE_NEAR(ne[0], 1.5, 1.0e-14);
        WCNS_REQUIRE_NEAR(sw[0], 0.138, 1.0e-14);
    }
    {
        InitialConditionConfig config;
        config.type = "isentropic_vortex";
        config.parameters = {{"x0", 5.0}, {"y0", 5.0}, {"beta", 5.0}};
        const auto center = FlowInitializer::evaluate(
            config, {5.0, 5.0, 0.0},
            gas, reference, floors, 2);
        constexpr Real pi = 3.141592653589793238462643383279502884;
        const Real expected_temperature = 1.0
            - (gas.gamma() - 1.0) * reference.mach() * reference.mach()
                * 25.0 * std::exp(1.0) / (8.0 * pi * pi);
        WCNS_REQUIRE_NEAR(center[4], expected_temperature, 1.0e-14);
        WCNS_REQUIRE_NEAR(
            center[0],
            std::pow(expected_temperature, 1.0 / (gas.gamma() - 1.0)),
            1.0e-14);
    }
    {
        InitialConditionConfig config;
        config.type = "couette";
        config.parameters = {
            {"y0", 0.0}, {"y1", 2.0},
            {"lower_velocity", -0.25}, {"upper_velocity", 0.75},
            {"lower_temperature", 1.0}, {"upper_temperature", 1.4},
            {"temperature_curvature", 0.2}, {"velocity_curvature", 0.4},
            {"pressure", 0.8},
        };
        const auto middle = FlowInitializer::evaluate(
            config, {0.0, 1.0, 0.0}, gas, reference, floors, 2);
        WCNS_REQUIRE_NEAR(middle[1], 0.35, 1.0e-14);
        WCNS_REQUIRE_NEAR(middle[4], 1.25, 1.0e-14);
        WCNS_REQUIRE_NEAR(
            pressure_primitive(middle, gas, reference, floors, 2)[4],
            0.8,
            1.0e-14);
    }
    {
        InitialConditionConfig config;
        config.type = "linear_conduction";
        config.parameters = {
            {"y0", -1.0}, {"y1", 1.0},
            {"lower_temperature", 1.0}, {"upper_temperature", 2.0},
            {"pressure", 0.6},
        };
        const auto middle = FlowInitializer::evaluate(
            config, {0.0, 0.0, 0.0}, gas, reference, floors, 2);
        WCNS_REQUIRE_NEAR(middle[1], 0.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(middle[4], 1.5, 1.0e-14);
        WCNS_REQUIRE_NEAR(
            pressure_primitive(middle, gas, reference, floors, 2)[4],
            0.6,
            1.0e-14);
    }
}
