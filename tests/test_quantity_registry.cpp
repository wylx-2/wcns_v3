#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/runtime/quantity_registry.hpp>

#include <memory>

namespace {

wcns::StructuredBlock cartesian_block()
{
    wcns::StructuredBlock block(0, "quantity", 0, 2, 2, {5, 5, 1}, 3);
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<wcns::Real>(i) / 4.0;
            block.coordinates.y(i, j, 0) = static_cast<wcns::Real>(j) / 4.0;
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    return block;
}

} // namespace

// 验收内建物理量只遍历真实单元，并按参考量一致地切换无量纲/有量纲输出。
void test_quantity_registry()
{
    auto block = cartesian_block();
    const auto profile = wcns::ProfileFactory::create(
        wcns::AlgorithmProfileKind::PhengleiWcns);
    auto metric = wcns::initialize_metric_field(block, profile).metric;
    wcns::GasModelInput gas_input;
    gas_input.molar_mass = 0.029;
    const auto gas = wcns::GasModel::from_input(gas_input);
    wcns::ReferenceInput reference_input;
    reference_input.velocity = 10.0;
    reference_input.density = 2.0;
    reference_input.temperature = 300.0;
    reference_input.length = 4.0;
    reference_input.viscosity = 1.0e-5;
    const auto reference = wcns::ReferenceScales::derive(reference_input, gas);
    wcns::InitialConditionConfig initial;
    initial.type = "uniform";
    initial.parameters = {
        {"rho", 1.25}, {"u", 0.2}, {"v", -0.1}, {"temperature", 1.1},
    };
    wcns::FlowInitializer::initialize_block(
        block, metric, initial, gas, reference);

    auto registry = wcns::FieldQuantityRegistry::create_builtin();
    wcns::QuantityContext context {
        gas,
        reference,
        {},
        wcns::TransportModel(wcns::TransportConfig {}),
        false,
    };
    const auto rho = registry.evaluate("rho", block, metric, context);
    const auto velocity = registry.evaluate("u", block, metric, context);
    const auto jacobian = registry.evaluate("jacobian", block, metric, context);
    WCNS_REQUIRE(rho.values.size() == 16);
    WCNS_REQUIRE_NEAR(rho.values.front(), 1.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(velocity.values.front(), 0.2, 1.0e-14);
    WCNS_REQUIRE_NEAR(jacobian.values.front(), 0.0625, 1.0e-13);

    context.dimensional = true;
    const auto dimensional_rho = registry.evaluate("rho", block, metric, context);
    const auto dimensional_u = registry.evaluate("u", block, metric, context);
    const auto dimensional_j = registry.evaluate(
        "jacobian", block, metric, context);
    WCNS_REQUIRE_NEAR(dimensional_rho.values.front(), 2.5, 1.0e-13);
    WCNS_REQUIRE_NEAR(dimensional_u.values.front(), 2.0, 1.0e-13);
    WCNS_REQUIRE_NEAR(dimensional_j.values.front(), 1.0, 1.0e-12);

    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        registry.validate_selection({"rho", "missing"}));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        registry.validate_selection({"rho", "rho"}));
}
