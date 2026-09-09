#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/runtime/quantity_registry.hpp>

#include <memory>
#include <string>
#include <utility>
#include <vector>

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

class ConstantFieldQuantity final : public wcns::IFieldQuantity {
public:
    ConstantFieldQuantity(
        std::string name,
        std::vector<std::string> dependencies,
        wcns::Real value)
        : value_(value)
    {
        descriptor_.name = std::move(name);
        descriptor_.dependencies = std::move(dependencies);
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate_cell(
        const wcns::StructuredBlock&,
        const wcns::MetricField&,
        wcns::Index3,
        const wcns::QuantityContext&) const override
    {
        return value_;
    }

private:
    wcns::QuantityDescriptor descriptor_;
    wcns::Real value_ = 0.0;
};

class ConstantStatisticQuantity final : public wcns::IStatisticQuantity {
public:
    explicit ConstantStatisticQuantity(std::string name)
    {
        descriptor_.name = std::move(name);
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate(const wcns::StatisticContext&) const override
    {
        return 7.0;
    }

private:
    wcns::QuantityDescriptor descriptor_;
};

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

    registry.register_quantity(std::make_shared<ConstantFieldQuantity>(
        "custom_cell", std::vector<std::string> {"rho"}, 7.0));
    registry.validate_selection({"custom_cell"});
    const auto custom = registry.evaluate(
        "custom_cell", block, metric, context);
    WCNS_REQUIRE(custom.values.size() == 16);
    WCNS_REQUIRE_NEAR(custom.values.front(), 7.0, 1.0e-15);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        registry.register_quantity(std::make_shared<ConstantFieldQuantity>(
            "custom_cell", std::vector<std::string> {}, 8.0)));

    wcns::FieldQuantityRegistry unknown_dependency;
    unknown_dependency.register_quantity(std::make_shared<ConstantFieldQuantity>(
        "unknown_root", std::vector<std::string> {"missing"}, 1.0));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        unknown_dependency.validate_selection({"unknown_root"}));

    wcns::FieldQuantityRegistry cyclic;
    cyclic.register_quantity(std::make_shared<ConstantFieldQuantity>(
        "cycle_a", std::vector<std::string> {"cycle_b"}, 1.0));
    cyclic.register_quantity(std::make_shared<ConstantFieldQuantity>(
        "cycle_b", std::vector<std::string> {"cycle_a"}, 1.0));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        cyclic.validate_selection({"cycle_a"}));

    auto statistics = wcns::StatisticRegistry::create_builtin();
    wcns::register_xz_plane_statistics(statistics, {0, 3});
    statistics.validate_selection({
        "xz_mean_u_j0", "xz_mass_flow_x_j0",
        "xz_mean_u_j3", "xz_mass_flow_x_j3"});
    WCNS_REQUIRE(wcns::xz_plane_statistic_names({2})
        == std::vector<std::string>({"xz_mean_u_j2", "xz_mass_flow_x_j2"}));
    wcns::PartitionConfig partition_config;
    partition_config.min_cells_per_active_direction = 1;
    const auto partition = wcns::StructuredPartitionPlan::build(
        {{0, "three-dimensional", 3, {4, 4, 4}}}, 1, partition_config);
    wcns::validate_xz_plane_statistics({0, 3}, partition);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        wcns::validate_xz_plane_statistics({4}, partition));
    const auto two_dimensional = wcns::StructuredPartitionPlan::build(
        {{0, "two-dimensional", 2, {4, 4, 1}}}, 1, partition_config);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        wcns::validate_xz_plane_statistics({0}, two_dimensional));
    wcns::register_yz_plane_statistics(statistics, {0.0, 3.141592653589793});
    const auto yz_names = wcns::yz_plane_statistic_names(2);
    WCNS_REQUIRE(yz_names == std::vector<std::string>({
        "yz_mean_u_plane0", "yz_mass_flow_x_plane0",
        "yz_mean_u_plane1", "yz_mass_flow_x_plane1"}));
    statistics.validate_selection(yz_names);
    wcns::validate_yz_plane_statistics({0.0, 3.141592653589793}, partition);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        wcns::validate_yz_plane_statistics({0.0}, two_dimensional));
    wcns::register_channel_wall_statistics(
        statistics, "bottom", "top", 1.0);
    statistics.validate_selection(wcns::channel_wall_statistic_names());
    WCNS_REQUIRE(wcns::channel_wall_statistic_names()
        == std::vector<std::string>({
            "channel_wall_shear_lower", "channel_wall_shear_upper",
            "channel_wall_shear_mean", "channel_friction_velocity",
            "channel_re_tau"}));
    statistics.register_quantity(
        std::make_shared<ConstantStatisticQuantity>("custom_statistic"));
    statistics.validate_selection({"custom_statistic"});
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        statistics.register_quantity(
            std::make_shared<ConstantStatisticQuantity>("custom_statistic")));
}
