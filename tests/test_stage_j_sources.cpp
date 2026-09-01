#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/solver/source_operator.hpp>
#include <wcns/solver/time_integrator.hpp>

#include <vector>

namespace {

wcns::StructuredBlock make_source_block()
{
    using namespace wcns;
    StructuredBlock block(0, "source", 0, 2, 2, {7, 7, 1}, 3);
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    const ConservativeState state {2.0, 4.0, 2.0, 0.0, 10.0};
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            store_state(block.flow.conservative, {i, j, 0}, state);
        }
    }
    return block;
}

} // namespace

// 验收常量、体力和制造源按配置顺序逐点累加，且二维 z 动量与体力能量功正确。
void test_stage_j_source_models()
{
    using namespace wcns;
    SourceTermConfig config;
    config.enable_source_terms = true;
    config.models = {
        SourceModelKind::UniformConservative,
        SourceModelKind::BodyForce,
        SourceModelKind::ManufacturedSolution,
    };
    config.uniform_conservative = {{1.0, 2.0, 3.0, 0.0, 4.0}};
    config.body_acceleration = {{0.5, -0.25, 0.0}};
    config.manufactured_amplitude = {{0.1, 0.2, 0.3, 0.0, 0.4}};
    const auto registry = SourceTermRegistry::create_stage_j(config);
    WCNS_REQUIRE(registry.size() == 3);
    const ConservativeState state {2.0, 4.0, 2.0, 0.0, 10.0};
    const auto source = registry.evaluate(state, {{2.0, 3.0, 0.0}}, 0.5, 2);
    const Real shape = 6.5;
    WCNS_REQUIRE_NEAR(source[0], 1.0 + 0.1 * shape, 1.0e-14);
    WCNS_REQUIRE_NEAR(source[1], 2.0 + 1.0 + 0.2 * shape, 1.0e-14);
    WCNS_REQUIRE_NEAR(source[2], 3.0 - 0.5 + 0.3 * shape, 1.0e-14);
    WCNS_REQUIRE(source[3] == 0.0);
    WCNS_REQUIRE_NEAR(source[4], 4.0 + 1.5 + 0.4 * shape, 1.0e-14);
    WCNS_REQUIRE(config.restart_signature().find("source_terms_v2;") == 0);

    auto invalid = config;
    invalid.body_acceleration[2] = 1.0;
    const auto invalid_registry = SourceTermRegistry::create_stage_j(invalid);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        invalid_registry.evaluate(state, {{0.0, 0.0, 0.0}}, 0.0, 2));
}

// 验收源项只累加到真实单元 dU/dt、关闭路径逐位不变且全局积分仅在诊断时乘 J。
void test_source_operator_balance()
{
    using namespace wcns;
    auto block = make_source_block();
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto metric = initialize_metric_field(block, profile).metric;
    block.flow.residual.fill(7.0);
    const auto before = block.flow.residual;
    const auto disabled = SourceTermRegistry::create_stage_j({});
    add_source_terms(block, metric, disabled, 0.0);
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            for (int component = 0; component < euler_components; ++component) {
                WCNS_REQUIRE(block.flow.residual(i, j, 0, component)
                    == before(i, j, 0, component));
            }
        }
    }

    SourceTermConfig config;
    config.enable_source_terms = true;
    config.models = {SourceModelKind::UniformConservative};
    config.uniform_conservative = {{1.0, 2.0, 3.0, 0.0, 5.0}};
    const auto registry = SourceTermRegistry::create_stage_j(config);
    block.flow.residual.fill(0.0);
    add_source_terms(block, metric, registry, 0.25);
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            WCNS_REQUIRE(block.flow.residual(i, j, 0, 0) == 1.0);
            WCNS_REQUIRE(block.flow.residual(i, j, 0, 4) == 5.0);
        }
    }
    const auto integral = volume_weighted_source(block, metric, registry, 0.25);
    WCNS_REQUIRE_NEAR(integral[0], static_cast<Real>(cells.size()), 2.0e-12);
    WCNS_REQUIRE_NEAR(integral[4], 5.0 * static_cast<Real>(cells.size()), 1.0e-11);
}

// 验收 SSPRK3 在三个残差调用中传递 t、t+dt、t+dt/2 并保持常量源三阶更新。
void test_timed_ssprk3_sources()
{
    using namespace wcns;
    StructuredBlock block(0, "rk-source", 0, 2, 2, {2, 2, 1}, 0);
    block.flow.conservative.fill(1.0);
    std::vector<Real> times;
    const std::vector<StructuredBlock*> blocks {&block};
    advance_ssprk3(blocks, 0.1, 2.0, [&](Real stage_time) {
        times.push_back(stage_time);
        block.flow.residual.fill(2.0);
    });
    WCNS_REQUIRE(times.size() == 3);
    WCNS_REQUIRE_NEAR(times[0], 2.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(times[1], 2.1, 1.0e-15);
    WCNS_REQUIRE_NEAR(times[2], 2.05, 1.0e-15);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            block.flow.conservative(0, 0, 0, component), 1.2, 2.0e-15);
    }
}
