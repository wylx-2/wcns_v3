#include "test_support.hpp"

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/solver/viscous_gradient.hpp>
#include <wcns/solver/viscous_halo.hpp>

#include <cmath>

namespace {

void fill_identity_coordinates(wcns::StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<wcns::Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<wcns::Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
}

void fill_linear_temperature_state(wcns::StructuredBlock& block)
{
    const auto cells = block.cell_extent();
    auto& state = block.flow.temperature_primitive;
    const int ghost = state.ghost_width();
    const auto store = [&](int i, int j) {
        const wcns::Real x = static_cast<wcns::Real>(i) + 0.5;
        const wcns::Real y = static_cast<wcns::Real>(j) + 0.5;
        state(i, j, 0, wcns::temperature_density) = 1.0;
        state(i, j, 0, wcns::temperature_velocity_x) = 1.0 + x + 2.0 * y;
        state(i, j, 0, wcns::temperature_velocity_y) = -2.0 - 0.5 * x + 0.25 * y;
        state(i, j, 0, wcns::temperature_velocity_z) = 0.0;
        state(i, j, 0, wcns::temperature_value) = 3.0 + 0.1 * x + 0.2 * y;
    };
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) store(i, j);
    }
    for (int layer = 1; layer <= ghost; ++layer) {
        for (int j = 0; j < cells.nj; ++j) {
            store(-layer, j);
            store(cells.ni - 1 + layer, j);
        }
        for (int i = 0; i < cells.ni; ++i) {
            store(i, -layer);
            store(i, cells.nj - 1 + layer);
        }
    }
}

void run_linear_gradient(wcns::AlgorithmProfileKind kind)
{
    using namespace wcns;
    StructuredBlock block(0, "linear-gradient", 0, 2, 2, {9, 9, 1}, 3);
    fill_identity_coordinates(block);
    fill_linear_temperature_state(block);
    const auto profile = ProfileFactory::create(kind);
    const auto metric = initialize_metric_field(block, profile).metric;
    const auto operands = compute_gradient_face_operands(block, metric, profile, 7);
    const auto gradients = compute_primitive_gradients(
        block, metric, operands, profile);
    WCNS_REQUIRE(gradients.version() == 7);
    WCNS_REQUIRE(gradients.halo_layers()
        == (kind == AlgorithmProfileKind::PhengleiWcns ? 2 : 3));
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const Index3 cell {i, j, 0};
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::VelocityX, 0), 1.0, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::VelocityX, 1), 2.0, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::VelocityY, 0), -0.5, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::VelocityY, 1), 0.25, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::Temperature, 0), 0.1, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                gradients(cell, ViscousPrimitive::Temperature, 1), 0.2, 2.0e-13);
            for (int variable = 0; variable < viscous_primitive_components; ++variable) {
                WCNS_REQUIRE(gradients(
                    cell, static_cast<ViscousPrimitive>(variable), 2) == 0.0);
            }
            for (int direction = 0; direction < 3; ++direction) {
                WCNS_REQUIRE(gradients(
                    cell, ViscousPrimitive::VelocityZ, direction) == 0.0);
            }
        }
    }
    WCNS_REQUIRE(std::isnan(gradients.values()(-1, -1, 0, 0)));
}

} // namespace

// 验收两套 profile 在二维仿射网格上以守恒形式精确恢复线性 primitive 梯度。
void test_viscous_linear_gradients()
{
    run_linear_gradient(wcns::AlgorithmProfileKind::PhengleiWcns);
    run_linear_gradient(wcns::AlgorithmProfileKind::Scmm6Wcns);
}

// 验收旋转周期的速度梯度二阶张量和温度梯度矢量接收侧变换。
void test_viscous_gradient_periodic_transform()
{
    using namespace wcns;
    GradientExchangeDescriptor descriptor;
    descriptor.connection = 0;
    descriptor.receiver_block = 0;
    descriptor.donor_block = 1;
    descriptor.dimension = 3;
    descriptor.periodic.rotation = {{{{0.0, -1.0, 0.0}},
        {{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
    PrimitiveGradients donor {{
        {{1.0, 2.0, 3.0}},
        {{4.0, 5.0, 6.0}},
        {{7.0, 8.0, 9.0}},
        {{10.0, 11.0, 12.0}},
    }};
    const auto received = transform_primitive_gradients_for_receiver(
        donor, descriptor);
    WCNS_REQUIRE_NEAR(received[0][0], 5.0, 0.0);
    WCNS_REQUIRE_NEAR(received[0][1], -4.0, 0.0);
    WCNS_REQUIRE_NEAR(received[1][0], -2.0, 0.0);
    WCNS_REQUIRE_NEAR(received[1][1], 1.0, 0.0);
    WCNS_REQUIRE_NEAR(received[3][0], 11.0, 0.0);
    WCNS_REQUIRE_NEAR(received[3][1], -10.0, 0.0);

    FaceFluxExchangeDescriptor face;
    face.connection = 0;
    face.receiver_block = 0;
    face.donor_block = 1;
    face.orientation = -1.0;
    face.periodic = descriptor.periodic;
    GradientOperandState operand {};
    for (int i = 0; i < gradient_operand_components; ++i) {
        operand[static_cast<std::size_t>(i)] = static_cast<Real>(i + 1);
    }
    const auto transformed = transform_gradient_operand_for_receiver(
        operand, face, 3);
    WCNS_REQUIRE_NEAR(transformed[0], -5.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[1], 4.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[3], 2.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[4], -1.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[9], -11.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[10], 10.0, 0.0);
}
