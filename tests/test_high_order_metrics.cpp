#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>

#include <stdexcept>

namespace {

void set_affine_2d(wcns::StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = 2.0 * i + 0.5 * j + 1.0;
            block.coordinates.y(i, j, 0) = 0.25 * i + 1.5 * j - 2.0;
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
}

} // namespace

// 验收 PH 加密交错度量对二维仿射映射逐点精确且只暴露真实几何域。
void test_phenglei_high_order_metrics()
{
    using namespace wcns;
    StructuredBlock block(0, "ph-affine", 0, 2, 2, {7, 6, 1}, 3);
    set_affine_2d(block);
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto result = initialize_metric_field(block, profile);
    const Real exact_jacobian = 2.0 * 1.5 - 0.5 * 0.25;
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            WCNS_REQUIRE_NEAR(result.metric.jacobian()(i, j, 0), exact_jacobian, 2.0e-13);
            WCNS_REQUIRE_NEAR(
                result.metric.cell_coordinates().x(i, j, 0),
                2.0 * (i + 0.5) + 0.5 * (j + 0.5) + 1.0,
                2.0e-14);
            WCNS_REQUIRE_NEAR(
                result.diagnostics.reference_volume(i, j, 0),
                exact_jacobian,
                2.0e-14);
        }
    }
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i <= cells.ni; ++i) {
            WCNS_REQUIRE_NEAR(result.metric.i_faces().x(i, j, 0), 1.5, 2.0e-14);
            WCNS_REQUIRE_NEAR(result.metric.i_faces().y(i, j, 0), -0.5, 2.0e-14);
        }
    }
    for (int j = 0; j <= cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            WCNS_REQUIRE_NEAR(result.metric.j_faces().x(i, j, 0), -0.25, 2.0e-14);
            WCNS_REQUIRE_NEAR(result.metric.j_faces().y(i, j, 0), 2.0, 2.0e-14);
        }
    }
    WCNS_REQUIRE(result.diagnostics.fallback_cell_count == 0);
    WCNS_REQUIRE_THROWS(std::logic_error, result.metric.k_faces());
    WCNS_REQUIRE_THROWS(std::out_of_range, result.metric.jacobian()(-1, 0, 0));

    StructuredBlock too_small(1, "ph-small", 0, 2, 2, {2, 3, 1}, 3);
    WCNS_REQUIRE_THROWS(
        GeometryError,
        initialize_metric_field(too_small, profile));
}

// 验收 PH 三维对称守恒公式对一般仿射映射给出解析余因子和 Jacobian。
void test_phenglei_high_order_metrics_3d()
{
    using namespace wcns;
    StructuredBlock block(2, "ph-affine-3d", 0, 3, 3, {4, 4, 4}, 3);
    const auto vertices = block.vertex_extent();
    for (int k = 0; k < vertices.nk; ++k) {
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                block.coordinates.x(i, j, k) = 2.0 * i - j + k;
                block.coordinates.y(i, j, k) = i + 3.0 * j;
                block.coordinates.z(i, j, k) = j + 2.0 * k;
            }
        }
    }
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto result = initialize_metric_field(block, profile);
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                WCNS_REQUIRE_NEAR(result.metric.jacobian()(i, j, k), 15.0, 2.0e-13);
            }
        }
    }
    WCNS_REQUIRE_NEAR(result.metric.i_faces().x(1, 1, 1), 6.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().y(1, 1, 1), 3.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().z(1, 1, 1), -3.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().x(1, 1, 1), -2.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().y(1, 1, 1), 4.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().z(1, 1, 1), 1.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.k_faces().x(1, 1, 1), 1.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.k_faces().y(1, 1, 1), -2.0, 2.0e-14);
    WCNS_REQUIRE_NEAR(result.metric.k_faces().z(1, 1, 1), 7.0, 2.0e-14);
}

// 验收 SCMM6 共同中心度量在二维仿射映射上保持统一 delta 和解析几何量。
void test_scmm6_high_order_metrics()
{
    using namespace wcns;
    StructuredBlock block(3, "scmm-affine", 0, 2, 2, {9, 8, 1}, 3);
    set_affine_2d(block);
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    const auto result = initialize_metric_field(block, profile);
    const Real exact_jacobian = 2.0 * 1.5 - 0.5 * 0.25;
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            WCNS_REQUIRE_NEAR(result.metric.jacobian()(i, j, 0), exact_jacobian, 3.0e-11);
            WCNS_REQUIRE_NEAR(
                result.metric.cell_coordinates().y(i, j, 0),
                0.25 * (i + 0.5) + 1.5 * (j + 0.5) - 2.0,
                3.0e-12);
        }
    }
    WCNS_REQUIRE_NEAR(result.metric.i_faces().x(4, 3, 0), 1.5, 3.0e-11);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().y(4, 3, 0), -0.5, 3.0e-11);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().x(4, 3, 0), -0.25, 3.0e-11);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().y(4, 3, 0), 2.0, 3.0e-11);

    MetricBuildOptions invalid_fallback;
    invalid_fallback.fallback = MetricFallback::PhengleiFiniteVolume;
    WCNS_REQUIRE_THROWS(
        ProfileError,
        initialize_metric_field(block, profile, invalid_fallback));
    StructuredBlock too_small(4, "scmm-small", 0, 2, 2, {5, 6, 1}, 3);
    WCNS_REQUIRE_THROWS(
        GeometryError,
        initialize_metric_field(too_small, profile));
}

// 验收 SCMM6 三维嵌套对称导数在一般仿射映射上给出解析余因子。
void test_scmm6_high_order_metrics_3d()
{
    using namespace wcns;
    StructuredBlock block(5, "scmm-affine-3d", 0, 3, 3, {7, 7, 7}, 3);
    const auto vertices = block.vertex_extent();
    for (int k = 0; k < vertices.nk; ++k) {
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                block.coordinates.x(i, j, k) = 2.0 * i - j + k;
                block.coordinates.y(i, j, k) = i + 3.0 * j;
                block.coordinates.z(i, j, k) = j + 2.0 * k;
            }
        }
    }
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    const auto result = initialize_metric_field(block, profile);
    WCNS_REQUIRE_NEAR(result.metric.jacobian()(2, 2, 2), 15.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().x(2, 2, 2), 6.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().y(2, 2, 2), 3.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.i_faces().z(2, 2, 2), -3.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().x(2, 2, 2), -2.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.j_faces().y(2, 2, 2), 4.0, 5.0e-10);
    WCNS_REQUIRE_NEAR(result.metric.k_faces().z(2, 2, 2), 7.0, 5.0e-10);
}
