#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

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

struct WarpedMetricError {
    wcns::Real jacobian_l2 = 0.0;
    wcns::Real identity_max = 0.0;
};

WarpedMetricError warped_metric_error(
    wcns::AlgorithmProfileKind kind,
    int cell_count)
{
    using namespace wcns;
    StructuredBlock block(
        10 + cell_count,
        "warped",
        0,
        2,
        2,
        {cell_count + 1, cell_count + 1, 1},
        3);
    constexpr Real amplitude = 0.01;
    constexpr Real two_pi = 6.283185307179586476925286766559;
    for (int j = 0; j <= cell_count; ++j) {
        for (int i = 0; i <= cell_count; ++i) {
            const Real s = static_cast<Real>(i) / cell_count;
            const Real t = static_cast<Real>(j) / cell_count;
            block.coordinates.x(i, j, 0)
                = s + amplitude * std::sin(two_pi * s) * std::sin(two_pi * t);
            block.coordinates.y(i, j, 0)
                = t + amplitude * std::cos(two_pi * s) * std::sin(two_pi * t);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    const auto profile = ProfileFactory::create(kind);
    MetricBuildOptions options;
    options.maximum_reference_relative_difference = 1.0;
    const auto result = initialize_metric_field(block, profile, options);
    const auto line = LineOperators::build(profile, cell_count);
    Real squared_error = 0.0;
    std::size_t sample_count = 0;
    const int margin = kind == AlgorithmProfileKind::PhengleiWcns ? 3 : 4;
    for (int j = margin; j < cell_count - margin; ++j) {
        for (int i = margin; i < cell_count - margin; ++i) {
            const Real s = (static_cast<Real>(i) + 0.5) / cell_count;
            const Real t = (static_cast<Real>(j) + 0.5) / cell_count;
            const Real x_s = 1.0
                + amplitude * two_pi * std::cos(two_pi * s) * std::sin(two_pi * t);
            const Real x_t
                = amplitude * two_pi * std::sin(two_pi * s) * std::cos(two_pi * t);
            const Real y_s
                = -amplitude * two_pi * std::sin(two_pi * s) * std::sin(two_pi * t);
            const Real y_t = 1.0
                + amplitude * two_pi * std::cos(two_pi * s) * std::cos(two_pi * t);
            const Real exact = (x_s * y_t - x_t * y_s)
                / static_cast<Real>(cell_count * cell_count);
            const Real error = result.metric.jacobian()(i, j, 0) - exact;
            squared_error += error * error;
            ++sample_count;
        }
    }

    Real identity_max = 0.0;
    for (int component = 0; component < 3; ++component) {
        const auto& i_component = component == 0 ? result.metric.i_faces().x
            : component == 1                      ? result.metric.i_faces().y
                                                  : result.metric.i_faces().z;
        const auto& j_component = component == 0 ? result.metric.j_faces().x
            : component == 1                      ? result.metric.j_faces().y
                                                  : result.metric.j_faces().z;
        std::vector<std::vector<Real>> di(
            static_cast<std::size_t>(cell_count),
            std::vector<Real>(static_cast<std::size_t>(cell_count)));
        for (int j = 0; j < cell_count; ++j) {
            std::vector<Real> faces(static_cast<std::size_t>(cell_count + 1));
            for (int i = 0; i <= cell_count; ++i) {
                faces[static_cast<std::size_t>(i)] = i_component(i, j, 0);
            }
            const auto derivative = line.differentiate(faces);
            for (int i = 0; i < cell_count; ++i) {
                di[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)]
                    = derivative[static_cast<std::size_t>(i)];
            }
        }
        for (int i = 0; i < cell_count; ++i) {
            std::vector<Real> faces(static_cast<std::size_t>(cell_count + 1));
            for (int j = 0; j <= cell_count; ++j) {
                faces[static_cast<std::size_t>(j)] = j_component(i, j, 0);
            }
            const auto derivative = line.differentiate(faces);
            for (int j = 0; j < cell_count; ++j) {
                identity_max = std::max(
                    identity_max,
                    std::abs(di[static_cast<std::size_t>(j)]
                               [static_cast<std::size_t>(i)]
                        + derivative[static_cast<std::size_t>(j)]));
            }
        }
    }
    return {std::sqrt(squared_error / sample_count), identity_max};
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
    const auto extracted = extract_metric_field(
        result.metric, {1, 1, 1}, {2, 2, 2});
    WCNS_REQUIRE_NEAR(extracted.jacobian()(0, 0, 0), 15.0, 2.0e-13);
    WCNS_REQUIRE_NEAR(extracted.i_faces().x(2, 1, 1), 6.0, 2.0e-14);
    const auto payload = pack_metric_field(extracted);
    WCNS_REQUIRE(
        payload.size() == metric_field_payload_size({2, 2, 2}, 3));
    const auto unpacked = unpack_metric_field(
        AlgorithmProfileKind::PhengleiWcns, {2, 2, 2}, 3, payload);
    WCNS_REQUIRE_NEAR(unpacked.jacobian()(0, 0, 0), 15.0, 2.0e-13);
    WCNS_REQUIRE_NEAR(unpacked.k_faces().z(1, 1, 2), 7.0, 2.0e-14);
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
    const auto i_faces = result.metric.i_faces().x.interior_extent();
    for (int j = 0; j < i_faces.nj; ++j) {
        for (int i = 0; i < i_faces.ni; ++i) {
            WCNS_REQUIRE(result.metric.i_faces().z(i, j, 0) == 0.0);
        }
    }
    const auto j_faces = result.metric.j_faces().x.interior_extent();
    for (int j = 0; j < j_faces.nj; ++j) {
        for (int i = 0; i < j_faces.ni; ++i) {
            WCNS_REQUIRE(result.metric.j_faces().z(i, j, 0) == 0.0);
        }
    }

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

// 验收光滑扭曲网格的 profile 收敛趋势及离散度量恒等式自由流残差。
void test_warped_metric_convergence()
{
    using namespace wcns;
    const auto ph_coarse
        = warped_metric_error(AlgorithmProfileKind::PhengleiWcns, 16);
    const auto ph_fine
        = warped_metric_error(AlgorithmProfileKind::PhengleiWcns, 32);
    WCNS_REQUIRE(ph_coarse.jacobian_l2 / ph_fine.jacobian_l2 > 3.0);
    WCNS_REQUIRE(ph_coarse.identity_max < 2.0e-12);
    WCNS_REQUIRE(ph_fine.identity_max < 2.0e-12);

    const auto scmm_coarse
        = warped_metric_error(AlgorithmProfileKind::Scmm6Wcns, 16);
    const auto scmm_fine
        = warped_metric_error(AlgorithmProfileKind::Scmm6Wcns, 32);
    WCNS_REQUIRE(scmm_coarse.jacobian_l2 / scmm_fine.jacobian_l2 > 20.0);
    WCNS_REQUIRE(scmm_coarse.identity_max < 2.0e-11);
    WCNS_REQUIRE(scmm_fine.identity_max < 2.0e-11);
}
