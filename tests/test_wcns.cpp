#include "test_support.hpp"

#include <wcns/mesh/metrics.hpp>
#include <wcns/solver/boundary_conditions.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <array>
#include <cmath>

namespace {

void set_cartesian_coordinates(wcns::StructuredBlock& block)
{
    const auto extent = block.vertex_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                block.coordinates.x(i, j, k) = static_cast<wcns::Real>(i);
                block.coordinates.y(i, j, k) = static_cast<wcns::Real>(j);
                block.coordinates.z(i, j, k) = static_cast<wcns::Real>(k);
            }
        }
    }
}

} // namespace

// 验收现有 WCNS 重构、物理边界 ghost 和守恒状态同步。
void test_wcns()
{
    using namespace wcns;

    const std::array<Real, 6> constant {{3.5, 3.5, 3.5, 3.5, 3.5, 3.5}};
    const auto constant_face = wcns5_reconstruct(constant);
    WCNS_REQUIRE_NEAR(constant_face.left, 3.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(constant_face.right, 3.5, 1.0e-14);

    const std::array<Real, 6> linear {{-3.0, -1.0, 1.0, 3.0, 5.0, 7.0}};
    const auto linear_face = wcns5_reconstruct(linear);
    WCNS_REQUIRE_NEAR(linear_face.left, 2.0, 1.0e-13);
    WCNS_REQUIRE_NEAR(linear_face.right, 2.0, 1.0e-13);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        wcns5_reconstruct(linear, WcnsParameters {0.0, 2}));

    const auto smooth_error = [](Real spacing) {
        std::array<Real, 6> values {};
        for (int point = 0; point < 6; ++point) {
            values[static_cast<std::size_t>(point)]
                = std::exp(static_cast<Real>(point - 2) * spacing);
        }
        return std::abs(
            wcns5_reconstruct(values).left - std::exp(0.5 * spacing));
    };
    WCNS_REQUIRE(smooth_error(0.1) < smooth_error(0.2));

    Field<Real> primitive({4, 1, 1}, euler_components, 3);
    for (int i = -3; i < 7; ++i) {
        for (int component = 0; component < euler_components; ++component) {
            primitive(i, 0, 0, component)
                = static_cast<Real>(component + 1) * static_cast<Real>(i);
        }
    }
    const auto reconstructed = reconstruct_euler_face(primitive, Axis::I, {2, 0, 0});
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            reconstructed.left[static_cast<std::size_t>(component)],
            1.5 * static_cast<Real>(component + 1),
            1.0e-12);
        WCNS_REQUIRE_NEAR(
            reconstructed.right[static_cast<std::size_t>(component)],
            1.5 * static_cast<Real>(component + 1),
            1.0e-12);
    }

    StructuredBlock block(0, "boundary", 0, 2, 2, {5, 4, 1}, 3);
    set_cartesian_coordinates(block);
    compute_metrics(block);
    const PrimitiveState state {1.0, 2.0, 1.0, 0.0, 1.0};
    const auto conservative = to_conservative(state);
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
    block.boundaries.push_back({
        "i-lower-wall",
        BoundaryType::SlipWall,
        {Axis::I, Side::Lower},
        {{0, 0, 0}, {0, 3, 0}},
        {{0, 0, 0}, {0, 2, 0}},
        {{0, 0, 0}, {0, 2, 0}},
        {},
    });
    block.boundaries.push_back({
        "i-upper-outflow",
        BoundaryType::Outflow,
        {Axis::I, Side::Upper},
        {{4, 0, 0}, {4, 3, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{4, 0, 0}, {4, 2, 0}},
        {},
    });
    update_primitive_interior(block);
    fill_physical_boundaries(block, state);

    for (int layer = 1; layer <= 3; ++layer) {
        const auto wall = load_primitive(block.flow.primitive, {-layer, 1, 0});
        WCNS_REQUIRE_NEAR(wall[primitive_density], 1.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(wall[velocity_x], -2.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(wall[velocity_y], 1.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(wall[pressure], 1.0, 1.0e-14);
        const auto wall_conservative
            = load_conservative(block.flow.conservative, {-layer, 1, 0});
        const auto expected_wall = to_conservative(wall);
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(
                wall_conservative[static_cast<std::size_t>(component)],
                expected_wall[static_cast<std::size_t>(component)],
                1.0e-14);
        }

        const auto outflow
            = load_primitive(block.flow.primitive, {cells.ni - 1 + layer, 1, 0});
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(
                outflow[static_cast<std::size_t>(component)],
                state[static_cast<std::size_t>(component)],
                1.0e-14);
        }
    }
}
