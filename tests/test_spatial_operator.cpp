#include "test_support.hpp"

#include <wcns/mesh/metrics.hpp>
#include <wcns/solver/boundary_conditions.hpp>
#include <wcns/solver/spatial_operator.hpp>
#include <wcns/solver/time_integrator.hpp>

#include <cmath>
#include <vector>

namespace {

wcns::StructuredBlock make_uniform_block()
{
    using namespace wcns;
    StructuredBlock block(0, "uniform", 0, 2, 2, {9, 7, 1}, 3);
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = 0.25 * static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = 0.2 * static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    compute_metrics(block);
    const auto cells = block.cell_extent();
    block.boundaries = {
        {"i-lower", BoundaryType::Farfield, {Axis::I, Side::Lower},
            {{0, 0, 0}, {0, vertices.nj - 1, 0}},
            {{0, 0, 0}, {0, cells.nj - 1, 0}}, {}},
        {"i-upper", BoundaryType::Farfield, {Axis::I, Side::Upper},
            {{vertices.ni - 1, 0, 0}, {vertices.ni - 1, vertices.nj - 1, 0}},
            {{cells.ni, 0, 0}, {cells.ni, cells.nj - 1, 0}}, {}},
        {"j-lower", BoundaryType::Farfield, {Axis::J, Side::Lower},
            {{0, 0, 0}, {vertices.ni - 1, 0, 0}},
            {{0, 0, 0}, {cells.ni - 1, 0, 0}}, {}},
        {"j-upper", BoundaryType::Farfield, {Axis::J, Side::Upper},
            {{0, vertices.nj - 1, 0}, {vertices.ni - 1, vertices.nj - 1, 0}},
            {{0, cells.nj, 0}, {cells.ni - 1, cells.nj, 0}}, {}},
    };
    return block;
}

void initialize(wcns::StructuredBlock& block, const wcns::PrimitiveState& state)
{
    const auto conservative = wcns::to_conservative(state);
    const auto extent = block.cell_extent();
    for (int j = 0; j < extent.nj; ++j) {
        for (int i = 0; i < extent.ni; ++i) {
            wcns::store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
}

} // namespace

void test_spatial_operator()
{
    using namespace wcns;
    auto block = make_uniform_block();
    const PrimitiveState freestream {1.1, 0.7, -0.2, 0.0, 1.0};
    const auto initial = to_conservative(freestream);
    initialize(block, freestream);

    const auto evaluate = [&] {
        update_primitive_interior(block);
        fill_physical_boundaries(block, freestream);
        compute_euler_residual(block);
    };
    evaluate();
    WCNS_REQUIRE(residual_l2(block) < 1.0e-12);
    const Real time_step = stable_time_step(block, 0.4);
    WCNS_REQUIRE(std::isfinite(time_step));
    WCNS_REQUIRE(time_step > 0.0);

    std::vector<StructuredBlock*> blocks {&block};
    advance_ssprk3(blocks, time_step, evaluate);
    const auto extent = block.cell_extent();
    for (int j = 0; j < extent.nj; ++j) {
        for (int i = 0; i < extent.ni; ++i) {
            const auto state = load_conservative(block.flow.conservative, {i, j, 0});
            for (int component = 0; component < euler_components; ++component) {
                WCNS_REQUIRE_NEAR(
                    state[static_cast<std::size_t>(component)],
                    initial[static_cast<std::size_t>(component)],
                    1.0e-13);
            }
        }
    }

    WCNS_REQUIRE_THROWS(std::invalid_argument, stable_time_step(block, 0.0));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        advance_ssprk3(blocks, -1.0, evaluate));
}
