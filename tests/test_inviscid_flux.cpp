#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/metrics.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/spatial_operator.hpp>

#include <cmath>
#include <utility>
#include <vector>

namespace {

wcns::GasModel flux_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales flux_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive(
        {340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
}

wcns::StructuredBlock make_flux_block(bool wall)
{
    using namespace wcns;
    StructuredBlock block(0, "flux", 0, 2, 2, {10, 9, 1}, 3);
    const auto vertices = block.vertex_extent();
    const auto cells = block.cell_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    compute_metrics(block);
    block.boundaries = {
        {"i-lower", wall ? BoundaryType::SlipWall : BoundaryType::Farfield,
            {Axis::I, Side::Lower}, {{0, 0, 0}, {0, vertices.nj - 1, 0}},
            {{0, 0, 0}, {0, cells.nj - 1, 0}},
            {{0, 0, 0}, {0, cells.nj - 1, 0}}, {}},
        {"i-upper", BoundaryType::Farfield, {Axis::I, Side::Upper},
            {{vertices.ni - 1, 0, 0}, {vertices.ni - 1, vertices.nj - 1, 0}},
            {{cells.ni - 1, 0, 0}, {cells.ni - 1, cells.nj - 1, 0}},
            {{cells.ni, 0, 0}, {cells.ni, cells.nj - 1, 0}}, {}},
        {"j-lower", BoundaryType::Farfield, {Axis::J, Side::Lower},
            {{0, 0, 0}, {vertices.ni - 1, 0, 0}},
            {{0, 0, 0}, {cells.ni - 1, 0, 0}},
            {{0, 0, 0}, {cells.ni - 1, 0, 0}}, {}},
        {"j-upper", BoundaryType::Farfield, {Axis::J, Side::Upper},
            {{0, vertices.nj - 1, 0}, {vertices.ni - 1, vertices.nj - 1, 0}},
            {{0, cells.nj - 1, 0}, {cells.ni - 1, cells.nj - 1, 0}},
            {{0, cells.nj, 0}, {cells.ni - 1, cells.nj, 0}}, {}},
    };
    return block;
}

wcns::BoundaryDataMap flux_boundaries(
    const wcns::StructuredBlock& block,
    const wcns::TemperaturePrimitiveState& target)
{
    wcns::BoundaryDataMap result;
    for (const auto& patch : block.boundaries) {
        wcns::BoundaryData data;
        if (patch.type == wcns::BoundaryType::Farfield
            || patch.type == wcns::BoundaryType::Inflow) {
            data.target_state = target;
        }
        result.emplace(patch.name, data);
    }
    return result;
}

void initialize_flux_state(
    wcns::StructuredBlock& block,
    const wcns::TemperaturePrimitiveState& state,
    const wcns::GasModel& gas,
    const wcns::ReferenceScales& reference,
    const wcns::NumericalFloors& floors)
{
    const auto conservative = wcns::thermodynamic_conservative(
        state, gas, reference, floors, block.cell_dimension());
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            wcns::store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
    wcns::update_temperature_primitive_interior(
        block, gas, reference, floors);
}

wcns::StructuredMesh make_flux_mesh()
{
    using namespace wcns;
    StructuredBlock left(0, "left", 0, 2, 2, {7, 7, 1}, 3);
    StructuredBlock right(1, "right", 0, 2, 2, {7, 7, 1}, 3);
    left.connectivities.push_back({
        "left-right", 0, 1, 0,
        {Axis::I, Side::Upper}, {Axis::I, Side::Lower},
        {{6, 0, 0}, {6, 6, 0}}, {{0, 0, 0}, {0, 6, 0}},
        {{5, 0, 0}, {5, 5, 0}}, {{0, 0, 0}, {0, 5, 0}},
        {{6, 0, 0}, {6, 5, 0}}, {{{1, 2, 3}}}, 3});
    right.connectivities.push_back({
        "right-left", 1, 0, 0,
        {Axis::I, Side::Lower}, {Axis::I, Side::Upper},
        {{0, 0, 0}, {0, 6, 0}}, {{6, 0, 0}, {6, 6, 0}},
        {{0, 0, 0}, {0, 5, 0}}, {{5, 0, 0}, {5, 5, 0}},
        {{0, 0, 0}, {0, 5, 0}}, {{{1, 2, 3}}}, 3});
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

} // namespace

// 验收两套 profile 的真实面 Rusanov 通量和高阶散度在笛卡尔自由流上保持常量。
void test_wcns_inviscid_freestream()
{
    using namespace wcns;
    const auto gas = flux_gas();
    const auto reference = flux_reference(gas);
    const NumericalFloors floors;
    const TemperaturePrimitiveState state {1.1, 0.7, -0.2, 0.0, 1.0};
    for (const auto kind : {
             AlgorithmProfileKind::PhengleiWcns,
             AlgorithmProfileKind::Scmm6Wcns}) {
        auto block = make_flux_block(false);
        initialize_flux_state(block, state, gas, reference, floors);
        const auto data = flux_boundaries(block, state);
        const auto ghost_result = PhysicalGhostStateOperator::fill(
            block, data, gas, reference, floors, 1);
        WCNS_REQUIRE(ghost_result.version == 1);
        const auto profile = ProfileFactory::create(kind);
        const auto metric = initialize_metric_field(block, profile).metric;
        ReconstructionConfig reconstruction;
        reconstruction.kind = ReconstructionKind::Linear5;
        ReconstructionDiagnostics diagnostics;
        const RiemannSolver riemann;
        const auto flux = compute_inviscid_face_fluxes(
            block, metric, profile, reconstruction, riemann, gas, reference,
            floors, data, {}, 1, diagnostics);
        compute_wcns_inviscid_residual(block, metric, flux, profile);
        WCNS_REQUIRE(residual_l2(block) < 2.0e-11);
        WCNS_REQUIRE(diagnostics.linear_faces
            == static_cast<std::size_t>(
                (block.cell_extent().ni + 1) * block.cell_extent().nj
                + block.cell_extent().ni * (block.cell_extent().nj + 1)));
    }
}

// 验收强固壁面状态经 Rusanov 后具有严格零质量和零切向动量通量。
void test_wcns_strong_wall_flux()
{
    using namespace wcns;
    const auto gas = flux_gas();
    const auto reference = flux_reference(gas);
    const NumericalFloors floors;
    const TemperaturePrimitiveState state {1.0, 0.6, 0.4, 0.0, 1.0};
    auto block = make_flux_block(true);
    initialize_flux_state(block, state, gas, reference, floors);
    const auto data = flux_boundaries(block, state);
    const auto ghost_result = PhysicalGhostStateOperator::fill(
        block, data, gas, reference, floors, 2);
    WCNS_REQUIRE(ghost_result.version == 2);
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto metric = initialize_metric_field(block, profile).metric;
    ReconstructionConfig reconstruction;
    reconstruction.kind = ReconstructionKind::Linear5;
    ReconstructionDiagnostics diagnostics;
    const RiemannSolver riemann;
    const auto flux = compute_inviscid_face_fluxes(
        block, metric, profile, reconstruction, riemann, gas, reference,
        floors, data, {}, 2, diagnostics);
    for (int j = 0; j < block.cell_extent().nj; ++j) {
        WCNS_REQUIRE_NEAR(flux.field(Axis::I)(0, j, 0, density), 0.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(flux.field(Axis::I)(0, j, 0, momentum_y), 0.0, 1.0e-14);
    }
}

// 验收 PH 层 0--1、SCMM6 层 0--2 的所有者方向、索引对和消息标签计划。
void test_face_flux_halo_plan()
{
    using namespace wcns;
    const auto mesh = make_flux_mesh();
    const auto ph = FaceFluxHaloPlan::build(
        mesh, ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns), 3);
    WCNS_REQUIRE(ph.exchanges().size() == 2);
    WCNS_REQUIRE(ph.exchanges()[0].pairs.size() == 6);
    WCNS_REQUIRE(ph.exchanges()[1].pairs.size() == 12);
    WCNS_REQUIRE(ph.exchanges()[0].shared_face_owner == 0);
    WCNS_REQUIRE(ph.exchanges()[0].message_tag()
        != ph.exchanges()[1].message_tag());

    const auto scmm = FaceFluxHaloPlan::build(
        mesh, ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns), 4);
    WCNS_REQUIRE(scmm.exchanges()[0].pairs.size() == 12);
    WCNS_REQUIRE(scmm.exchanges()[1].pairs.size() == 18);
    WCNS_REQUIRE(scmm.exchanges()[0].pairs.front().layer == 1);
    WCNS_REQUIRE(scmm.exchanges()[1].pairs.front().layer == 0);
    WCNS_REQUIRE_THROWS(
        TopologyError,
        FaceFluxHaloPlan::build(
            mesh, ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns), 0));
}
