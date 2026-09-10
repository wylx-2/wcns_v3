#include "test_support.hpp"

#include <wcns/solver/robustness.hpp>

#include <limits>
#include <stdexcept>
#include <vector>

namespace {

wcns::StructuredBlock make_robustness_block()
{
    return wcns::StructuredBlock(0, "robustness", 0, 2, 2, {9, 7, 1}, 3);
}

wcns::ConnectivityPatch lower_i_connection()
{
    using namespace wcns;
    return {
        "lower-i", 0, 1, 0,
        {Axis::I, Side::Lower}, {Axis::I, Side::Upper},
        {{0, 0, 0}, {0, 6, 0}}, {{8, 0, 0}, {8, 6, 0}},
        {{0, 0, 0}, {0, 5, 0}}, {{7, 0, 0}, {7, 5, 0}},
        {{0, 0, 0}, {0, 5, 0}}, {{{1, 2, 3}}}, 3};
}

void initialize_positive_state(wcns::StructuredBlock& block)
{
    const wcns::ConservativeState state {1.0, 0.0, 0.0, 0.0, 2.5};
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            wcns::store_state(block.flow.conservative, {i, j, 0}, state);
        }
    }
}

} // namespace

// 验收 Q 阶段降阶梯、候选态事务语义，以及 troubled-cell 到真实离散支持域的传播。
void test_robustness()
{
    using namespace wcns;

    RobustnessConfig config;
    config.validate();
    config.time_step_reduction = 1.0;
    WCNS_REQUIRE_THROWS(std::invalid_argument, config.validate());

    ReconstructionConfig reconstruction;
    reconstruction.scheme = "weno_z";
    reconstruction.variables = ReconstructionVariables::Characteristic;
    RiemannConfig riemann;
    riemann.scheme = "hllc";
    const auto ladder = RobustnessLadder::build(reconstruction, riemann);
    WCNS_REQUIRE(ladder.size() == 4);
    WCNS_REQUIRE(ladder.strategy(0).reconstruction.scheme == "weno_z");
    WCNS_REQUIRE(
        ladder.strategy(1).reconstruction.variables
        == ReconstructionVariables::Primitive);
    WCNS_REQUIRE(ladder.strategy(2).reconstruction.scheme == "linear5");
    WCNS_REQUIRE(ladder.strategy(3).reconstruction.scheme == "zero_order");
    WCNS_REQUIRE(ladder.strategy(3).force_rusanov);

    ReconstructionConfig already_linear;
    already_linear.scheme = "linear5";
    already_linear.variables = ReconstructionVariables::Primitive;
    RiemannConfig rusanov;
    const auto compact = RobustnessLadder::build(already_linear, rusanov);
    WCNS_REQUIRE(compact.size() == 2);
    WCNS_REQUIRE(!compact.strategy(1).force_rusanov);

    auto block = make_robustness_block();
    initialize_positive_state(block);
    std::vector<StructuredBlock*> blocks {&block};
    const auto initial = capture_conservative_state(blocks);
    block.flow.residual.fill(0.0);
    block.flow.residual(2, 3, 0, density) = -2.0;
    const auto invalid_candidate
        = form_ssprk_candidate(blocks, initial, 1.0, 0.0, 1.0);
    GasModelInput gas_input;
    gas_input.specific_gas_constant = 287.0;
    const auto gas = GasModel::from_input(gas_input);
    const auto reference = ReferenceScales::derive(
        {340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
    const NumericalFloors floors;
    const auto validation = validate_candidate_state(
        invalid_candidate, blocks, gas, reference, floors, 1, 0.0);
    WCNS_REQUIRE(!validation.valid());
    WCNS_REQUIRE(validation.troubled_cells.size() == 1);
    WCNS_REQUIRE(validation.troubled_cells.front().cell == (Index3 {2, 3, 0}));
    WCNS_REQUIRE(validation.troubled_cells.front().reason == "density_floor");
    WCNS_REQUIRE_NEAR(block.flow.conservative(2, 3, 0, density), 1.0, 0.0);

    auto corrupted = initial;
    corrupted.front().values.front() = std::numeric_limits<Real>::quiet_NaN();
    const auto non_finite = validate_candidate_state(
        corrupted, blocks, gas, reference, floors, 2, 0.1);
    WCNS_REQUIRE(!non_finite.valid());
    WCNS_REQUIRE(
        non_finite.troubled_cells.front().reason
        == "non_finite_conservative");

    block.connectivities.push_back(lower_i_connection());
    TroubledCell troubled;
    troubled.block = block.id();
    troubled.rank = block.owner_rank();
    troubled.cell = {0, 2, 0};
    troubled.rk_stage = 1;
    troubled.stage_time = 0.0;
    const auto ph = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    FaceRobustnessField ph_levels(
        block.cell_extent(), block.cell_dimension(), ph.kind());
    WCNS_REQUIRE(request_troubled_cell_support(
        block, ph, FluxDifferenceMode::Profile, {troubled}, ph_levels, 3));
    for (const int face : {-1, 0, 1, 2}) {
        WCNS_REQUIRE(ph_levels.level(Axis::I, {face, 2, 0}) == 1);
    }

    const auto scmm = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    FaceRobustnessField scmm_levels(
        block.cell_extent(), block.cell_dimension(), scmm.kind());
    WCNS_REQUIRE(request_troubled_cell_support(
        block, scmm, FluxDifferenceMode::Profile, {troubled}, scmm_levels, 3));
    for (const int face : {-2, -1, 0, 1, 2, 3}) {
        WCNS_REQUIRE(scmm_levels.level(Axis::I, {face, 2, 0}) == 1);
    }

    FaceRobustnessField two_point_levels(
        block.cell_extent(), block.cell_dimension(), ph.kind());
    TroubledCell interior = troubled;
    interior.cell = {3, 2, 0};
    WCNS_REQUIRE(request_troubled_cell_support(
        block, ph, FluxDifferenceMode::ConservativeTwoPoint,
        {interior}, two_point_levels, 3));
    WCNS_REQUIRE(two_point_levels.level(Axis::I, {3, 2, 0}) == 1);
    WCNS_REQUIRE(two_point_levels.level(Axis::I, {4, 2, 0}) == 1);
    WCNS_REQUIRE(two_point_levels.level(Axis::I, {2, 2, 0}) == 1);
    WCNS_REQUIRE(two_point_levels.level(Axis::I, {5, 2, 0}) == 1);
    WCNS_REQUIRE(two_point_levels.level(Axis::I, {1, 2, 0}) == 0);
}
