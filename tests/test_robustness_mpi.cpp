#include "test_support.hpp"

#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/robustness.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace {

wcns::StructuredMesh make_two_rank_topology()
{
    using namespace wcns;
    StructuredBlock left(0, "left", 0, 2, 2, {7, 7, 1}, 3);
    StructuredBlock right(1, "right", 1, 2, 2, {7, 7, 1}, 3);
    left.connectivities.push_back({"left-right",
                                   0,
                                   1,
                                   1,
                                   {Axis::I, Side::Upper},
                                   {Axis::I, Side::Lower},
                                   {{6, 0, 0}, {6, 6, 0}},
                                   {{0, 0, 0}, {0, 6, 0}},
                                   {{5, 0, 0}, {5, 5, 0}},
                                   {{0, 0, 0}, {0, 5, 0}},
                                   {{6, 0, 0}, {6, 5, 0}},
                                   {{{1, 2, 3}}},
                                   3});
    right.connectivities.push_back({"right-left",
                                    1,
                                    0,
                                    0,
                                    {Axis::I, Side::Lower},
                                    {Axis::I, Side::Upper},
                                    {{0, 0, 0}, {0, 6, 0}},
                                    {{6, 0, 0}, {6, 6, 0}},
                                    {{0, 0, 0}, {0, 5, 0}},
                                    {{5, 0, 0}, {5, 5, 0}},
                                    {{0, 0, 0}, {0, 5, 0}},
                                    {{{1, 2, 3}}},
                                    3});
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

void test_reverse_owner_request_exchange(const wcns::MpiRuntime& mpi)
{
    using namespace wcns;
    if (mpi.size() < 2) return;
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto topology = make_two_rank_topology();
    const auto plan = FaceFluxHaloPlan::build(topology, profile, 77);
    FaceRobustnessRegistry registry;
    std::unique_ptr<FaceRobustnessField> levels;
    if (mpi.rank() == 0 || mpi.rank() == 1) {
        levels = std::make_unique<FaceRobustnessField>(
            topology.block(mpi.rank()).cell_extent(), 2, profile.kind());
        registry.add(mpi.rank(), *levels);
    }
    if (mpi.rank() == 1) {
        WCNS_REQUIRE(levels->request_next(Axis::I, {0, 2, 0}, 3));
        WCNS_REQUIRE(levels->request_next(Axis::I, {0, 2, 0}, 3));
    }
    FaceRobustnessExchanger exchanger(mpi, plan);
    WCNS_REQUIRE(exchanger.exchange_requests(registry, 1, 1));
    if (mpi.rank() == 0) {
        WCNS_REQUIRE(levels->level(Axis::I, {6, 2, 0}) == 2);
    }
    WCNS_REQUIRE(!exchanger.exchange_requests(registry, 1, 2));
}

} // namespace

int main(int argc, char** argv)
{
    try {
        using namespace wcns;
        MpiRuntime mpi(argc, argv);
        test_reverse_owner_request_exchange(mpi);
        StructuredBlock block(mpi.rank(), "robust-retry", mpi.rank(), 2, 2, {5, 5, 1}, 3);
        const ConservativeState state {1.0, 0.0, 0.0, 0.0, 2.5};
        const auto cells = block.cell_extent();
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                store_state(block.flow.conservative, {i, j, 0}, state);
            }
        }
        std::vector<StructuredBlock*> blocks {&block};

        GasModelInput gas_input;
        gas_input.specific_gas_constant = 287.0;
        const auto gas = GasModel::from_input(gas_input);
        const auto reference
            = ReferenceScales::derive({340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
        const NumericalFloors floors;
        const auto profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
        ReconstructionConfig reconstruction;
        reconstruction.scheme = "weno_z";
        reconstruction.variables = ReconstructionVariables::Characteristic;
        RiemannConfig riemann;
        riemann.scheme = "hllc";
        const auto ladder = RobustnessLadder::build(reconstruction, riemann);
        RobustnessConfig config;
        config.enabled = true;
        RobustnessDiagnostics diagnostics;
        std::size_t evaluations = 0;
        const auto evaluator = [&](Real, int, const BlockFaceRobustnessMap&) {
            ++evaluations;
            block.flow.residual.fill(0.0);
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    block.flow.residual(i, j, 0, density) = -2.0;
                }
            }
        };
        const StructuredMesh topology;
        const Real accepted = advance_ssprk3_with_robustness(mpi,
                                                             blocks,
                                                             topology,
                                                             profile,
                                                             FluxDifferenceMode::Profile,
                                                             gas,
                                                             reference,
                                                             floors,
                                                             config,
                                                             ladder,
                                                             1.0,
                                                             0.0,
                                                             evaluator,
                                                             diagnostics);
        WCNS_REQUIRE_NEAR(accepted, 0.25, 0.0);
        WCNS_REQUIRE(diagnostics.step_retries == 2);
        WCNS_REQUIRE(diagnostics.local_recomputations == 6);
        WCNS_REQUIRE(diagnostics.troubled_cells > 0);
        WCNS_REQUIRE(evaluations == 11);
        WCNS_REQUIRE(diagnostics.face_levels[0] == 120);
        WCNS_REQUIRE(diagnostics.face_levels[1] == 0);
        WCNS_REQUIRE_NEAR(block.flow.conservative(0, 0, 0, density), 0.5, 1.0e-14);

        if (mpi.rank() == 0) {
            std::cout << "WCNS robustness retry test passed\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
