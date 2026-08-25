#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/metrics.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/euler_solver.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

namespace {

bool occupied(const wcns::StructuredBlock& block, wcns::FaceLocation face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face == face) {
            return true;
        }
    }
    for (const auto& boundary : block.boundaries) {
        if (boundary.face == face) {
            return true;
        }
    }
    return false;
}

void add_farfield(wcns::StructuredBlock& block, wcns::FaceLocation face)
{
    using namespace wcns;
    if (occupied(block, face)) {
        return;
    }
    const auto vertices = block.vertex_extent();
    const auto cells = block.cell_extent();
    BoundaryPatch patch;
    patch.name = "test-farfield";
    patch.type = BoundaryType::Farfield;
    patch.face = face;
    if (face.axis == Axis::I) {
        const int vertex_i = face.side == Side::Lower ? 0 : vertices.ni - 1;
        const int cell_i = face.side == Side::Lower ? 0 : cells.ni;
        patch.vertex_range = {{vertex_i, 0, 0}, {vertex_i, vertices.nj - 1, 0}};
        patch.cell_face_range = {{cell_i, 0, 0}, {cell_i, cells.nj - 1, 0}};
    } else {
        const int vertex_j = face.side == Side::Lower ? 0 : vertices.nj - 1;
        const int cell_j = face.side == Side::Lower ? 0 : cells.nj;
        patch.vertex_range = {{0, vertex_j, 0}, {vertices.ni - 1, vertex_j, 0}};
        patch.cell_face_range = {{0, cell_j, 0}, {cells.ni - 1, cell_j, 0}};
    }
    block.boundaries.push_back(std::move(patch));
}

void prepare_block(wcns::StructuredBlock& block, const wcns::PrimitiveState& state)
{
    using namespace wcns;
    compute_metrics(block);
    add_farfield(block, {Axis::I, Side::Lower});
    add_farfield(block, {Axis::I, Side::Upper});
    add_farfield(block, {Axis::J, Side::Lower});
    add_farfield(block, {Axis::J, Side::Upper});
    const auto conservative = to_conservative(state);
    const auto extent = block.cell_extent();
    for (int j = 0; j < extent.nj; ++j) {
        for (int i = 0; i < extent.ni; ++i) {
            store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: wcns_euler_solver_tests <multi-2d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        using namespace wcns;
        MpiRuntime mpi(argc, argv);
        CgnsReader reader;
        auto global_mesh = reader.read_mesh(argv[1], 0, 3);
        std::vector<BlockLoad> loads;
        for (const auto& block : global_mesh.blocks()) {
            loads.push_back({block.id(), block.cell_extent().size()});
        }
        const auto distribution = BlockDistribution::balanced(std::move(loads), mpi.size());
        distribution.apply(global_mesh);
        const auto topology = DistributedTopology::build(global_mesh, distribution);

        const auto metadata = reader.read_metadata(argv[1]);
        std::vector<StructuredBlock> local_storage;
        for (const auto& zone : metadata.zones) {
            if (distribution.owner(zone.block_id) == mpi.rank()) {
                local_storage.push_back(reader.read_block(argv[1], zone, mpi.rank(), 3));
            }
        }
        LocalBlockSet local(mpi.rank(), std::move(local_storage), distribution);
        const PrimitiveState freestream {1.0, 0.35, -0.15, 0.0, 1.0};
        const auto expected = to_conservative(freestream);
        for (auto& block : local.blocks()) {
            prepare_block(block, freestream);
        }

        EulerSolver solver(
            mpi,
            local,
            topology,
            distribution.rank_count(),
            freestream,
            SpatialParameters {{}, {}, 0.25});
        solver.compute_residuals();
        WCNS_REQUIRE(solver.global_residual_l2() < 1.0e-11);
        const Real time_step = solver.advance_cfl();
        WCNS_REQUIRE(std::isfinite(time_step));
        WCNS_REQUIRE(time_step > 0.0);

        Real local_error = 0.0;
        for (const auto& block : local.blocks()) {
            const auto extent = block.cell_extent();
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const auto actual = load_conservative(
                        block.flow.conservative, {i, j, 0});
                    for (int component = 0; component < euler_components; ++component) {
                        local_error = std::max(
                            local_error,
                            std::abs(actual[static_cast<std::size_t>(component)]
                                - expected[static_cast<std::size_t>(component)]));
                    }
                }
            }
        }
        WCNS_REQUIRE(mpi.max(local_error) < 1.0e-11);
        if (mpi.rank() == 0) {
            std::cout << "Euler multiblock constant-state test passed with "
                      << mpi.size() << " ranks\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
