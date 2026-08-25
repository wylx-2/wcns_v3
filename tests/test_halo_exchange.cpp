#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/parallel/mpi_runtime.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

namespace {

wcns::Real encoded_value(wcns::BlockId block, wcns::Index3 index, int component)
{
    return static_cast<wcns::Real>(
        100000 * block + 10000 * index.k + 1000 * index.j + 10 * index.i
        + component + 1);
}

void initialize_field(wcns::StructuredBlock& block)
{
    auto& field = block.flow.conservative;
    field.fill(-1.0);
    const auto extent = field.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < field.components(); ++component) {
                    field(i, j, k, component)
                        = encoded_value(block.id(), {i, j, k}, component);
                }
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: wcns_halo_exchange_tests <multi-2d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        wcns::MpiRuntime mpi(argc, argv);
        wcns::CgnsReader reader;
        auto global_mesh = reader.read_mesh(argv[1], 0, 3);

        std::vector<wcns::BlockLoad> loads;
        for (const auto& block : global_mesh.blocks()) {
            loads.push_back({block.id(), block.cell_extent().size()});
        }
        const auto distribution = wcns::BlockDistribution::balanced(
            std::move(loads), mpi.size());
        distribution.apply(global_mesh);
        const auto topology = wcns::DistributedTopology::build(
            global_mesh, distribution);

        const auto metadata = reader.read_metadata(argv[1]);
        std::vector<wcns::StructuredBlock> local_blocks;
        for (const auto& zone : metadata.zones) {
            if (distribution.owner(zone.block_id) == mpi.rank()) {
                local_blocks.push_back(
                    reader.read_block(argv[1], zone, mpi.rank(), 3));
            }
        }
        wcns::LocalBlockSet local(
            mpi.rank(), std::move(local_blocks), distribution);
        wcns::BlockFieldRegistry fields(wcns::euler_components);
        for (auto& block : local.blocks()) {
            initialize_field(block);
            fields.add(block.id(), block.flow.conservative);
        }

        wcns::HaloExchanger exchanger(mpi, topology, distribution.rank_count());
        exchanger.exchange(fields);

        bool correct = true;
        for (const auto& exchange : topology.exchanges()) {
            if (exchange.receiver_rank != mpi.rank()) {
                continue;
            }
            const auto& receiver = fields.field(exchange.halo.receiver_block);
            for (const auto& pair : exchange.halo.cell_pairs) {
                for (int component = 0; component < fields.components(); ++component) {
                    const auto actual = receiver(
                        pair.receiver_ghost.i,
                        pair.receiver_ghost.j,
                        pair.receiver_ghost.k,
                        component);
                    const auto expected = encoded_value(
                        exchange.halo.donor_block, pair.donor_interior, component);
                    correct = correct && actual == expected;
                }
            }
        }
        WCNS_REQUIRE(mpi.all_true(correct));
        if (mpi.rank() == 0) {
            std::cout << "halo exchange test passed with " << mpi.size() << " ranks\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

