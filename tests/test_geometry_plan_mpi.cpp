#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/geometry_halo.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/mpi_runtime.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <set>

// 验收各 MPI rank 从同一多块拓扑生成一致的分阶段几何消息和 donor rank。
int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: wcns_geometry_plan_tests <multi-2d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        using namespace wcns;
        MpiRuntime mpi(argc, argv);
        CgnsReader reader;
        auto mesh = reader.read_mesh(argv[1], 0, 3);
        std::vector<BlockLoad> loads;
        for (const auto& block : mesh.blocks()) {
            loads.push_back({block.id(), block.cell_extent().size()});
        }
        const auto distribution
            = BlockDistribution::balanced(std::move(loads), mpi.size());
        distribution.apply(mesh);
        const auto profile
            = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
        const auto plan = GeometryHaloPlan::build(mesh, profile);
        WCNS_REQUIRE(plan.exchanges().size() == 12);
        std::set<int> tags;
        bool valid = true;
        for (const auto& exchange : plan.exchanges()) {
            valid = valid && tags.insert(exchange.message_tag()).second;
            valid = valid
                && exchange.donor_rank == distribution.owner(exchange.donor_block);
            valid = valid && exchange.shared_face_owner == 0;
        }
        WCNS_REQUIRE(mpi.all_true(valid));
        std::cout << "geometry plan MPI tests passed on rank " << mpi.rank() << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
