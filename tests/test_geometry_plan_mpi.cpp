#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/geometry_halo.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/inviscid_flux.hpp>

#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <set>
#include <unordered_map>
#include <tuple>
#include <utility>

namespace {

wcns::Real synthetic_flux(
    wcns::BlockId block, wcns::Axis axis, wcns::Index3 index, int component)
{
    return 10000.0 * static_cast<wcns::Real>(block)
        + 1000.0 * static_cast<wcns::Real>(static_cast<int>(axis))
        + 100.0 * index.i + 10.0 * index.j + index.k + 0.01 * component;
}

wcns::ConservativeState expected_received(
    const wcns::FaceFluxExchangeDescriptor& descriptor,
    const wcns::FaceFluxPair& pair)
{
    wcns::ConservativeState donor {};
    for (int component = 0; component < wcns::euler_components; ++component) {
        donor[static_cast<std::size_t>(component)] = synthetic_flux(
            descriptor.donor_block, descriptor.donor_axis, pair.donor, component);
    }
    const auto momentum = descriptor.periodic.inverse().apply_vector(
        {{donor[1], donor[2], donor[3]}});
    return {
        descriptor.orientation * donor[0],
        descriptor.orientation * momentum[0],
        descriptor.orientation * momentum[1],
        descriptor.orientation * momentum[2],
        descriptor.orientation * donor[4],
    };
}

} // namespace

// 验收各 MPI rank 的几何计划及 profile 专属真实面通量层交换与变换。
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

        const std::uint64_t flux_version = 11;
        const auto flux_plan = FaceFluxHaloPlan::build(
            mesh, profile, flux_version);
        std::unordered_map<BlockId, InviscidFaceFluxField> flux_fields;
        FaceFluxFieldRegistry registry;
        for (const auto& block : mesh.blocks()) {
            if (block.owner_rank() != mpi.rank()) {
                continue;
            }
            auto [iterator, inserted] = flux_fields.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(
                    block.cell_extent(), block.cell_dimension(),
                    profile.kind(), flux_version));
            WCNS_REQUIRE(inserted);
            for (const Axis axis : {Axis::I, Axis::J}) {
                auto& field = iterator->second.field(axis);
                const auto extent = field.interior_extent();
                for (int k = 0; k < extent.nk; ++k) {
                    for (int j = 0; j < extent.nj; ++j) {
                        for (int i = 0; i < extent.ni; ++i) {
                            for (int component = 0;
                                 component < euler_components; ++component) {
                                field(i, j, k, component) = synthetic_flux(
                                    block.id(), axis, {i, j, k}, component);
                            }
                        }
                    }
                }
            }
            registry.add(block.id(), iterator->second);
        }
        FaceFluxHaloExchanger(mpi, flux_plan).exchange(registry);
        bool flux_valid = true;
        for (const auto& descriptor : flux_plan.exchanges()) {
            if (descriptor.receiver_rank != mpi.rank()) {
                continue;
            }
            const auto& receiver = registry.field(descriptor.receiver_block)
                                       .field(descriptor.receiver_axis);
            for (const auto& pair : descriptor.pairs) {
                const auto expected = expected_received(descriptor, pair);
                for (int component = 0; component < euler_components; ++component) {
                    flux_valid = flux_valid
                        && std::abs(receiver(
                            pair.receiver.i, pair.receiver.j, pair.receiver.k,
                            component)
                            - expected[static_cast<std::size_t>(component)])
                            < 1.0e-12;
                }
            }
        }
        WCNS_REQUIRE(mpi.all_true(flux_valid));
        std::cout << "geometry plan MPI tests passed on rank " << mpi.rank() << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
