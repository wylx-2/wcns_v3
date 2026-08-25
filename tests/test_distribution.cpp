#include "test_support.hpp"

#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

namespace {

wcns::StructuredMesh make_two_block_mesh()
{
    using namespace wcns;
    StructuredBlock left(0, "left", 0, 2, 2, {5, 4, 1}, 3);
    StructuredBlock right(1, "right", 0, 2, 2, {4, 5, 1}, 3);
    left.connectivities.push_back({
        "left-to-right",
        0,
        1,
        0,
        {Axis::I, Side::Upper},
        {Axis::J, Side::Lower},
        {{4, 0, 0}, {4, 3, 0}},
        {{3, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{{2, -1, 3}}},
        3,
    });
    right.connectivities.push_back({
        "right-to-left",
        1,
        0,
        0,
        {Axis::J, Side::Lower},
        {Axis::I, Side::Upper},
        {{3, 0, 0}, {0, 0, 0}},
        {{4, 0, 0}, {4, 3, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{{-2, 1, 3}}},
        3,
    });
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

} // namespace

void test_distribution()
{
    using namespace wcns;

    const auto distribution = BlockDistribution::balanced(
        {{0, 100}, {1, 60}, {2, 40}, {3, 20}}, 2);
    WCNS_REQUIRE(distribution.owner(0) == 0);
    WCNS_REQUIRE(distribution.owner(1) == 1);
    WCNS_REQUIRE(distribution.owner(2) == 1);
    WCNS_REQUIRE(distribution.owner(3) == 0);
    WCNS_REQUIRE(distribution.rank_loads() == (std::vector<std::size_t> {120, 100}));
    WCNS_REQUIRE(distribution.local_blocks(0) == (std::vector<BlockId> {0, 3}));
    WCNS_REQUIRE_THROWS(std::out_of_range, distribution.owner(9));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, BlockDistribution::balanced({{0, 1}}, 0));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, BlockDistribution::balanced({{0, 1}, {0, 2}}, 2));

    auto mesh = make_two_block_mesh();
    const auto split = BlockDistribution::balanced({{0, 12}, {1, 12}}, 2);
    split.apply(mesh);
    WCNS_REQUIRE(mesh.block(0).owner_rank() == 0);
    WCNS_REQUIRE(mesh.block(1).owner_rank() == 1);
    WCNS_REQUIRE(mesh.block(0).connectivities.front().donor_rank == 1);

    const auto topology = DistributedTopology::build(mesh, split);
    WCNS_REQUIRE(topology.exchanges().size() == 2);
    WCNS_REQUIRE(topology.exchanges()[0].connection == 0);
    WCNS_REQUIRE(topology.exchanges()[1].connection == 0);
    WCNS_REQUIRE(topology.exchanges()[0].message_tag() == 1024);
    WCNS_REQUIRE(topology.exchanges()[1].message_tag() == 1025);
    WCNS_REQUIRE(topology.receives(0).size() == 1);
    WCNS_REQUIRE(topology.sends(0).size() == 1);
    WCNS_REQUIRE(topology.local_copies(0).empty());

    auto local_mesh = make_two_block_mesh();
    const auto one_rank = BlockDistribution::balanced({{0, 12}, {1, 12}}, 1);
    one_rank.apply(local_mesh);
    const auto local_topology = DistributedTopology::build(local_mesh, one_rank);
    WCNS_REQUIRE(local_topology.local_copies(0).size() == 2);

    std::vector<StructuredBlock> local_blocks;
    local_blocks.emplace_back(0, "local", 0, 2, 2, Extent3 {3, 3, 1}, 3);
    LocalBlockSet local_set(0, std::move(local_blocks), split);
    WCNS_REQUIRE(local_set.contains(0));
    WCNS_REQUIRE(local_set.block(0).owner_rank() == 0);
    WCNS_REQUIRE_THROWS(std::out_of_range, local_set.block(1));
}

