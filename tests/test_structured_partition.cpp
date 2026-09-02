#include "test_support.hpp"

#include <wcns/runtime/structured_partition.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

void test_structured_partition()
{
    using namespace wcns;
    PartitionConfig config;
    config.mode = PartitionMode::AutoSplit;
    config.allow_idle_ranks = false;
    config.max_load_ratio = 1.2;
    config.min_cells_per_active_direction = 8;

    const auto plan = StructuredPartitionPlan::build(
        {{7, "single", 2, {64, 32, 1}}},
        8,
        config);
    WCNS_REQUIRE(plan.leaves().size() >= 8);
    WCNS_REQUIRE(plan.maximum_feasible_leaf_count() == 32);
    WCNS_REQUIRE(plan.digest() != 0);
    std::set<RankId> owners;
    std::size_t covered = 0;
    for (const auto& leaf : plan.leaves()) {
        WCNS_REQUIRE(leaf.cell_extent().ni >= 8);
        WCNS_REQUIRE(leaf.cell_extent().nj >= 8);
        WCNS_REQUIRE(leaf.cell_extent().nk == 1);
        WCNS_REQUIRE(leaf.vertex_extent().ni == leaf.cell_extent().ni + 1);
        WCNS_REQUIRE(leaf.vertex_extent().nj == leaf.cell_extent().nj + 1);
        owners.insert(leaf.owner);
        covered += leaf.cell_count();
    }
    WCNS_REQUIRE(covered == static_cast<std::size_t>(64 * 32));
    WCNS_REQUIRE(owners.size() == 8);

    const auto repeated = StructuredPartitionPlan::build(
        {{7, "single", 2, {64, 32, 1}}},
        8,
        config);
    WCNS_REQUIRE(repeated.summary() == plan.summary());
    WCNS_REQUIRE(repeated.digest() == plan.digest());

    WCNS_REQUIRE_THROWS(
        CaseConfigurationError,
        StructuredPartitionPlan::build(
            {{0, "too-small", 2, {12, 8, 1}}},
            2,
            config));

    config.allow_idle_ranks = true;
    const auto idle = StructuredPartitionPlan::build(
        {{0, "too-small", 2, {12, 8, 1}}},
        2,
        config);
    WCNS_REQUIRE(idle.leaves().size() == 1);
    WCNS_REQUIRE(idle.distribution().rank_loads().at(1) == 0);

    config.mode = PartitionMode::ZonesOnly;
    config.allow_idle_ranks = false;
    WCNS_REQUIRE_THROWS(
        CaseConfigurationError,
        StructuredPartitionPlan::build(
            {{0, "left", 2, {16, 16, 1}}},
            2,
            config));

    const auto native_zones = StructuredPartitionPlan::build(
        {
            {0, "left", 2, {16, 16, 1}},
            {1, "right", 2, {16, 16, 1}},
        },
        2,
        config);
    WCNS_REQUIRE(native_zones.leaves().size() == 2);
    WCNS_REQUIRE(
        *std::max_element(
            native_zones.distribution().rank_loads().begin(),
            native_zones.distribution().rank_loads().end())
        == 256);
}
