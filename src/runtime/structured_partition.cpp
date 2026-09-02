#include <wcns/runtime/structured_partition.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace wcns {
namespace {

struct SplitCandidate {
    std::size_t leaf_index = 0;
    int axis = 0;
    int cut = 0;
    std::size_t largest_child = 0;
    std::size_t interface_cells = 0;
    int axis_length = 0;
    bool valid = false;
};

bool candidate_less(const SplitCandidate& lhs, const SplitCandidate& rhs)
{
    return std::make_tuple(
        lhs.largest_child,
        lhs.interface_cells,
        -lhs.axis_length,
        lhs.axis,
        lhs.cut)
        < std::make_tuple(
            rhs.largest_child,
            rhs.interface_cells,
            -rhs.axis_length,
            rhs.axis,
            rhs.cut);
}

SplitCandidate best_split(
    const PartitionLeaf& leaf,
    std::size_t leaf_index,
    int minimum_cells)
{
    SplitCandidate best;
    const auto extent = leaf.cell_extent();
    for (int axis = 0; axis < leaf.cell_dimension; ++axis) {
        const auto axis_index = static_cast<std::size_t>(axis);
        const int length = extent[axis_index];
        if (length < 2 * minimum_cells) continue;
        const int local_cut = length / 2;
        const int cut = leaf.cells.begin[axis_index] + local_cut;
        auto left_extent = extent;
        auto right_extent = extent;
        left_extent[axis_index] = local_cut;
        right_extent[axis_index] = length - local_cut;
        auto interface_extent = extent;
        interface_extent[axis_index] = 1;
        SplitCandidate candidate {
            leaf_index,
            axis,
            cut,
            std::max(left_extent.size(), right_extent.size()),
            interface_extent.size(),
            length,
            true,
        };
        if (!best.valid || candidate_less(candidate, best)) {
            best = candidate;
        }
    }
    return best;
}

bool boxes_overlap(const CellBox& lhs, const CellBox& rhs, int dimension)
{
    for (int axis = 0; axis < dimension; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        if (lhs.end[a] <= rhs.begin[a] || rhs.end[a] <= lhs.begin[a]) {
            return false;
        }
    }
    return true;
}

std::size_t maximum_leaves_for_zone(
    const PartitionZone& zone,
    int minimum_cells)
{
    std::size_t result = 1;
    for (int axis = 0; axis < zone.cell_dimension; ++axis) {
        const int chunks = zone.cell_extent[static_cast<std::size_t>(axis)]
            / minimum_cells;
        if (chunks <= 0) return 0;
        if (result > std::numeric_limits<std::size_t>::max()
                / static_cast<std::size_t>(chunks)) {
            throw std::overflow_error("maximum partition leaf count exceeds size_t");
        }
        result *= static_cast<std::size_t>(chunks);
    }
    return result;
}

std::uint64_t fnv1a(const std::string& text)
{
    std::uint64_t result = 14695981039346656037ull;
    for (const unsigned char byte : text) {
        result ^= static_cast<std::uint64_t>(byte);
        result *= 1099511628211ull;
    }
    return result;
}

std::vector<BlockLoad> leaf_loads(const std::vector<PartitionLeaf>& leaves)
{
    std::vector<BlockLoad> loads;
    loads.reserve(leaves.size());
    for (const auto& leaf : leaves) {
        loads.push_back({leaf.block, leaf.cell_count()});
    }
    return loads;
}

Real load_ratio(const BlockDistribution& distribution)
{
    const auto& loads = distribution.rank_loads();
    if (loads.empty()) return 1.0;
    const auto maximum = *std::max_element(loads.begin(), loads.end());
    const auto total = std::accumulate(
        loads.begin(), loads.end(), static_cast<std::size_t>(0));
    if (total == 0) return 1.0;
    return static_cast<Real>(maximum)
        / (static_cast<Real>(total) / static_cast<Real>(loads.size()));
}

void assign_stable_ids(std::vector<PartitionLeaf>& leaves)
{
    std::sort(
        leaves.begin(),
        leaves.end(),
        [](const PartitionLeaf& lhs, const PartitionLeaf& rhs) {
            return std::tie(
                lhs.source_zone,
                lhs.cells.begin.i,
                lhs.cells.begin.j,
                lhs.cells.begin.k,
                lhs.cells.end.i,
                lhs.cells.end.j,
                lhs.cells.end.k)
                < std::tie(
                    rhs.source_zone,
                    rhs.cells.begin.i,
                    rhs.cells.begin.j,
                    rhs.cells.begin.k,
                    rhs.cells.end.i,
                    rhs.cells.end.j,
                    rhs.cells.end.k);
        });
    for (std::size_t index = 0; index < leaves.size(); ++index) {
        if (index > static_cast<std::size_t>(std::numeric_limits<BlockId>::max())) {
            throw std::overflow_error("partition leaf id exceeds BlockId range");
        }
        leaves[index].block = static_cast<BlockId>(index);
    }
}

void split_leaf(
    std::vector<PartitionLeaf>& leaves,
    const SplitCandidate& split)
{
    auto left = leaves.at(split.leaf_index);
    auto right = left;
    const auto axis = static_cast<std::size_t>(split.axis);
    left.cells.end[axis] = split.cut;
    right.cells.begin[axis] = split.cut;
    ++left.generation;
    ++right.generation;
    leaves[split.leaf_index] = std::move(left);
    leaves.push_back(std::move(right));
}

SplitCandidate choose_split(
    const std::vector<PartitionLeaf>& leaves,
    int minimum_cells)
{
    std::vector<std::size_t> order(leaves.size());
    std::iota(order.begin(), order.end(), static_cast<std::size_t>(0));
    std::sort(
        order.begin(),
        order.end(),
        [&](std::size_t lhs, std::size_t rhs) {
            if (leaves[lhs].cell_count() != leaves[rhs].cell_count()) {
                return leaves[lhs].cell_count() > leaves[rhs].cell_count();
            }
            return std::tie(
                leaves[lhs].source_zone,
                leaves[lhs].cells.begin.i,
                leaves[lhs].cells.begin.j,
                leaves[lhs].cells.begin.k)
                < std::tie(
                    leaves[rhs].source_zone,
                    leaves[rhs].cells.begin.i,
                    leaves[rhs].cells.begin.j,
                    leaves[rhs].cells.begin.k);
        });
    for (const auto index : order) {
        const auto split = best_split(leaves[index], index, minimum_cells);
        if (split.valid) return split;
    }
    return {};
}

} // namespace

Extent3 CellBox::extent() const
{
    return {end.i - begin.i, end.j - begin.j, end.k - begin.k};
}

std::size_t CellBox::cell_count() const
{
    return extent().size();
}

bool CellBox::valid(int dimension) const
{
    if (dimension != 2 && dimension != 3) return false;
    for (int axis = 0; axis < dimension; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        if (begin[a] < 0 || end[a] <= begin[a]) return false;
    }
    return dimension == 3
        ? end.k > begin.k
        : begin.k == 0 && end.k == 1;
}

Extent3 PartitionLeaf::vertex_extent() const
{
    const auto cell = cell_extent();
    return {
        cell.ni + 1,
        cell.nj + 1,
        cell_dimension == 3 ? cell.nk + 1 : 1,
    };
}

StructuredPartitionPlan StructuredPartitionPlan::build(
    std::vector<PartitionZone> zones,
    int rank_count,
    const PartitionConfig& config)
{
    if (rank_count <= 0) {
        throw CaseConfigurationError("partition rank count must be positive");
    }
    if (zones.empty()) {
        throw CaseConfigurationError("partition requires at least one zone");
    }
    std::sort(
        zones.begin(),
        zones.end(),
        [](const PartitionZone& lhs, const PartitionZone& rhs) {
            return lhs.source_zone < rhs.source_zone;
        });
    for (std::size_t index = 0; index < zones.size(); ++index) {
        const auto& zone = zones[index];
        if (zone.source_zone < 0 || zone.name.empty()
            || (zone.cell_dimension != 2 && zone.cell_dimension != 3)
            || zone.cell_extent.ni <= 0 || zone.cell_extent.nj <= 0
            || (zone.cell_dimension == 3
                ? zone.cell_extent.nk <= 0 : zone.cell_extent.nk != 1)) {
            throw CaseConfigurationError("partition zone metadata are invalid");
        }
        if (index > 0 && zones[index - 1].source_zone == zone.source_zone) {
            throw CaseConfigurationError("partition source zone ids must be unique");
        }
    }

    StructuredPartitionPlan result;
    result.zones_ = std::move(zones);
    for (const auto& zone : result.zones_) {
        const auto feasible = maximum_leaves_for_zone(
            zone, config.min_cells_per_active_direction);
        if (feasible == 0) {
            throw CaseConfigurationError(
                "zone " + zone.name
                + " is smaller than partition/profile minimum");
        }
        if (result.maximum_feasible_leaf_count_
            > std::numeric_limits<std::size_t>::max() - feasible) {
            throw std::overflow_error("maximum partition leaf count exceeds size_t");
        }
        result.maximum_feasible_leaf_count_ += feasible;
        result.leaves_.push_back({
            invalid_block_id,
            zone.source_zone,
            zone.name,
            zone.cell_dimension,
            {{0, 0, 0}, {
                zone.cell_extent.ni,
                zone.cell_extent.nj,
                zone.cell_extent.nk,
            }},
            invalid_rank_id,
            0,
        });
    }

    if (config.mode == PartitionMode::ZonesOnly
        && !config.allow_idle_ranks
        && result.leaves_.size() < static_cast<std::size_t>(rank_count)) {
        throw CaseConfigurationError(
            "zones_only partition has fewer zones than ranks");
    }
    if (!config.allow_idle_ranks
        && result.maximum_feasible_leaf_count_
            < static_cast<std::size_t>(rank_count)) {
        throw CaseConfigurationError(
            "requested ranks exceed maximum feasible leaves: "
            + std::to_string(result.maximum_feasible_leaf_count_));
    }

    const auto required_leaves = config.mode == PartitionMode::ZonesOnly
        ? result.leaves_.size()
        : std::min(
            static_cast<std::size_t>(rank_count),
            result.maximum_feasible_leaf_count_);
    while (result.leaves_.size() < required_leaves) {
        const auto split = choose_split(
            result.leaves_, config.min_cells_per_active_direction);
        if (!split.valid) break;
        split_leaf(result.leaves_, split);
    }

    if (config.mode != PartitionMode::ZonesOnly) {
        const auto leaf_limit = std::min(
            result.maximum_feasible_leaf_count_,
            std::max(
                required_leaves,
                static_cast<std::size_t>(rank_count) * 4));
        while (result.leaves_.size() < leaf_limit) {
            assign_stable_ids(result.leaves_);
            const auto trial = BlockDistribution::balanced(
                leaf_loads(result.leaves_), rank_count);
            if (load_ratio(trial) <= config.max_load_ratio) break;
            const auto split = choose_split(
                result.leaves_, config.min_cells_per_active_direction);
            if (!split.valid) break;
            split_leaf(result.leaves_, split);
        }
    }

    assign_stable_ids(result.leaves_);
    result.distribution_ = BlockDistribution::balanced(
        leaf_loads(result.leaves_), rank_count);
    for (auto& leaf : result.leaves_) {
        leaf.owner = result.distribution_.owner(leaf.block);
    }
    if (!config.allow_idle_ranks) {
        for (const auto load : result.distribution_.rank_loads()) {
            if (load == 0) {
                throw CaseConfigurationError(
                    "partition produced an idle rank despite idle ranks being disabled");
            }
        }
    }
    result.validate();
    result.digest_ = fnv1a(result.summary());
    return result;
}

void StructuredPartitionPlan::validate() const
{
    if (zones_.empty() || leaves_.empty()) {
        throw CaseConfigurationError("partition plan is empty");
    }
    std::unordered_map<BlockId, const PartitionZone*> zones;
    for (const auto& zone : zones_) zones.emplace(zone.source_zone, &zone);
    std::unordered_map<BlockId, std::size_t> covered;
    for (std::size_t first = 0; first < leaves_.size(); ++first) {
        const auto& leaf = leaves_[first];
        if (leaf.block != static_cast<BlockId>(first)) {
            throw CaseConfigurationError("partition leaf ids are not stable and contiguous");
        }
        const auto zone_iterator = zones.find(leaf.source_zone);
        if (zone_iterator == zones.end() || !leaf.cells.valid(leaf.cell_dimension)) {
            throw CaseConfigurationError("partition leaf references invalid zone metadata");
        }
        const auto& zone = *zone_iterator->second;
        for (int axis = 0; axis < zone.cell_dimension; ++axis) {
            const auto a = static_cast<std::size_t>(axis);
            if (leaf.cells.end[a] > zone.cell_extent[a]) {
                throw CaseConfigurationError("partition leaf lies outside source zone");
            }
        }
        covered[leaf.source_zone] += leaf.cell_count();
        for (std::size_t second = first + 1; second < leaves_.size(); ++second) {
            if (leaf.source_zone == leaves_[second].source_zone
                && boxes_overlap(leaf.cells, leaves_[second].cells, leaf.cell_dimension)) {
                throw CaseConfigurationError("partition leaves overlap");
            }
        }
        if (distribution_.owner(leaf.block) != leaf.owner) {
            throw CaseConfigurationError("partition leaf owner disagrees with distribution");
        }
    }
    for (const auto& zone : zones_) {
        if (covered[zone.source_zone] != zone.cell_extent.size()) {
            throw CaseConfigurationError(
                "partition leaves do not exactly cover source zone " + zone.name);
        }
    }
}

std::string StructuredPartitionPlan::summary() const
{
    std::ostringstream result;
    result << "partition_plan(maximum_feasible=" << maximum_feasible_leaf_count_;
    for (const auto& leaf : leaves_) {
        result << ";leaf=" << leaf.block << ",zone=" << leaf.source_zone
               << ",range=[" << leaf.cells.begin.i << ',' << leaf.cells.begin.j
               << ',' << leaf.cells.begin.k << ":" << leaf.cells.end.i << ','
               << leaf.cells.end.j << ',' << leaf.cells.end.k << "),owner="
               << leaf.owner << ",generation=" << leaf.generation;
    }
    result << ')';
    return result.str();
}

} // namespace wcns
