#include <wcns/parallel/block_distribution.hpp>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace wcns {

BlockDistribution BlockDistribution::balanced(
    std::vector<BlockLoad> loads,
    int rank_count)
{
    if (rank_count <= 0) {
        throw std::invalid_argument("rank count must be positive");
    }
    std::sort(
        loads.begin(),
        loads.end(),
        [](const BlockLoad& lhs, const BlockLoad& rhs) {
            if (lhs.cell_count != rhs.cell_count) {
                return lhs.cell_count > rhs.cell_count;
            }
            return lhs.block < rhs.block;
        });

    BlockDistribution result;
    result.rank_loads_.assign(static_cast<std::size_t>(rank_count), 0);
    result.assignments_.reserve(loads.size());
    result.owners_.reserve(loads.size());
    for (const auto& load : loads) {
        if (load.block < 0 || load.cell_count == 0) {
            throw std::invalid_argument("block loads require a valid id and positive cell count");
        }
        if (result.owners_.find(load.block) != result.owners_.end()) {
            throw std::invalid_argument("block load ids must be unique");
        }
        const auto target = static_cast<RankId>(std::distance(
            result.rank_loads_.begin(),
            std::min_element(result.rank_loads_.begin(), result.rank_loads_.end())));
        auto& rank_load = result.rank_loads_[static_cast<std::size_t>(target)];
        if (load.cell_count > std::numeric_limits<std::size_t>::max() - rank_load) {
            throw std::overflow_error("rank cell load exceeds size_t range");
        }
        rank_load += load.cell_count;
        result.assignments_.push_back({load.block, target, load.cell_count});
        result.owners_.emplace(load.block, target);
    }
    std::sort(
        result.assignments_.begin(),
        result.assignments_.end(),
        [](const BlockAssignment& lhs, const BlockAssignment& rhs) {
            return lhs.block < rhs.block;
        });
    return result;
}

RankId BlockDistribution::owner(BlockId block) const
{
    const auto iterator = owners_.find(block);
    if (iterator == owners_.end()) {
        throw std::out_of_range("block id is absent from the distribution");
    }
    return iterator->second;
}

std::vector<BlockId> BlockDistribution::local_blocks(RankId rank) const
{
    if (rank < 0 || rank >= rank_count()) {
        throw std::out_of_range("rank is outside the distribution");
    }
    std::vector<BlockId> result;
    for (const auto& assignment : assignments_) {
        if (assignment.owner == rank) {
            result.push_back(assignment.block);
        }
    }
    return result;
}

void BlockDistribution::apply(StructuredMesh& mesh) const
{
    if (mesh.block_count() != assignments_.size()) {
        throw TopologyError("mesh and distribution contain different block counts");
    }
    for (const auto& assignment : assignments_) {
        mesh.block(assignment.block).set_owner_rank(assignment.owner);
    }
    for (const auto& descriptor : mesh.blocks()) {
        auto& block = mesh.block(descriptor.id());
        for (auto& connection : block.connectivities) {
            connection.donor_rank = owner(connection.donor_block);
        }
    }
    mesh.validate_connectivities();
}

LocalBlockSet::LocalBlockSet(
    RankId rank,
    std::vector<StructuredBlock> blocks,
    const BlockDistribution& distribution)
    : rank_(rank)
    , blocks_(std::move(blocks))
{
    if (rank < 0 || rank >= distribution.rank_count()) {
        throw std::invalid_argument("local block rank is outside the distribution");
    }
    indices_.reserve(blocks_.size());
    for (std::size_t index = 0; index < blocks_.size(); ++index) {
        auto& block = blocks_[index];
        if (distribution.owner(block.id()) != rank_) {
            throw std::invalid_argument("local block is assigned to another rank");
        }
        block.set_owner_rank(rank_);
        for (auto& connection : block.connectivities) {
            connection.donor_rank = distribution.owner(connection.donor_block);
        }
        if (!indices_.emplace(block.id(), index).second) {
            throw std::invalid_argument("local block ids must be unique");
        }
    }
}

bool LocalBlockSet::contains(BlockId id) const noexcept
{
    return indices_.find(id) != indices_.end();
}

StructuredBlock& LocalBlockSet::block(BlockId id)
{
    const auto iterator = indices_.find(id);
    if (iterator == indices_.end()) {
        throw std::out_of_range("block id is not local to this rank");
    }
    return blocks_[iterator->second];
}

const StructuredBlock& LocalBlockSet::block(BlockId id) const
{
    const auto iterator = indices_.find(id);
    if (iterator == indices_.end()) {
        throw std::out_of_range("block id is not local to this rank");
    }
    return blocks_[iterator->second];
}

} // namespace wcns

