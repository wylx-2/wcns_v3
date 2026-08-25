#pragma once

#include <wcns/core/types.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/mesh/structured_mesh.hpp>

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace wcns {

struct BlockLoad {
    BlockId block = invalid_block_id;
    std::size_t cell_count = 0;
};

struct BlockAssignment {
    BlockId block = invalid_block_id;
    RankId owner = invalid_rank_id;
    std::size_t cell_count = 0;
};

class BlockDistribution {
public:
    [[nodiscard]] static BlockDistribution balanced(
        std::vector<BlockLoad> loads,
        int rank_count);

    [[nodiscard]] RankId owner(BlockId block) const;
    [[nodiscard]] int rank_count() const noexcept
    {
        return static_cast<int>(rank_loads_.size());
    }
    [[nodiscard]] const std::vector<std::size_t>& rank_loads() const noexcept
    {
        return rank_loads_;
    }
    [[nodiscard]] const std::vector<BlockAssignment>& assignments() const noexcept
    {
        return assignments_;
    }
    [[nodiscard]] std::vector<BlockId> local_blocks(RankId rank) const;

    void apply(StructuredMesh& mesh) const;

private:
    std::vector<BlockAssignment> assignments_;
    std::vector<std::size_t> rank_loads_;
    std::unordered_map<BlockId, RankId> owners_;
};

class LocalBlockSet {
public:
    LocalBlockSet(
        RankId rank,
        std::vector<StructuredBlock> blocks,
        const BlockDistribution& distribution);

    [[nodiscard]] RankId rank() const noexcept { return rank_; }
    [[nodiscard]] bool contains(BlockId id) const noexcept;
    [[nodiscard]] StructuredBlock& block(BlockId id);
    [[nodiscard]] const StructuredBlock& block(BlockId id) const;
    [[nodiscard]] std::vector<StructuredBlock>& blocks() noexcept { return blocks_; }
    [[nodiscard]] const std::vector<StructuredBlock>& blocks() const noexcept
    {
        return blocks_;
    }

private:
    RankId rank_ = invalid_rank_id;
    std::vector<StructuredBlock> blocks_;
    std::unordered_map<BlockId, std::size_t> indices_;
};

} // namespace wcns

