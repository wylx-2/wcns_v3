#pragma once

#include <wcns/core/types.hpp>
#include <wcns/mesh/structured_block.hpp>

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace wcns {

class StructuredMesh {
public:
    StructuredMesh() = default;
    explicit StructuredMesh(std::vector<StructuredBlock> blocks);

    [[nodiscard]] std::size_t block_count() const noexcept { return blocks_.size(); }
    [[nodiscard]] bool contains(BlockId id) const noexcept;

    [[nodiscard]] StructuredBlock& block(BlockId id);
    [[nodiscard]] const StructuredBlock& block(BlockId id) const;

    [[nodiscard]] const std::vector<StructuredBlock>& blocks() const noexcept
    {
        return blocks_;
    }

    void validate_connectivities(bool validate_coordinates = true) const;

private:
    std::vector<StructuredBlock> blocks_;
    std::unordered_map<BlockId, std::size_t> block_index_;
};

} // namespace wcns
