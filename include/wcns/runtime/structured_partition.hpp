#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/types.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/runtime/case_config.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wcns {

struct CellBox {
    Index3 begin {};
    Index3 end {};

    [[nodiscard]] Extent3 extent() const;
    [[nodiscard]] std::size_t cell_count() const;
    [[nodiscard]] bool valid(int dimension) const;
};

struct PartitionZone {
    BlockId source_zone = invalid_block_id;
    std::string name;
    int cell_dimension = 0;
    Extent3 cell_extent {};
};

struct PartitionLeaf {
    BlockId block = invalid_block_id;
    BlockId source_zone = invalid_block_id;
    std::string source_zone_name;
    int cell_dimension = 0;
    CellBox cells {};
    RankId owner = invalid_rank_id;
    int generation = 0;

    [[nodiscard]] Extent3 cell_extent() const { return cells.extent(); }
    [[nodiscard]] Extent3 vertex_extent() const;
    [[nodiscard]] std::size_t cell_count() const { return cells.cell_count(); }
};

class StructuredPartitionPlan {
public:
    [[nodiscard]] static StructuredPartitionPlan build(
        std::vector<PartitionZone> zones,
        int rank_count,
        const PartitionConfig& config);

    [[nodiscard]] const std::vector<PartitionZone>& zones() const noexcept
    {
        return zones_;
    }
    [[nodiscard]] const std::vector<PartitionLeaf>& leaves() const noexcept
    {
        return leaves_;
    }
    [[nodiscard]] const BlockDistribution& distribution() const noexcept
    {
        return distribution_;
    }
    [[nodiscard]] std::size_t maximum_feasible_leaf_count() const noexcept
    {
        return maximum_feasible_leaf_count_;
    }
    [[nodiscard]] std::uint64_t digest() const noexcept { return digest_; }
    [[nodiscard]] std::string summary() const;

    void validate() const;

private:
    std::vector<PartitionZone> zones_;
    std::vector<PartitionLeaf> leaves_;
    BlockDistribution distribution_;
    std::size_t maximum_feasible_leaf_count_ = 0;
    std::uint64_t digest_ = 0;
};

} // namespace wcns
