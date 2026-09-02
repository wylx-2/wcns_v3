#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/types.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/mesh/structured_mesh.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace wcns {

class CgnsError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct CgnsBaseMetadata {
    int file_index = 0;
    std::string name;
    int cell_dimension = 0;
    int physical_dimension = 0;
};

struct CgnsZoneMetadata {
    BlockId block_id = invalid_block_id;
    int base_file_index = 0;
    int zone_file_index = 0;
    std::string base_name;
    std::string name;
    int cell_dimension = 0;
    int physical_dimension = 0;
    Extent3 vertex_extent;
    Extent3 cell_extent;
};

struct CgnsMeshMetadata {
    std::vector<CgnsBaseMetadata> bases;
    std::vector<CgnsZoneMetadata> zones;
};

struct CgnsPartitionLeaf {
    BlockId block_id = invalid_block_id;
    BlockId source_zone = invalid_block_id;
    Index3 cell_begin {};
    Index3 cell_end {};
    RankId owner_rank = invalid_rank_id;
};

struct CgnsPartitionedMesh {
    StructuredMesh global_mesh;
    std::vector<StructuredBlock> local_blocks;
};

class CgnsReader {
public:
    [[nodiscard]] CgnsMeshMetadata read_metadata(const std::string& path) const;

    [[nodiscard]] StructuredBlock read_block(
        const std::string& path,
        const CgnsZoneMetadata& zone,
        RankId owner_rank,
        int ghost_width) const;

    [[nodiscard]] StructuredMesh read_mesh(
        const std::string& path,
        RankId owner_rank,
        int ghost_width) const;

    [[nodiscard]] CgnsPartitionedMesh read_partitioned_mesh(
        const std::string& path,
        const std::vector<CgnsPartitionLeaf>& leaves,
        RankId local_rank,
        int ghost_width) const;
};

} // namespace wcns
