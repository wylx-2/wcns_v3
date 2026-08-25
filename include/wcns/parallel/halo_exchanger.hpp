#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>

#include <unordered_map>

namespace wcns {

class BlockFieldRegistry {
public:
    explicit BlockFieldRegistry(int components);

    void add(BlockId block, Field<Real>& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] Field<Real>& field(BlockId block) const;
    [[nodiscard]] int components() const noexcept { return components_; }

private:
    int components_ = 0;
    std::unordered_map<BlockId, Field<Real>*> fields_;
};

class HaloExchanger {
public:
    HaloExchanger(
        const MpiRuntime& mpi,
        const DistributedTopology& topology,
        int distribution_rank_count);

    void exchange(const BlockFieldRegistry& fields) const;

private:
    const MpiRuntime& mpi_;
    const DistributedTopology& topology_;
};

} // namespace wcns

