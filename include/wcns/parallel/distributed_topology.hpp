#pragma once

#include <wcns/mesh/halo_exchange.hpp>
#include <wcns/mesh/structured_mesh.hpp>
#include <wcns/parallel/block_distribution.hpp>

#include <vector>

namespace wcns {

struct DirectedExchange {
    ConnectionId connection = invalid_connection_id;
    RankId receiver_rank = invalid_rank_id;
    RankId donor_rank = invalid_rank_id;
    HaloExchangePlan halo;

    [[nodiscard]] int message_tag(int tag_base = 1024) const;
};

class DistributedTopology {
public:
    [[nodiscard]] static DistributedTopology build(
        const StructuredMesh& mesh,
        const BlockDistribution& distribution);

    [[nodiscard]] const std::vector<DirectedExchange>& exchanges() const noexcept
    {
        return exchanges_;
    }
    [[nodiscard]] std::vector<const DirectedExchange*> receives(RankId rank) const;
    [[nodiscard]] std::vector<const DirectedExchange*> sends(RankId rank) const;
    [[nodiscard]] std::vector<const DirectedExchange*> local_copies(RankId rank) const;

private:
    std::vector<DirectedExchange> exchanges_;
};

} // namespace wcns

