#pragma once

#include <wcns/mesh/topology.hpp>

#include <cstdint>

namespace wcns {

// Identifies the owned face and residual evaluation that produced a local event.
// rk_stage is zero for a direct residual evaluation and 1--3 inside SSPRK3.
struct FaceDiagnosticLocation {
    BlockId block = invalid_block_id;
    RankId rank = invalid_rank_id;
    Axis axis = Axis::I;
    Index3 face {};
    std::uint64_t residual_evaluation = 0;
    int rk_stage = 0;
};

} // namespace wcns
