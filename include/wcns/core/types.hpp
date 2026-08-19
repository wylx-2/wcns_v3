#pragma once

#include <cstdint>

namespace wcns {

using Real = double;
using BlockId = std::int32_t;
using RankId = std::int32_t;

inline constexpr BlockId invalid_block_id = static_cast<BlockId>(-1);
inline constexpr RankId invalid_rank_id = static_cast<RankId>(-1);

} // namespace wcns

