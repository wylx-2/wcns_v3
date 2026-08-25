#pragma once

#include <wcns/mesh/structured_block.hpp>

#include <functional>
#include <vector>

namespace wcns {

using ResidualEvaluator = std::function<void()>;

// Advances all supplied blocks synchronously. The evaluator must refresh
// primitive/halo/boundary state and then compute every block residual.
void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals);

} // namespace wcns
