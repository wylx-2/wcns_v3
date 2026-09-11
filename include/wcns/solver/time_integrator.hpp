#pragma once

#include <wcns/mesh/structured_block.hpp>

#include <functional>
#include <vector>

namespace wcns {

using ResidualEvaluator = std::function<void()>;
using TimedResidualEvaluator = std::function<void(Real)>;

class SsprkWorkspace {
public:
    void prepare(const std::vector<StructuredBlock*>& blocks);

private:
    std::vector<std::vector<Real>> initial_;

    friend void advance_ssprk3(
        const std::vector<StructuredBlock*>&, SsprkWorkspace&, Real,
        const ResidualEvaluator&);
    friend void advance_ssprk3(
        const std::vector<StructuredBlock*>&, SsprkWorkspace&, Real, Real,
        const TimedResidualEvaluator&);
};

// Advances all supplied blocks synchronously. The evaluator must refresh
// primitive/halo/boundary state and then compute every block residual.
void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals);

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    SsprkWorkspace& workspace,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals);

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    Real initial_time,
    const TimedResidualEvaluator& evaluate_residuals);

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    SsprkWorkspace& workspace,
    Real time_step,
    Real initial_time,
    const TimedResidualEvaluator& evaluate_residuals);

} // namespace wcns
