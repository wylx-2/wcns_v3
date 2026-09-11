#include <wcns/solver/time_integrator.hpp>

#include <wcns/solver/flow_fields.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace wcns {
namespace {

using StateBuffer = std::vector<Real>;

void capture_interior(const StructuredBlock& block, StateBuffer& result)
{
    const auto extent = block.cell_extent();
    result.resize(extent.size() * static_cast<std::size_t>(euler_components));
    std::size_t offset = 0;
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < euler_components; ++component) {
                    result[offset++]
                        = block.flow.conservative(i, j, k, component);
                }
            }
        }
    }
}

void update_stage(
    StructuredBlock& block,
    const StateBuffer& initial,
    Real initial_weight,
    Real stage_weight,
    Real residual_weight)
{
    const auto extent = block.cell_extent();
    std::size_t offset = 0;
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < euler_components; ++component) {
                    auto& value = block.flow.conservative(i, j, k, component);
                    value = initial_weight * initial[offset] + stage_weight * value
                        + residual_weight * block.flow.residual(i, j, k, component);
                    ++offset;
                }
            }
        }
    }
}

} // namespace

void SsprkWorkspace::prepare(
    const std::vector<StructuredBlock*>& blocks)
{
    initial_.resize(blocks.size());
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        if (blocks[b] == nullptr) {
            throw std::invalid_argument("SSPRK3 block pointer must not be null");
        }
        capture_interior(*blocks[b], initial_[b]);
    }
}

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals)
{
    SsprkWorkspace workspace;
    advance_ssprk3(blocks, workspace, time_step, evaluate_residuals);
}

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    SsprkWorkspace& workspace,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals)
{
    if (blocks.empty() || !evaluate_residuals) {
        throw std::invalid_argument("SSPRK3 requires blocks and a residual evaluator");
    }
    if (!std::isfinite(time_step) || time_step <= 0.0) {
        throw std::invalid_argument("SSPRK3 time step must be positive and finite");
    }
    workspace.prepare(blocks);

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(*blocks[b], workspace.initial_[b], 1.0, 0.0, time_step);
    }

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(
            *blocks[b], workspace.initial_[b], 0.75, 0.25, 0.25 * time_step);
    }

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(
            *blocks[b], workspace.initial_[b], 1.0 / 3.0, 2.0 / 3.0,
            2.0 * time_step / 3.0);
    }
}

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    Real initial_time,
    const TimedResidualEvaluator& evaluate_residuals)
{
    SsprkWorkspace workspace;
    advance_ssprk3(
        blocks, workspace, time_step, initial_time, evaluate_residuals);
}

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    SsprkWorkspace& workspace,
    Real time_step,
    Real initial_time,
    const TimedResidualEvaluator& evaluate_residuals)
{
    if (blocks.empty() || !evaluate_residuals) {
        throw std::invalid_argument("timed SSPRK3 requires blocks and a residual evaluator");
    }
    if (!std::isfinite(time_step) || time_step <= 0.0
        || !std::isfinite(initial_time)) {
        throw std::invalid_argument("timed SSPRK3 time inputs are invalid");
    }
    workspace.prepare(blocks);

    evaluate_residuals(initial_time);
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(*blocks[b], workspace.initial_[b], 1.0, 0.0, time_step);
    }
    evaluate_residuals(initial_time + time_step);
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(
            *blocks[b], workspace.initial_[b], 0.75, 0.25, 0.25 * time_step);
    }
    evaluate_residuals(initial_time + 0.5 * time_step);
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(
            *blocks[b], workspace.initial_[b], 1.0 / 3.0, 2.0 / 3.0,
            2.0 * time_step / 3.0);
    }
}

} // namespace wcns
