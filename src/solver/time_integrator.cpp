#include <wcns/solver/time_integrator.hpp>

#include <wcns/solver/flow_fields.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace wcns {
namespace {

using StateBuffer = std::vector<Real>;

StateBuffer capture_interior(const StructuredBlock& block)
{
    const auto extent = block.cell_extent();
    StateBuffer result;
    result.reserve(extent.size() * static_cast<std::size_t>(euler_components));
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < euler_components; ++component) {
                    result.push_back(block.flow.conservative(i, j, k, component));
                }
            }
        }
    }
    return result;
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

void advance_ssprk3(
    const std::vector<StructuredBlock*>& blocks,
    Real time_step,
    const ResidualEvaluator& evaluate_residuals)
{
    if (blocks.empty() || !evaluate_residuals) {
        throw std::invalid_argument("SSPRK3 requires blocks and a residual evaluator");
    }
    if (!std::isfinite(time_step) || time_step <= 0.0) {
        throw std::invalid_argument("SSPRK3 time step must be positive and finite");
    }
    std::vector<StateBuffer> initial;
    initial.reserve(blocks.size());
    for (const auto* block : blocks) {
        if (block == nullptr) {
            throw std::invalid_argument("SSPRK3 block pointer must not be null");
        }
        initial.push_back(capture_interior(*block));
    }

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(*blocks[b], initial[b], 1.0, 0.0, time_step);
    }

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(*blocks[b], initial[b], 0.75, 0.25, 0.25 * time_step);
    }

    evaluate_residuals();
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        update_stage(
            *blocks[b], initial[b], 1.0 / 3.0, 2.0 / 3.0, 2.0 * time_step / 3.0);
    }
}

} // namespace wcns
