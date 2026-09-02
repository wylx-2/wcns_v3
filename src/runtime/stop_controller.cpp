#include <wcns/runtime/stop_controller.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace wcns {
namespace {

using LineWeights = std::array<LineConservationWeights, 3>;

const PartitionLeaf& leaf_for(
    const std::unordered_map<BlockId, const PartitionLeaf*>& leaves,
    BlockId block)
{
    const auto iterator = leaves.find(block);
    if (iterator == leaves.end()) {
        throw CaseConfigurationError(
            "residual norm cannot find partition leaf for block "
            + std::to_string(block));
    }
    return *iterator->second;
}

} // namespace

const char* stop_reason_name(StopReason reason)
{
    switch (reason) {
    case StopReason::Running: return "running";
    case StopReason::SteadyConverged: return "steady_converged";
    case StopReason::PhysicalTimeReached: return "physical_time_reached";
    case StopReason::MaximumSteps: return "maximum_steps";
    case StopReason::WallTimeCheckpoint: return "wall_time_checkpoint";
    case StopReason::UserSignalCheckpoint: return "user_signal_checkpoint";
    case StopReason::NumericalFailure: return "numerical_failure";
    }
    throw std::invalid_argument("invalid stop reason");
}

int stop_reason_exit_code(StopReason reason)
{
    switch (reason) {
    case StopReason::SteadyConverged:
    case StopReason::PhysicalTimeReached:
        return 0;
    case StopReason::MaximumSteps:
    case StopReason::WallTimeCheckpoint:
    case StopReason::UserSignalCheckpoint:
        return 2;
    case StopReason::NumericalFailure:
        return 3;
    case StopReason::Running:
        return 4;
    }
    return 4;
}

Real ResidualNorms::total_l2() const
{
    Real sum = 0.0;
    for (const Real value : l2) sum += value * value;
    return std::sqrt(sum / static_cast<Real>(euler_components));
}

ResidualNorms compute_global_residual_norms(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    const AlgorithmProfile& profile)
{
    std::unordered_map<BlockId, const PartitionLeaf*> leaves;
    for (const auto& leaf : partition.leaves()) {
        leaves.emplace(leaf.block, &leaf);
    }
    std::unordered_map<BlockId, LineWeights> zone_weights;
    for (const auto& zone : partition.zones()) {
        LineWeights lines;
        lines[0] = build_line_conservation_weights(
            profile, zone.cell_extent.ni);
        lines[1] = build_line_conservation_weights(
            profile, zone.cell_extent.nj);
        lines[2] = build_line_conservation_weights(
            profile,
            zone.cell_dimension == 3 ? zone.cell_extent.nk : 1,
            zone.cell_dimension == 2);
        zone_weights.emplace(zone.source_zone, std::move(lines));
    }

    std::array<Real, euler_components> local_squared {{}};
    std::array<Real, euler_components> local_maximum {{}};
    Real local_weight = 0.0;
    bool local_finite = true;
    for (const auto& block : local_blocks.blocks()) {
        const auto& leaf = leaf_for(leaves, block.id());
        const auto metric_iterator = metrics.find(block.id());
        const auto line_iterator = zone_weights.find(leaf.source_zone);
        if (metric_iterator == metrics.end() || line_iterator == zone_weights.end()) {
            throw CaseConfigurationError(
                "residual norm is missing metric or source-zone weights");
        }
        if (block.cell_extent() != leaf.cell_extent()) {
            throw CaseConfigurationError(
                "residual norm block and partition extents differ");
        }
        const auto& jacobian = metric_iterator->second.jacobian();
        const auto& lines = line_iterator->second;
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real weight
                        = lines[0].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.i + i)]
                        * lines[1].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.j + j)]
                        * lines[2].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.k + k)]
                        * jacobian(i, j, k);
                    if (!std::isfinite(weight) || weight <= 0.0) {
                        local_finite = false;
                        continue;
                    }
                    local_weight += weight;
                    for (int component = 0; component < euler_components;
                         ++component) {
                        const Real value = block.flow.residual(
                            i, j, k, component);
                        if (!std::isfinite(value)) {
                            local_finite = false;
                            continue;
                        }
                        const auto index = static_cast<std::size_t>(component);
                        local_squared[index] += weight * value * value;
                        local_maximum[index] = std::max(
                            local_maximum[index], std::abs(value));
                    }
                }
            }
        }
    }

    ResidualNorms result;
    const Real global_weight = mpi.sum(local_weight);
    result.finite = mpi.all_true(local_finite)
        && std::isfinite(global_weight) && global_weight > 0.0;
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        const Real squared = mpi.sum(local_squared[index]);
        result.linf[index] = mpi.max(local_maximum[index]);
        if (!std::isfinite(squared) || squared < 0.0
            || !std::isfinite(result.linf[index])) {
            result.finite = false;
        }
        result.l2[index] = result.finite
            ? std::sqrt(squared / global_weight)
            : std::numeric_limits<Real>::quiet_NaN();
    }
    if (!result.finite) {
        result.l2.fill(std::numeric_limits<Real>::quiet_NaN());
        result.linf.fill(std::numeric_limits<Real>::quiet_NaN());
    }
    return result;
}

StopController::StopController(CaseRunConfig config)
    : config_(std::move(config))
{
    config_.validate();
}

bool StopController::steady_passed(const ResidualNorms& residuals)
{
    if (!steady_state_.reference_initialized) {
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            steady_state_.reference_l2[index] = std::max(
                residuals.l2[index], config_.steady.reference_floor);
            steady_state_.reference_linf[index] = std::max(
                residuals.linf[index], config_.steady.reference_floor);
        }
        steady_state_.reference_initialized = true;
    }
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        const bool l2 = residuals.l2[index] <= config_.steady.l2_absolute
            || residuals.l2[index] / steady_state_.reference_l2[index]
                <= config_.steady.l2_relative;
        const bool linf = !config_.steady.linf_enabled
            || residuals.linf[index] <= config_.steady.linf_absolute
            || residuals.linf[index] / steady_state_.reference_linf[index]
                <= config_.steady.linf_relative;
        if (!l2 || !linf) return false;
    }
    return true;
}

StopDecision StopController::evaluate(const SimulationProgress& progress)
{
    StopDecision decision;
    if (progress.numerical_failure || !progress.residuals.finite
        || !std::isfinite(progress.time) || progress.time < 0.0
        || !std::isfinite(progress.time_step) || progress.time_step <= 0.0) {
        decision.reason = StopReason::NumericalFailure;
        return decision;
    }

    if (config_.mode == RunMode::Steady
        && progress.step % config_.steady.check_interval_steps == 0) {
        decision.residual_checked = true;
        const bool passed = steady_passed(progress.residuals);
        if (progress.step >= config_.steady.min_steps && passed) {
            ++steady_state_.consecutive_passes;
        } else {
            steady_state_.consecutive_passes = 0;
        }
        if (steady_state_.consecutive_passes
            >= config_.steady.consecutive_checks) {
            decision.reason = StopReason::SteadyConverged;
            return decision;
        }
    }

    if (config_.mode == RunMode::Unsteady) {
        const Real tolerance = 64.0 * std::numeric_limits<Real>::epsilon()
            * std::max(Real {1.0}, config_.end_time);
        if (progress.time + tolerance >= config_.end_time) {
            decision.reason = StopReason::PhysicalTimeReached;
            return decision;
        }
    }
    if (progress.user_signal) {
        decision.reason = StopReason::UserSignalCheckpoint;
        return decision;
    }
    if (config_.max_wall_time > 0.0
        && progress.wall_time >= config_.max_wall_time) {
        decision.reason = StopReason::WallTimeCheckpoint;
        return decision;
    }
    if (progress.step >= config_.max_steps) {
        decision.reason = StopReason::MaximumSteps;
    }
    return decision;
}

void StopController::restore_steady_state(SteadyConvergenceState state)
{
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        if (state.reference_initialized
            && (!std::isfinite(state.reference_l2[index])
                || state.reference_l2[index] <= 0.0
                || !std::isfinite(state.reference_linf[index])
                || state.reference_linf[index] <= 0.0)) {
            throw CaseConfigurationError(
                "restored steady residual reference is invalid");
        }
    }
    steady_state_ = std::move(state);
}

} // namespace wcns
