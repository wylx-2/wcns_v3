#pragma once

#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/runtime/case_config.hpp>
#include <wcns/runtime/structured_partition.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>

#include <array>
#include <cstddef>
#include <string>

namespace wcns {

enum class StopReason {
    Running,
    SteadyConverged,
    PhysicalTimeReached,
    MaximumSteps,
    WallTimeCheckpoint,
    UserSignalCheckpoint,
    NumericalFailure,
};

[[nodiscard]] const char* stop_reason_name(StopReason reason);
[[nodiscard]] int stop_reason_exit_code(StopReason reason);

struct ResidualNorms {
    std::array<Real, euler_components> l2 {{}};
    std::array<Real, euler_components> linf {{}};
    bool finite = true;

    [[nodiscard]] Real total_l2() const;
};

[[nodiscard]] ResidualNorms compute_global_residual_norms(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    const AlgorithmProfile& profile);

struct SteadyConvergenceState {
    bool reference_initialized = false;
    std::array<Real, euler_components> reference_l2 {{}};
    std::array<Real, euler_components> reference_linf {{}};
    std::size_t consecutive_passes = 0;
};

struct SimulationProgress {
    std::size_t step = 0;
    Real time = 0.0;
    Real time_step = 0.0;
    Real wall_time = 0.0;
    ResidualNorms residuals;
    bool numerical_failure = false;
    bool user_signal = false;
};

struct StopDecision {
    StopReason reason = StopReason::Running;
    bool residual_checked = false;
};

class StopController {
public:
    explicit StopController(CaseRunConfig config);

    [[nodiscard]] StopDecision evaluate(const SimulationProgress& progress);
    [[nodiscard]] const SteadyConvergenceState& steady_state() const noexcept
    {
        return steady_state_;
    }
    void restore_steady_state(SteadyConvergenceState state);

private:
    [[nodiscard]] bool steady_passed(const ResidualNorms& residuals);

    CaseRunConfig config_;
    SteadyConvergenceState steady_state_;
};

} // namespace wcns
