#include <wcns/runtime/simulation_driver.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace wcns {
namespace {

Real default_wall_clock()
{
    using Clock = std::chrono::steady_clock;
    static const auto origin = Clock::now();
    return std::chrono::duration<Real>(Clock::now() - origin).count();
}

template <class Solver>
SolverDiagnostics collect_diagnostics(const Solver& solver)
{
    return {
        solver.global_reconstruction_fallback_count(),
        solver.global_riemann_face_count(),
        solver.global_riemann_fallback_count(),
    };
}

} // namespace

InviscidSimulationSolver::InviscidSimulationSolver(
    InviscidWcnsSolver& solver,
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    AlgorithmProfile profile)
    : solver_(solver)
    , mpi_(mpi)
    , local_blocks_(local_blocks)
    , metrics_(metrics)
    , partition_(partition)
    , profile_(std::move(profile))
{
}

Real InviscidSimulationSolver::global_time_step(Real cfl)
{
    return solver_.global_time_step(cfl);
}

void InviscidSimulationSolver::advance(Real time_step, Real initial_time)
{
    solver_.advance(time_step, initial_time);
}

void InviscidSimulationSolver::refresh_residuals(Real time)
{
    solver_.compute_residuals(time);
}

ResidualNorms InviscidSimulationSolver::residual_norms() const
{
    return compute_global_residual_norms(
        mpi_, local_blocks_, metrics_, partition_, profile_);
}

SolverDiagnostics InviscidSimulationSolver::diagnostics() const
{
    return collect_diagnostics(solver_);
}

ViscousSimulationSolver::ViscousSimulationSolver(
    ViscousWcnsSolver& solver,
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    AlgorithmProfile profile)
    : solver_(solver)
    , mpi_(mpi)
    , local_blocks_(local_blocks)
    , metrics_(metrics)
    , partition_(partition)
    , profile_(std::move(profile))
{
}

Real ViscousSimulationSolver::global_time_step(Real cfl)
{
    return solver_.global_time_step(cfl);
}

void ViscousSimulationSolver::advance(Real time_step, Real initial_time)
{
    solver_.advance(time_step, initial_time);
}

void ViscousSimulationSolver::refresh_residuals(Real time)
{
    solver_.compute_residuals(time);
}

ResidualNorms ViscousSimulationSolver::residual_norms() const
{
    return compute_global_residual_norms(
        mpi_, local_blocks_, metrics_, partition_, profile_);
}

SolverDiagnostics ViscousSimulationSolver::diagnostics() const
{
    return collect_diagnostics(solver_);
}

SimulationDriver::SimulationDriver(
    const MpiRuntime& mpi,
    ISimulationSolver& solver,
    CaseRunConfig config,
    ISimulationObserver& observer,
    StopRequest stop_requested,
    WallClock wall_clock)
    : mpi_(mpi)
    , solver_(solver)
    , config_(std::move(config))
    , observer_(observer)
    , stop_requested_(std::move(stop_requested))
    , wall_clock_(std::move(wall_clock))
{
    config_.validate();
    if (!stop_requested_) stop_requested_ = [] { return false; };
    if (!wall_clock_) wall_clock_ = default_wall_clock;
}

bool SimulationDriver::all_ranks_succeeded(bool local_success) const
{
    return mpi_.all_true(local_success);
}

Real SimulationDriver::limited_time_step(
    Real proposed,
    const SimulationState& state) const
{
    Real result = proposed;
    if (config_.mode == RunMode::Unsteady) {
        result = std::min(result, config_.end_time - state.time);
    }
    const Real event = observer_.next_time_event(state);
    if (std::isfinite(event)) {
        const Real tolerance = 64.0 * std::numeric_limits<Real>::epsilon()
            * std::max({Real {1.0}, std::abs(state.time), std::abs(event)});
        if (event > state.time + tolerance) {
            result = std::min(result, event - state.time);
        }
    }
    return result;
}

SimulationState SimulationDriver::run(SimulationInitialState initial)
{
    if (!std::isfinite(initial.time) || initial.time < 0.0) {
        throw CaseConfigurationError(
            "simulation initial time must be finite and non-negative");
    }
    StopController controller(config_);
    controller.restore_steady_state(std::move(initial.steady));
    SimulationState state;
    state.step = initial.step;
    state.time = initial.time;
    const Real wall_origin = wall_clock_();

    bool local_success = true;
    try {
        solver_.refresh_residuals(state.time);
    } catch (const std::exception&) {
        local_success = false;
    }
    const bool initial_success = all_ranks_succeeded(local_success);
    if (initial_success) {
        state.residuals = solver_.residual_norms();
        state.diagnostics = solver_.diagnostics();
        observer_.on_initial(state);
    } else {
        state.residuals.finite = false;
        state.stop_reason = StopReason::NumericalFailure;
        state.steady = controller.steady_state();
        observer_.on_final(state);
        return state;
    }

    if (config_.mode == RunMode::Unsteady) {
        const Real tolerance = 64.0 * std::numeric_limits<Real>::epsilon()
            * std::max(Real {1.0}, config_.end_time);
        if (state.time + tolerance >= config_.end_time) {
            state.stop_reason = StopReason::PhysicalTimeReached;
            observer_.on_final(state);
            return state;
        }
    }

    while (state.stop_reason == StopReason::Running) {
        local_success = true;
        Real proposed = std::numeric_limits<Real>::quiet_NaN();
        try {
            proposed = solver_.global_time_step(config_.cfl);
            state.time_step = limited_time_step(proposed, state);
            if (!std::isfinite(state.time_step) || state.time_step <= 0.0) {
                local_success = false;
            }
        } catch (const std::exception&) {
            local_success = false;
        }
        if (all_ranks_succeeded(local_success)) {
            try {
                solver_.advance(state.time_step, state.time);
            } catch (const std::exception&) {
                local_success = false;
            }
        }
        const bool advance_success = all_ranks_succeeded(local_success);
        if (advance_success) {
            ++state.step;
            state.time += state.time_step;
            try {
                solver_.refresh_residuals(state.time);
            } catch (const std::exception&) {
                local_success = false;
            }
        }
        const bool residual_success = all_ranks_succeeded(
            advance_success && local_success);
        if (residual_success) {
            state.residuals = solver_.residual_norms();
            state.diagnostics = solver_.diagnostics();
        } else {
            state.residuals.finite = false;
        }
        state.wall_time = wall_clock_() - wall_origin;
        const bool local_signal = stop_requested_();
        const bool global_signal = !mpi_.all_true(!local_signal);
        const SimulationProgress progress {
            state.step,
            state.time,
            state.time_step,
            state.wall_time,
            state.residuals,
            !residual_success,
            global_signal,
        };
        const auto decision = controller.evaluate(progress);
        state.stop_reason = decision.reason;
        state.steady = controller.steady_state();
        observer_.on_step(state, decision.residual_checked);
    }
    observer_.on_final(state);
    return state;
}

} // namespace wcns
