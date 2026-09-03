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

void CompositeSimulationObserver::add(ISimulationObserver& observer)
{
    if (std::find(observers_.begin(), observers_.end(), &observer)
        != observers_.end()) {
        throw std::invalid_argument("simulation observer was added twice");
    }
    observers_.push_back(&observer);
}

Real CompositeSimulationObserver::next_time_event(
    const SimulationState& state) const
{
    Real result = std::numeric_limits<Real>::infinity();
    for (const auto* observer : observers_) {
        result = std::min(result, observer->next_time_event(state));
    }
    return result;
}

void CompositeSimulationObserver::on_initial(const SimulationState& state)
{
    for (auto* observer : observers_) observer->on_initial(state);
}

void CompositeSimulationObserver::on_step(
    const SimulationState& state,
    bool residual_checked)
{
    for (auto* observer : observers_) {
        observer->on_step(state, residual_checked);
    }
}

void CompositeSimulationObserver::on_final(const SimulationState& state)
{
    for (auto* observer : observers_) observer->on_final(state);
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
    const auto notify = [&](const std::function<void()>& callback) {
        bool local_observer_success = true;
        std::string local_error;
        try {
            callback();
        } catch (const std::exception& error) {
            local_observer_success = false;
            local_error = error.what();
        } catch (...) {
            local_observer_success = false;
            local_error = "unknown non-standard exception";
        }
        if (!mpi_.all_true(local_observer_success)) {
            const Real failure_candidate = local_observer_success
                ? static_cast<Real>(mpi_.size())
                : static_cast<Real>(mpi_.rank());
            const auto failure_rank = static_cast<RankId>(
                mpi_.min(failure_candidate));
            const std::string collective_error = mpi_.broadcast_string(
                mpi_.rank() == failure_rank ? std::move(local_error) : std::string {},
                failure_rank);
            throw std::runtime_error(
                "simulation output/observer failed on MPI rank "
                + std::to_string(failure_rank) + ": " + collective_error);
        }
    };
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
        notify([&] { observer_.on_initial(state); });
    } else {
        state.residuals.finite = false;
        state.stop_reason = StopReason::NumericalFailure;
        state.steady = controller.steady_state();
        notify([&] { observer_.on_final(state); });
        return state;
    }

    if (config_.mode == RunMode::Unsteady) {
        const Real tolerance = 64.0 * std::numeric_limits<Real>::epsilon()
            * std::max(Real {1.0}, config_.end_time);
        if (state.time + tolerance >= config_.end_time) {
            state.stop_reason = StopReason::PhysicalTimeReached;
            notify([&] { observer_.on_final(state); });
            return state;
        }
    } else if (state.steady.reference_initialized
        && state.step >= config_.steady.min_steps
        && state.steady.consecutive_passes
            >= config_.steady.consecutive_checks) {
        state.stop_reason = StopReason::SteadyConverged;
        notify([&] { observer_.on_final(state); });
        return state;
    }
    if (state.step >= config_.max_steps) {
        state.stop_reason = StopReason::MaximumSteps;
        notify([&] { observer_.on_final(state); });
        return state;
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
        notify([&] { observer_.on_step(state, decision.residual_checked); });
    }
    notify([&] { observer_.on_final(state); });
    return state;
}

} // namespace wcns
