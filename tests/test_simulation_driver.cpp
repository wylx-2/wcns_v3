#include "test_support.hpp"

#include <wcns/runtime/simulation_driver.hpp>

#include <cmath>
#include <vector>

namespace {

class FakeSolver final : public wcns::ISimulationSolver {
public:
    wcns::Real proposed_time_step = 0.3;
    wcns::Real residual_value = 1.0;
    std::size_t advance_count = 0;
    std::vector<wcns::Real> time_steps;

    wcns::Real global_time_step(wcns::Real) override
    {
        return proposed_time_step;
    }

    void advance(wcns::Real time_step, wcns::Real) override
    {
        ++advance_count;
        time_steps.push_back(time_step);
        residual_value *= 1.0e-4;
    }

    void refresh_residuals(wcns::Real) override {}

    wcns::ResidualNorms residual_norms() const override
    {
        wcns::ResidualNorms result;
        result.l2.fill(residual_value);
        result.linf.fill(2.0 * residual_value);
        return result;
    }

    wcns::SolverDiagnostics diagnostics() const override
    {
        return {};
    }
};

class TimeEventObserver final : public wcns::ISimulationObserver {
public:
    wcns::Real next_time_event(const wcns::SimulationState& state) const override
    {
        return state.time < 0.5 - 1.0e-14
            ? 0.5 : std::numeric_limits<wcns::Real>::infinity();
    }

    void on_step(const wcns::SimulationState& state, bool) override
    {
        times.push_back(state.time);
    }

    void on_final(const wcns::SimulationState&) override
    {
        ++final_count;
    }

    std::vector<wcns::Real> times;
    int final_count = 0;
};

} // namespace

// 验收唯一驱动精确命中时间事件/t_end，并在最终时刻与 max_steps 重合时报告正常完成。
void test_simulation_driver()
{
    int argc = 0;
    char** argv = nullptr;
    wcns::MpiRuntime mpi(argc, argv);
    {
        FakeSolver solver;
        TimeEventObserver observer;
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Unsteady;
        config.end_time = 1.0;
        config.max_steps = 4;
        wcns::Real clock = 0.0;
        wcns::SimulationDriver driver(
            mpi,
            solver,
            config,
            observer,
            {},
            [&clock] { clock += 0.01; return clock; });
        const auto final = driver.run();
        WCNS_REQUIRE(final.stop_reason == wcns::StopReason::PhysicalTimeReached);
        WCNS_REQUIRE(final.step == 4);
        WCNS_REQUIRE_NEAR(final.time, 1.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(solver.time_steps[0], 0.3, 1.0e-14);
        WCNS_REQUIRE_NEAR(solver.time_steps[1], 0.2, 1.0e-14);
        WCNS_REQUIRE_NEAR(solver.time_steps[2], 0.3, 1.0e-14);
        WCNS_REQUIRE_NEAR(solver.time_steps[3], 0.2, 1.0e-14);
        WCNS_REQUIRE(observer.final_count == 1);
    }
    {
        FakeSolver solver;
        TimeEventObserver observer;
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Steady;
        config.max_steps = 10;
        config.steady.min_steps = 2;
        config.steady.consecutive_checks = 1;
        config.steady.l2_absolute = 1.0e-20;
        config.steady.l2_relative = 1.0e-3;
        config.steady.linf_absolute = 1.0e-20;
        config.steady.linf_relative = 1.0e-3;
        wcns::SimulationDriver driver(mpi, solver, config, observer);
        const auto final = driver.run();
        WCNS_REQUIRE(final.stop_reason == wcns::StopReason::SteadyConverged);
        WCNS_REQUIRE(final.step == 2);
        WCNS_REQUIRE(final.steady.consecutive_passes == 1);
    }
}
