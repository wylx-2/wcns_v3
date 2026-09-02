#include "test_support.hpp"

#include <wcns/runtime/output_manager.hpp>

// 验收步/时间/显式/初末事件取并集，并在同一 step/time 只消费一次。
void test_output_manager()
{
    wcns::OutputScheduleConfig config;
    config.every_steps = 2;
    config.every_time = 0.5;
    config.explicit_times = {0.75};
    config.write_initial = true;
    config.write_final = true;
    wcns::OutputSchedule schedule(config);

    wcns::SimulationState state;
    WCNS_REQUIRE(schedule.consume(state, true, false));
    WCNS_REQUIRE(!schedule.consume(state, true, false));
    WCNS_REQUIRE_NEAR(schedule.next_time(0.0), 0.5, 1.0e-15);

    state.step = 1;
    state.time = 0.5;
    WCNS_REQUIRE(schedule.consume(state, false, false));
    WCNS_REQUIRE(!schedule.consume(state, false, true));
    WCNS_REQUIRE_NEAR(schedule.next_time(0.5), 0.75, 1.0e-15);

    state.step = 2;
    state.time = 0.75;
    WCNS_REQUIRE(schedule.consume(state, false, false));
    WCNS_REQUIRE(!schedule.consume(state, false, true));

    state.step = 3;
    state.time = 0.8;
    WCNS_REQUIRE(schedule.consume(state, false, true));
    WCNS_REQUIRE(!schedule.consume(state, false, true));
}
