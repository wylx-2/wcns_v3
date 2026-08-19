#include "test_support.hpp"

#include <wcns/core/field.hpp>

#include <stdexcept>

void test_field()
{
    using wcns::Extent3;
    using wcns::Field;

    Field<double> state({2, 1, 1}, 5, 1, 0.0);
    WCNS_REQUIRE(state.interior_extent() == (Extent3 {2, 1, 1}));
    WCNS_REQUIRE(state.storage_extent() == (Extent3 {4, 3, 3}));
    WCNS_REQUIRE(state.components() == 5);
    WCNS_REQUIRE(state.ghost_width() == 1);
    WCNS_REQUIRE(state.cell_count() == 36);
    WCNS_REQUIRE(state.size() == 180);

    for (int component = 0; component < state.components(); ++component) {
        state(0, 0, 0, component) = 10.0 + component;
    }
    for (int component = 0; component < state.components(); ++component) {
        WCNS_REQUIRE(state(0, 0, 0, component) == 10.0 + component);
    }

    WCNS_REQUIRE(
        state.linear_index(0, 0, 0, 0) + 1 == state.linear_index(0, 0, 0, 1));
    WCNS_REQUIRE(
        state.linear_index(0, 0, 0, 4) + 1 == state.linear_index(1, 0, 0, 0));

    state(-1, -1, -1, 0) = -10.0;
    state(2, 1, 1, 4) = 20.0;
    WCNS_REQUIRE(state(-1, -1, -1, 0) == -10.0);
    WCNS_REQUIRE(state(2, 1, 1, 4) == 20.0);

    const auto& read_only = state;
    WCNS_REQUIRE(read_only(0, 0, 0, 3) == 13.0);
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(-2, 0, 0, 0));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(3, 0, 0, 0));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(0, 0, 0, -1));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(0, 0, 0, 5));

    state.fill(3.5);
    for (std::size_t index = 0; index < state.size(); ++index) {
        WCNS_REQUIRE(state.data()[index] == 3.5);
    }

    WCNS_REQUIRE_THROWS(std::invalid_argument, (Field<double>({1, 1, 1}, 0, 0)));
    WCNS_REQUIRE_THROWS(std::invalid_argument, (Field<double>({1, 1, 1}, 5, -1)));
    WCNS_REQUIRE_THROWS(std::invalid_argument, (Field<double>({-1, 1, 1}, 5, 0)));
}

