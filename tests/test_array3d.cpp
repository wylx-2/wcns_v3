#include "test_support.hpp"

#include <wcns/core/array3d.hpp>

#include <limits>
#include <stdexcept>

void test_array3d()
{
    using wcns::Array3D;
    using wcns::Extent3;

    Array3D<double> values({2, 3, 1}, 2, -1.0);
    WCNS_REQUIRE(values.interior_extent() == (Extent3 {2, 3, 1}));
    WCNS_REQUIRE(values.storage_extent() == (Extent3 {6, 7, 5}));
    WCNS_REQUIRE(values.ghost_width() == 2);
    WCNS_REQUIRE(values.size() == 210);

    values(0, 0, 0) = 10.0;
    values(1, 2, 0) = 20.0;
    values(-2, -2, -2) = 30.0;
    values(3, 4, 2) = 40.0;

    WCNS_REQUIRE(values(0, 0, 0) == 10.0);
    WCNS_REQUIRE(values(1, 2, 0) == 20.0);
    WCNS_REQUIRE(values(-2, -2, -2) == 30.0);
    WCNS_REQUIRE(values(3, 4, 2) == 40.0);

    WCNS_REQUIRE(values.linear_index(-2, -2, -2) == 0);
    WCNS_REQUIRE(values.linear_index(3, 4, 2) == values.size() - 1);
    WCNS_REQUIRE(values.linear_index(0, 0, 0) + 1 == values.linear_index(1, 0, 0));

    const auto& read_only = values;
    WCNS_REQUIRE(read_only(-2, -2, -2) == 30.0);
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(-3, 0, 0));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(4, 0, 0));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(0, 5, 0));
    WCNS_REQUIRE_THROWS(std::out_of_range, read_only(0, 0, 3));

    values.fill(7.0);
    for (std::size_t index = 0; index < values.size(); ++index) {
        WCNS_REQUIRE(values.data()[index] == 7.0);
    }

    WCNS_REQUIRE_THROWS(std::invalid_argument, (Array3D<double>({-1, 2, 3}, 0)));
    WCNS_REQUIRE_THROWS(std::invalid_argument, (Array3D<double>({1, 2, 3}, -1)));
    WCNS_REQUIRE_THROWS(
        std::overflow_error,
        (Array3D<double>({1, 1, 1}, std::numeric_limits<int>::max())));
}

