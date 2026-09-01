#include "test_support.hpp"

#include <wcns/core/index.hpp>

#include <stdexcept>
#include <limits>

// 验收三维索引、范围方向、计数和越界检查。
void test_index()
{
    using namespace wcns;

    Index3 index {1, 2, 3};
    WCNS_REQUIRE(index[0] == 1);
    WCNS_REQUIRE(index[1] == 2);
    WCNS_REQUIRE(index[2] == 3);
    WCNS_REQUIRE_THROWS(std::out_of_range, index[3]);

    const Extent3 extent {4, 3, 2};
    WCNS_REQUIRE(extent.valid());
    WCNS_REQUIRE(!extent.empty());
    WCNS_REQUIRE(extent.size() == 24);
    WCNS_REQUIRE((Extent3 {4, 3, 0}.empty()));
    WCNS_REQUIRE_THROWS(std::invalid_argument, (Extent3 {-1, 2, 3}.size()));
    WCNS_REQUIRE_THROWS(
        std::overflow_error,
        (Extent3 {
             std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max()}
             .size()));

    const IndexRange3 forward {{1, 2, 3}, {3, 5, 3}};
    WCNS_REQUIRE(forward.step() == (Index3 {1, 1, 0}));
    WCNS_REQUIRE(forward.counts() == (Extent3 {3, 4, 1}));
    WCNS_REQUIRE(forward.size() == 12);
    WCNS_REQUIRE(forward.at({2, 3, 0}) == (Index3 {3, 5, 3}));

    const IndexRange3 reversed {{4, 7, 1}, {2, 3, 1}};
    WCNS_REQUIRE(reversed.step() == (Index3 {-1, -1, 0}));
    WCNS_REQUIRE(reversed.counts() == (Extent3 {3, 5, 1}));
    WCNS_REQUIRE(reversed.at({1, 2, 0}) == (Index3 {3, 5, 1}));
    WCNS_REQUIRE_THROWS(std::out_of_range, reversed.at({3, 0, 0}));
}
