#include "test_support.hpp"

#include <wcns/mesh/structured_mesh.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

namespace {

wcns::ConnectivityPatch receiver_connection()
{
    return {
        "left-to-right",
        0,
        1,
        1,
        {wcns::Axis::I, wcns::Side::Upper},
        {wcns::Axis::J, wcns::Side::Lower},
        {{4, 0, 0}, {4, 3, 0}},
        {{3, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{{2, -1, 3}}},
        3,
    };
}

wcns::ConnectivityPatch donor_connection()
{
    return {
        "right-to-left",
        1,
        0,
        0,
        {wcns::Axis::J, wcns::Side::Lower},
        {wcns::Axis::I, wcns::Side::Upper},
        {{3, 0, 0}, {0, 0, 0}},
        {{4, 0, 0}, {4, 3, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{{-2, 1, 3}}},
        3,
    };
}

} // namespace

void test_topology()
{
    using namespace wcns;

    const IndexTransform transform {{{2, -1, 3}}};
    WCNS_REQUIRE(transform.valid(2));
    WCNS_REQUIRE(transform.map({4, 3, 0}, {4, 0, 0}, {3, 0, 0}, 2)
        == (Index3 {0, 0, 0}));
    WCNS_REQUIRE(transform.inverse(2) == (IndexTransform {{{-2, 1, 3}}}));
    WCNS_REQUIRE((!IndexTransform {{{2, 1, 3}}}.valid(2)));
    WCNS_REQUIRE((!IndexTransform {{{1, 2, -3}}}.valid(2)));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (IndexTransform {{{2, 1, 3}}}.map({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 2)));

    StructuredBlock left(0, "left", 0, 2, 2, {5, 4, 1}, 3);
    StructuredBlock right(1, "right", 1, 2, 2, {4, 5, 1}, 3);
    left.connectivities.push_back(receiver_connection());
    right.connectivities.push_back(donor_connection());

    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    StructuredMesh mesh(std::move(blocks));
    WCNS_REQUIRE(mesh.block_count() == 2);
    WCNS_REQUIRE(mesh.contains(0));
    WCNS_REQUIRE(mesh.block(1).name() == "right");
    mesh.validate_connectivities();

    StructuredBlock incomplete_left(0, "left", 0, 2, 2, {5, 4, 1}, 3);
    StructuredBlock incomplete_right(1, "right", 1, 2, 2, {4, 5, 1}, 3);
    incomplete_left.connectivities.push_back(receiver_connection());
    std::vector<StructuredBlock> incomplete_blocks;
    incomplete_blocks.push_back(std::move(incomplete_left));
    incomplete_blocks.push_back(std::move(incomplete_right));
    StructuredMesh incomplete(std::move(incomplete_blocks));
    WCNS_REQUIRE_THROWS(TopologyError, incomplete.validate_connectivities());
    WCNS_REQUIRE_THROWS(std::out_of_range, incomplete.block(7));
}

