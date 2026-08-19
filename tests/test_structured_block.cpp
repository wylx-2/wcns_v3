#include "test_support.hpp"

#include <wcns/mesh/structured_block.hpp>

#include <stdexcept>

void test_structured_block()
{
    using namespace wcns;

    StructuredBlock block(7, "zone-7", 1, 2, 2, {5, 4, 1}, 3);
    WCNS_REQUIRE(block.id() == 7);
    WCNS_REQUIRE(block.name() == "zone-7");
    WCNS_REQUIRE(block.owner_rank() == 1);
    WCNS_REQUIRE(block.cell_dimension() == 2);
    WCNS_REQUIRE(block.physical_dimension() == 2);
    WCNS_REQUIRE(block.vertex_extent() == (Extent3 {5, 4, 1}));
    WCNS_REQUIRE(block.cell_extent() == (Extent3 {4, 3, 1}));
    WCNS_REQUIRE(block.ghost_width() == 3);

    WCNS_REQUIRE(block.coordinates.x.interior_extent() == block.vertex_extent());
    WCNS_REQUIRE(block.coordinates.x.ghost_width() == 0);
    WCNS_REQUIRE(block.cell_metrics.volume.interior_extent() == block.cell_extent());
    WCNS_REQUIRE(block.cell_metrics.volume.ghost_width() == 3);

    WCNS_REQUIRE(
        block.face_metrics.i_faces.area.interior_extent() == (Extent3 {5, 3, 1}));
    WCNS_REQUIRE(
        block.face_metrics.j_faces.area.interior_extent() == (Extent3 {4, 4, 1}));
    WCNS_REQUIRE(
        block.face_metrics.k_faces.area.interior_extent() == (Extent3 {4, 3, 2}));

    WCNS_REQUIRE(block.flow.conservative.components() == euler_components);
    WCNS_REQUIRE(block.flow.conservative.ghost_width() == 3);
    WCNS_REQUIRE(block.flow.primitive.ghost_width() == 3);
    WCNS_REQUIRE(block.flow.residual.ghost_width() == 0);

    block.set_owner_rank(3);
    WCNS_REQUIRE(block.owner_rank() == 3);
    block.set_owner_rank(invalid_rank_id);
    WCNS_REQUIRE(block.owner_rank() == invalid_rank_id);

    BoundaryPatch wall;
    wall.name = "wall";
    wall.type = BoundaryType::SlipWall;
    wall.face = {Axis::J, Side::Lower};
    block.boundaries.push_back(wall);
    WCNS_REQUIRE(block.boundaries.size() == 1);
    WCNS_REQUIRE(block.boundaries.front().face == (FaceLocation {Axis::J, Side::Lower}));

    ConnectivityPatch connection;
    connection.name = "to-zone-8";
    connection.receiver_block = block.id();
    connection.donor_block = 8;
    connection.transform.receiver_to_donor = {{2, -1, 3}};
    connection.ghost_width = block.ghost_width();
    block.connectivities.push_back(connection);
    WCNS_REQUIRE(block.connectivities.front().transform.valid());
    WCNS_REQUIRE((!IndexTransform {{{1, 1, 3}}}.valid()));

    StructuredBlock volume(8, "volume", 0, 3, 3, {4, 3, 2}, 2);
    WCNS_REQUIRE(volume.cell_extent() == (Extent3 {3, 2, 1}));

    WCNS_REQUIRE_THROWS(
        std::invalid_argument, (StructuredBlock(-1, "bad", 0, 2, 2, {2, 2, 1}, 3)));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, (StructuredBlock(0, "", 0, 2, 2, {2, 2, 1}, 3)));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, (StructuredBlock(0, "bad", 0, 2, 2, {2, 2, 2}, 3)));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, (StructuredBlock(0, "bad", 0, 3, 3, {2, 2, 1}, 3)));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument, (StructuredBlock(0, "bad", -2, 2, 2, {2, 2, 1}, 3)));
}
