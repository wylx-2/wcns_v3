#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/metrics.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void test_2d(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.bases.size() == 1);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    WCNS_REQUIRE(metadata.bases.front().name == "Base2D");
    const auto& zone = metadata.zones.front();
    WCNS_REQUIRE(zone.block_id == 0);
    WCNS_REQUIRE(zone.name == "Zone2D");
    WCNS_REQUIRE(zone.cell_dimension == 2);
    WCNS_REQUIRE(zone.physical_dimension == 2);
    WCNS_REQUIRE(zone.vertex_extent == (wcns::Extent3 {5, 4, 1}));
    WCNS_REQUIRE(zone.cell_extent == (wcns::Extent3 {4, 3, 1}));

    auto block = reader.read_block(path, zone, 0, 3);
    WCNS_REQUIRE(block.id() == 0);
    WCNS_REQUIRE(block.owner_rank() == 0);
    WCNS_REQUIRE(block.ghost_width() == 3);
    WCNS_REQUIRE(block.coordinates.x(0, 0, 0) == 0.0);
    WCNS_REQUIRE(block.coordinates.x(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.y(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(4, 3, 0) == 0.0);
    WCNS_REQUIRE(block.boundaries.size() == 4);

    const auto& imin = block.boundaries[0];
    WCNS_REQUIRE(imin.name == "imin");
    WCNS_REQUIRE(imin.type == wcns::BoundaryType::Inflow);
    WCNS_REQUIRE(imin.face == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Lower}));
    WCNS_REQUIRE(imin.vertex_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 3, 0}}));
    WCNS_REQUIRE(imin.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 2, 0}}));
    WCNS_REQUIRE(imin.boundary_face_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 2, 0}}));

    const auto& imax = block.boundaries[1];
    WCNS_REQUIRE(imax.type == wcns::BoundaryType::Outflow);
    WCNS_REQUIRE(imax.face == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(imax.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{3, 0, 0}, {3, 2, 0}}));
    WCNS_REQUIRE(imax.boundary_face_range.untyped() == (wcns::IndexRange3 {{4, 0, 0}, {4, 2, 0}}));

    const auto& jmax = block.boundaries[3];
    WCNS_REQUIRE(jmax.type == wcns::BoundaryType::SlipWall);
    WCNS_REQUIRE(jmax.face == (wcns::FaceLocation {wcns::Axis::J, wcns::Side::Upper}));
    WCNS_REQUIRE(jmax.vertex_range.untyped() == (wcns::IndexRange3 {{4, 3, 0}, {0, 3, 0}}));
    WCNS_REQUIRE(jmax.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{3, 2, 0}, {0, 2, 0}}));
    WCNS_REQUIRE(jmax.boundary_face_range.untyped() == (wcns::IndexRange3 {{3, 3, 0}, {0, 3, 0}}));

    wcns::compute_metrics(block);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_x(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_y(0, 0, 0), 1.0 / 6.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.volume(0, 0, 0), 1.0 / 12.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.jacobian(3, 2, 0), 1.0 / 12.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.area(0, 0, 0), 1.0 / 3.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_x(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_y(0, 0, 0), 0.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.area(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.normal_y(0, 0, 0), 1.0, 1.0e-14);

    double total_area = 0.0;
    for (int j = 0; j < block.cell_extent().nj; ++j) {
        for (int i = 0; i < block.cell_extent().ni; ++i) {
            total_area += block.cell_metrics.volume(i, j, 0);
        }
    }
    WCNS_REQUIRE_NEAR(total_area, 1.0, 1.0e-13);
}

void test_3d(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.bases.size() == 1);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    const auto& zone = metadata.zones.front();
    WCNS_REQUIRE(zone.name == "Zone3D");
    WCNS_REQUIRE(zone.cell_dimension == 3);
    WCNS_REQUIRE(zone.physical_dimension == 3);
    WCNS_REQUIRE(zone.vertex_extent == (wcns::Extent3 {4, 3, 3}));
    WCNS_REQUIRE(zone.cell_extent == (wcns::Extent3 {3, 2, 2}));

    auto block = reader.read_block(path, zone, 2, 2);
    WCNS_REQUIRE(block.owner_rank() == 2);
    WCNS_REQUIRE(block.coordinates.x(3, 2, 2) == 3.0);
    WCNS_REQUIRE(block.coordinates.y(3, 2, 2) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(3, 2, 2) == 0.5);
    WCNS_REQUIRE(block.boundaries.size() == 6);
    WCNS_REQUIRE(block.boundaries[0].type == wcns::BoundaryType::Farfield);
    WCNS_REQUIRE(
        block.boundaries[3].adjacent_cell_range.untyped()
        == (wcns::IndexRange3 {{2, 1, 0}, {0, 1, 1}}));
    WCNS_REQUIRE(
        block.boundaries[5].face
        == (wcns::FaceLocation {wcns::Axis::K, wcns::Side::Upper}));

    wcns::compute_metrics(block);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_x(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_y(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_z(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.volume(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.area(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_x(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.area(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.normal_y(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.k_faces.area(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.k_faces.normal_z(0, 0, 0), 1.0, 1.0e-14);

    double total_volume = 0.0;
    for (int k = 0; k < block.cell_extent().nk; ++k) {
        for (int j = 0; j < block.cell_extent().nj; ++j) {
            for (int i = 0; i < block.cell_extent().ni; ++i) {
                total_volume += block.cell_metrics.volume(i, j, k);
            }
        }
    }
    WCNS_REQUIRE_NEAR(total_volume, 1.5, 1.0e-13);
}

void test_out_of_extent_boundary(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    try {
        static_cast<void>(reader.read_block(path, metadata.zones.front(), 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find(
                "boundary vertex PointRange end J index 4 is outside [0, 3]")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("out-of-extent CGNS boundary was accepted");
}

void test_multiblock_2d(const char* path)
{
    wcns::CgnsReader reader;
    auto mesh = reader.read_mesh(path, 2, 3);
    WCNS_REQUIRE(mesh.block_count() == 2);
    const auto& left = mesh.block(0);
    const auto& right = mesh.block(1);
    WCNS_REQUIRE(left.name() == "Left2D");
    WCNS_REQUIRE(right.name() == "Right2D");
    WCNS_REQUIRE(left.connectivities.size() == 1);
    WCNS_REQUIRE(right.connectivities.size() == 1);

    const auto& connection = left.connectivities.front();
    WCNS_REQUIRE(connection.receiver_block == 0);
    WCNS_REQUIRE(connection.donor_block == 1);
    WCNS_REQUIRE(connection.donor_rank == 2);
    WCNS_REQUIRE(
        connection.receiver_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(
        connection.donor_face
        == (wcns::FaceLocation {wcns::Axis::J, wcns::Side::Lower}));
    WCNS_REQUIRE(
        connection.receiver_vertex_range.untyped()
        == (wcns::IndexRange3 {{4, 0, 0}, {4, 3, 0}}));
    WCNS_REQUIRE(
        connection.donor_vertex_range.untyped()
        == (wcns::IndexRange3 {{3, 0, 0}, {0, 0, 0}}));
    WCNS_REQUIRE(connection.transform == (wcns::IndexTransform {{{2, -1, 3}}}));
    WCNS_REQUIRE(
        connection.shared_face_range.untyped()
        == (wcns::IndexRange3 {{4, 0, 0}, {4, 2, 0}}));
    WCNS_REQUIRE(connection.ghost_width == 3);
}

void test_multiblock_3d(const char* path)
{
    wcns::CgnsReader reader;
    auto mesh = reader.read_mesh(path, 0, 2);
    WCNS_REQUIRE(mesh.block_count() == 2);
    const auto& connection = mesh.block(0).connectivities.front();
    WCNS_REQUIRE(connection.receiver_block == 0);
    WCNS_REQUIRE(connection.donor_block == 1);
    WCNS_REQUIRE(
        connection.receiver_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(
        connection.donor_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Lower}));
    WCNS_REQUIRE(connection.transform == (wcns::IndexTransform {{{1, -2, 3}}}));
    WCNS_REQUIRE(
        connection.donor_vertex_range.untyped()
        == (wcns::IndexRange3 {{0, 2, 0}, {0, 0, 2}}));
    WCNS_REQUIRE(
        connection.donor_adjacent_cell_range.untyped()
        == (wcns::IndexRange3 {{0, 1, 0}, {0, 0, 1}}));
}

void test_one_sided_connectivity(const char* path)
{
    wcns::CgnsReader reader;
    try {
        static_cast<void>(reader.read_mesh(path, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find("has no reciprocal donor record")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("one-sided CGNS connectivity was accepted");
}

void test_unknown_donor(const char* path)
{
    wcns::CgnsReader reader;
    try {
        static_cast<void>(reader.read_mesh(path, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find("unknown donor zone: Missing2D")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("unknown CGNS donor zone was accepted");
}

void test_out_of_extent_connectivity(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    auto reduced_zone = metadata.zones.front();
    reduced_zone.vertex_extent.nj = 3;
    reduced_zone.cell_extent.nj = 2;
    try {
        static_cast<void>(reader.read_block(path, reduced_zone, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find(
                "connectivity receiver vertex range end J index 3 is outside [0, 2]")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("out-of-extent CGNS connectivity was accepted");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 8) {
        std::cerr
            << "usage: wcns_cgns_reader_tests <2d.cgns> <3d.cgns> "
               "<invalid-2d.cgns> <multi-2d.cgns> <multi-3d.cgns> "
               "<one-sided-2d.cgns> <unknown-donor-2d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        test_2d(argv[1]);
        test_3d(argv[2]);
        test_out_of_extent_boundary(argv[3]);
        test_multiblock_2d(argv[4]);
        test_multiblock_3d(argv[5]);
        test_one_sided_connectivity(argv[6]);
        test_unknown_donor(argv[7]);
        test_out_of_extent_connectivity(argv[4]);
        std::cout << "CGNS reader tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
