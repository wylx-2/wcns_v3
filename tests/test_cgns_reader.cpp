#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>

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

    const auto block = reader.read_block(path, zone, 0, 3);
    WCNS_REQUIRE(block.id() == 0);
    WCNS_REQUIRE(block.owner_rank() == 0);
    WCNS_REQUIRE(block.ghost_width() == 3);
    WCNS_REQUIRE(block.coordinates.x(0, 0, 0) == 0.0);
    WCNS_REQUIRE(block.coordinates.x(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.y(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(4, 3, 0) == 0.0);
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

    const auto block = reader.read_block(path, zone, 2, 2);
    WCNS_REQUIRE(block.owner_rank() == 2);
    WCNS_REQUIRE(block.coordinates.x(3, 2, 2) == 3.0);
    WCNS_REQUIRE(block.coordinates.y(3, 2, 2) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(3, 2, 2) == 0.5);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: wcns_cgns_reader_tests <2d.cgns> <3d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        test_2d(argv[1]);
        test_3d(argv[2]);
        std::cout << "CGNS reader tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

