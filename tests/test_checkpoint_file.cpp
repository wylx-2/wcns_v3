#include "test_support.hpp"

#include <cgnslib.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(
            std::string(operation) + ": " + cg_get_error());
    }
}

std::map<std::string, std::string> descriptors(int file)
{
    check_cgns(cg_goto(file, 1, "end"), "cg_goto checkpoint test");
    int count = 0;
    check_cgns(cg_ndescriptors(&count), "cg_ndescriptors checkpoint test");
    std::map<std::string, std::string> result;
    for (int index = 1; index <= count; ++index) {
        char name[33] = {};
        char* value = nullptr;
        check_cgns(
            cg_descriptor_read(index, name, &value),
            "cg_descriptor_read checkpoint test");
        result.emplace(name, value == nullptr ? "" : value);
        if (value != nullptr) cg_free(value);
    }
    return result;
}

} // namespace

// 独立重读检查点的网格、五个守恒场及严格重启元数据。
int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: wcns_checkpoint_file_tests <checkpoint.cgns>\n";
        return EXIT_FAILURE;
    }
    int file = 0;
    try {
        check_cgns(cg_open(argv[1], CG_MODE_READ, &file), "cg_open checkpoint test");
        const auto metadata = descriptors(file);
        WCNS_REQUIRE(metadata.at("WCNS_Version") == "1");
        WCNS_REQUIRE(metadata.at("WCNS_Step") == "1");
        WCNS_REQUIRE(!metadata.at("WCNS_MeshSignature").empty());
        WCNS_REQUIRE(!metadata.at("WCNS_RestartSignature").empty());
        int coordinates = 0;
        int fields = 0;
        check_cgns(cg_ncoords(file, 1, 1, &coordinates), "cg_ncoords checkpoint test");
        check_cgns(cg_nfields(file, 1, 1, 1, &fields), "cg_nfields checkpoint test");
        WCNS_REQUIRE(coordinates == 2);
        WCNS_REQUIRE(fields == 5);
        char zone_name[33] = {};
        cgsize_t size[9] = {};
        check_cgns(cg_zone_read(file, 1, 1, zone_name, size), "cg_zone_read checkpoint test");
        std::vector<double> density(static_cast<std::size_t>(size[2] * size[3]));
        cgsize_t lower[2] = {1, 1};
        cgsize_t upper[2] = {size[2], size[3]};
        check_cgns(
            cg_field_read(
                file, 1, 1, 1, "Density", RealDouble,
                lower, upper, density.data()),
            "cg_field_read checkpoint Density");
        for (const double value : density) WCNS_REQUIRE_NEAR(value, 1.0, 1.0e-12);
        check_cgns(cg_close(file), "cg_close checkpoint test");
        file = 0;
        std::cout << "checkpoint independently re-read\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        if (file > 0) cg_close(file);
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
