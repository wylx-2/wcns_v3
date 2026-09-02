#include "test_support.hpp"

#include <cgnslib.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string join(const std::string& directory, const std::string& name)
{
    return directory + "/" + name;
}

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(
            std::string(operation) + ": " + cg_get_error());
    }
}

std::string read_text(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open output text file: " + path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

void validate_cgns(const std::string& path)
{
    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_READ, &file), "cg_open field output");
    int bases = 0;
    int zones = 0;
    int solutions = 0;
    int fields = 0;
    check_cgns(cg_nbases(file, &bases), "cg_nbases field output");
    check_cgns(cg_nzones(file, 1, &zones), "cg_nzones field output");
    check_cgns(cg_nsols(file, 1, 1, &solutions), "cg_nsols field output");
    check_cgns(cg_nfields(file, 1, 1, 1, &fields), "cg_nfields field output");
    WCNS_REQUIRE(bases == 1);
    WCNS_REQUIRE(zones == 1);
    WCNS_REQUIRE(solutions == 1);
    WCNS_REQUIRE(fields == 12);

    char zone_name[33] = {};
    cgsize_t size[9] = {};
    check_cgns(cg_zone_read(file, 1, 1, zone_name, size), "cg_zone_read field output");
    WCNS_REQUIRE(size[0] == 33);
    WCNS_REQUIRE(size[1] == 17);
    WCNS_REQUIRE(size[2] == 32);
    WCNS_REQUIRE(size[3] == 16);
    std::vector<double> values(static_cast<std::size_t>(size[2] * size[3]));
    cgsize_t lower[2] = {1, 1};
    cgsize_t upper[2] = {size[2], size[3]};
    check_cgns(
        cg_field_read(
            file, 1, 1, 1, "Density", RealDouble,
            lower, upper, values.data()),
        "cg_field_read Density");
    for (const double value : values) WCNS_REQUIRE_NEAR(value, 1.0, 1.0e-12);
    check_cgns(
        cg_field_read(
            file, 1, 1, 1, "VelocityX", RealDouble,
            lower, upper, values.data()),
        "cg_field_read VelocityX");
    for (const double value : values) WCNS_REQUIRE_NEAR(value, 0.2, 1.0e-12);
    check_cgns(cg_close(file), "cg_close field output");
}

} // namespace

// 独立重读 CGNS/Tecplot/历史/统计文件，防止只验证“写文件未报错”。
int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr
            << "usage: wcns_field_output_tests <output-directory> "
               "<rank-count> <history-format>\n";
        return EXIT_FAILURE;
    }
    try {
        const std::string directory = argv[1];
        const std::string ranks = argv[2];
        const std::string stem
            = "output-freestream.field.step00000000.time0p000000000eP00";
        validate_cgns(join(directory, stem + ".cgns"));
        const auto tecplot = read_text(join(directory, stem + ".dat"));
        WCNS_REQUIRE(tecplot.find("DATAPACKING=POINT") != std::string::npos);
        WCNS_REQUIRE(tecplot.find("\"rho\"") != std::string::npos);
        const std::string history_format = argv[3];
        const auto history = read_text(join(
            directory, "output-freestream.history.r" + ranks
                + (history_format == "tecplot" ? ".dat" : ".txt")));
        if (history_format == "tecplot") {
            WCNS_REQUIRE(history.find("stop_reason_code") != std::string::npos);
            WCNS_REQUIRE(history.find("STOP_REASON_CODES") != std::string::npos);
        } else {
            WCNS_REQUIRE(history.find("steady_converged") != std::string::npos);
        }
        const auto statistics = read_text(join(
            directory, "output-freestream.statistics.r" + ranks + ".dat"));
        WCNS_REQUIRE(statistics.find("total_mass") != std::string::npos);
        std::cout << "production output files independently re-read\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
