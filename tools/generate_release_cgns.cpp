#include <cgnslib.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(std::string(operation) + ": " + cg_get_error());
    }
}

int parse_positive(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const long value = std::stol(text, &consumed);
    if (consumed != std::string(text).size() || value <= 0 || value > 1000000) {
        throw std::invalid_argument(std::string(name) + " must be in [1,1000000]");
    }
    return static_cast<int>(value);
}

double parse_warp(const char* text)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)
        || std::abs(value) >= 0.2) {
        throw std::invalid_argument("warp must be finite with absolute value < 0.2");
    }
    return value;
}

bool parse_bool(const char* text)
{
    const std::string value(text);
    if (value == "true") return true;
    if (value == "false") return false;
    throw std::invalid_argument("periodic_x must be true or false");
}

void write_boundary(
    int file,
    int base,
    int zone,
    const std::string& name,
    const std::vector<cgsize_t>& range)
{
    int boundary = 0;
    check_cgns(
        cg_boco_write(
            file, base, zone, name.c_str(), BCFarfield,
            PointRange, 2, range.data(), &boundary),
        "cg_boco_write release grid");
    check_cgns(
        cg_boco_gridlocation_write(file, base, zone, boundary, Vertex),
        "cg_boco_gridlocation_write release grid");
}

int write_connection(
    int file,
    int base,
    int zone,
    const std::string& name,
    const std::string& donor,
    const std::vector<cgsize_t>& range,
    const std::vector<cgsize_t>& donor_range,
    int dimension)
{
    std::array<int, 3> transform {{1, 2, 3}};
    int connection = 0;
    check_cgns(
        cg_1to1_write(
            file, base, zone, name.c_str(), donor.c_str(),
            range.data(), donor_range.data(), transform.data(), &connection),
        "cg_1to1_write release grid");
    static_cast<void>(dimension);
    return connection;
}

std::vector<cgsize_t> face_range(
    int dimension,
    int i,
    int nj,
    int nk)
{
    if (dimension == 2) return {i, 1, i, nj};
    return {i, 1, 1, i, nj, nk};
}

std::string zone_name(int index)
{
    return "Zone" + std::to_string(index + 1);
}

void generate(
    const std::string& path,
    int dimension,
    int cells_i,
    int cells_j,
    int cells_k,
    int zones_i,
    double warp,
    bool periodic_x)
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("dimension must be 2 or 3");
    }
    if (dimension == 2 && cells_k != 1) {
        throw std::invalid_argument("two-dimensional grids require cells_k=1");
    }
    if (cells_i % zones_i != 0) {
        throw std::invalid_argument("cells_i must be divisible by zones_i");
    }
    if (periodic_x && zones_i < 2) {
        throw std::invalid_argument(
            "periodic_x requires at least two zones; self-connectivity is not emitted");
    }
    const int local_cells_i = cells_i / zones_i;
    if (local_cells_i < 1) {
        throw std::invalid_argument("each zone must contain an i cell");
    }

    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &file), "cg_open release grid");
    try {
        int base = 0;
        check_cgns(
            cg_base_write(file, "WCNSReleaseGrid", dimension, dimension, &base),
            "cg_base_write release grid");
        const int ni = local_cells_i + 1;
        const int nj = cells_j + 1;
        const int nk = dimension == 3 ? cells_k + 1 : 1;
        for (int zone_index = 0; zone_index < zones_i; ++zone_index) {
            std::array<cgsize_t, 9> size {{}};
            size[0] = ni;
            size[1] = nj;
            size[2] = dimension == 3 ? nk : local_cells_i;
            size[dimension] = local_cells_i;
            size[dimension + 1] = cells_j;
            if (dimension == 3) size[dimension + 2] = cells_k;
            int zone = 0;
            check_cgns(
                cg_zone_write(
                    file, base, zone_name(zone_index).c_str(),
                    size.data(), Structured, &zone),
                "cg_zone_write release grid");

            const std::size_t count = static_cast<std::size_t>(ni)
                * static_cast<std::size_t>(nj) * static_cast<std::size_t>(nk);
            std::vector<double> x(count);
            std::vector<double> y(count);
            std::vector<double> z(count);
            constexpr double pi = 3.141592653589793238462643383279502884;
            for (int k = 0; k < nk; ++k) {
                const double zeta = dimension == 3
                    ? static_cast<double>(k) / static_cast<double>(cells_k) : 0.0;
                for (int j = 0; j < nj; ++j) {
                    const double eta
                        = static_cast<double>(j) / static_cast<double>(cells_j);
                    for (int i = 0; i < ni; ++i) {
                        const int global_i = zone_index * local_cells_i + i;
                        const double xi
                            = static_cast<double>(global_i) / static_cast<double>(cells_i);
                        const double envelope = dimension == 3
                            ? std::sin(pi * zeta) : 1.0;
                        const auto index = static_cast<std::size_t>((k * nj + j) * ni + i);
                        x[index] = xi + warp * std::sin(pi * xi)
                            * std::sin(pi * eta) * envelope;
                        y[index] = eta;
                        z[index] = zeta;
                    }
                }
            }
            int coordinate = 0;
            check_cgns(
                cg_coord_write(
                    file, base, zone, RealDouble,
                    "CoordinateX", x.data(), &coordinate),
                "cg_coord_write release X");
            check_cgns(
                cg_coord_write(
                    file, base, zone, RealDouble,
                    "CoordinateY", y.data(), &coordinate),
                "cg_coord_write release Y");
            if (dimension == 3) {
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateZ", z.data(), &coordinate),
                    "cg_coord_write release Z");
            }

            if (dimension == 2) {
                write_boundary(file, base, zone, "bottom", {1, 1, ni, 1});
                write_boundary(file, base, zone, "top", {1, nj, ni, nj});
            } else {
                write_boundary(file, base, zone, "bottom", {1, 1, 1, ni, 1, nk});
                write_boundary(file, base, zone, "top", {1, nj, 1, ni, nj, nk});
                write_boundary(file, base, zone, "front", {1, 1, 1, ni, nj, 1});
                write_boundary(file, base, zone, "back", {1, 1, nk, ni, nj, nk});
            }
            if (!periodic_x && zone_index == 0) {
                write_boundary(
                    file, base, zone, "left",
                    face_range(dimension, 1, nj, nk));
            }
            if (!periodic_x && zone_index + 1 == zones_i) {
                write_boundary(
                    file, base, zone, "right",
                    face_range(dimension, ni, nj, nk));
            }
        }

        for (int zone_index = 0; zone_index + 1 < zones_i; ++zone_index) {
            const int left = zone_index + 1;
            const int right = zone_index + 2;
            write_connection(
                file, base, left,
                "to-" + zone_name(zone_index + 1), zone_name(zone_index + 1),
                face_range(dimension, ni, nj, nk),
                face_range(dimension, 1, nj, nk), dimension);
            write_connection(
                file, base, right,
                "to-" + zone_name(zone_index), zone_name(zone_index),
                face_range(dimension, 1, nj, nk),
                face_range(dimension, ni, nj, nk), dimension);
        }
        if (periodic_x) {
            const int forward = write_connection(
                file, base, 1, "periodic-left", zone_name(zones_i - 1),
                face_range(dimension, 1, nj, nk),
                face_range(dimension, ni, nj, nk), dimension);
            const int reverse = write_connection(
                file, base, zones_i, "periodic-right", zone_name(0),
                face_range(dimension, ni, nj, nk),
                face_range(dimension, 1, nj, nk), dimension);
            std::array<float, 3> center {{0.0F, 0.0F, 0.0F}};
            std::array<float, 3> angle {{0.0F, 0.0F, 0.0F}};
            std::array<float, 3> translation {{1.0F, 0.0F, 0.0F}};
            check_cgns(
                cg_1to1_periodic_write(
                    file, base, 1, forward,
                    center.data(), angle.data(), translation.data()),
                "cg_1to1_periodic_write forward");
            translation[0] = -1.0F;
            check_cgns(
                cg_1to1_periodic_write(
                    file, base, zones_i, reverse,
                    center.data(), angle.data(), translation.data()),
                "cg_1to1_periodic_write reverse");
        }
        check_cgns(cg_close(file), "cg_close release grid");
        file = 0;
    } catch (...) {
        if (file != 0) cg_close(file);
        throw;
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 9) {
        std::cerr
            << "usage: wcns_generate_release_cgns <output.cgns> <dimension> "
               "<cells_i> <cells_j> <cells_k> <zones_i> <warp> <periodic_x>\n";
        return EXIT_FAILURE;
    }
    try {
        check_cgns(cg_set_file_type(CG_FILE_ADF), "cg_set_file_type release grid");
        generate(
            argv[1], parse_positive(argv[2], "dimension"),
            parse_positive(argv[3], "cells_i"),
            parse_positive(argv[4], "cells_j"),
            parse_positive(argv[5], "cells_k"),
            parse_positive(argv[6], "zones_i"),
            parse_warp(argv[7]), parse_bool(argv[8]));
        std::cout << "generated release CGNS grid: " << argv[1] << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "release grid generation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
