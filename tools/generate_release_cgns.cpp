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

double parse_positive_real(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)
        || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive and finite");
    }
    return value;
}

double parse_nonnegative_real(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)
        || value < 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and nonnegative");
    }
    return value;
}

double parse_finite_real(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
    return value;
}

double clustered_unit_coordinate(double logical, double center, double strength)
{
    if (!(center > 0.0 && center < 1.0) || !(strength > 0.0)) {
        throw std::invalid_argument(
            "cluster center must be inside the interval and strength positive");
    }
    const double denominator = std::sinh(strength);
    if (logical <= center) {
        return center * (1.0
            - std::sinh(strength * (center - logical) / center) / denominator);
    }
    return center + (1.0 - center)
        * std::sinh(strength * (logical - center) / (1.0 - center))
        / denominator;
}

double wall_clustered_unit_coordinate(double logical, double strength)
{
    if (strength == 0.0) return logical;
    if (!(strength > 0.0) || !std::isfinite(strength)) {
        throw std::invalid_argument("wall-cluster strength must be finite and nonnegative");
    }
    return 0.5 * (1.0
        + std::tanh(strength * (2.0 * logical - 1.0)) / std::tanh(strength));
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
    const std::vector<cgsize_t>& range,
    BCType_t type = BCFarfield)
{
    int boundary = 0;
    check_cgns(
        cg_boco_write(
            file, base, zone, name.c_str(), type,
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
    bool periodic_x,
    double length_x,
    double length_y,
    double length_z,
    bool omit_reciprocal = false,
    double cluster_x = 0.0,
    double cluster_y = 0.0,
    double cluster_strength = 0.0)
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("dimension must be 2 or 3");
    }
    if (dimension == 2 && cells_k != 1) {
        throw std::invalid_argument("two-dimensional grids require cells_k=1");
    }
    if (!(length_x > 0.0) || !(length_y > 0.0) || !(length_z > 0.0)) {
        throw std::invalid_argument("release grid lengths must be positive");
    }
    const bool clustered = cluster_strength > 0.0;
    if (clustered
        && (!(cluster_x > 0.0 && cluster_x < length_x)
            || !(cluster_y > 0.0 && cluster_y < length_y)
            || !std::isfinite(cluster_strength))) {
        throw std::invalid_argument(
            "clustered rectangle centers must be interior and strength finite");
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
                        const double base_x = clustered
                            ? length_x * clustered_unit_coordinate(
                                xi, cluster_x / length_x, cluster_strength)
                            : length_x * xi;
                        const double base_y = clustered
                            ? length_y * clustered_unit_coordinate(
                                eta, cluster_y / length_y, cluster_strength)
                            : length_y * eta;
                        x[index] = base_x + length_x * warp * std::sin(pi * xi)
                            * std::sin(pi * eta) * envelope;
                        y[index] = base_y;
                        z[index] = length_z * zeta;
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
            if (omit_reciprocal) continue;
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

void generate_periodic_square(
    const std::string& path,
    int cells_i,
    int cells_j,
    double length,
    double x_warp_amplitude = 0.0,
    double y_warp_amplitude = 0.0)
{
    if (cells_i % 2 != 0 || cells_j % 2 != 0
        || !std::isfinite(length) || length <= 0.0) {
        throw std::invalid_argument(
            "periodic-square requires even cell counts and positive finite length");
    }
    if (!std::isfinite(x_warp_amplitude)
        || !std::isfinite(y_warp_amplitude)
        || x_warp_amplitude < 0.0 || y_warp_amplitude < 0.0) {
        throw std::invalid_argument(
            "periodic-square warp amplitudes must be finite and nonnegative");
    }
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double maximum_cross_product = 8.0 * pi * pi
        * x_warp_amplitude * y_warp_amplitude / (length * length);
    if (!(maximum_cross_product < 1.0)) {
        throw std::invalid_argument(
            "periodic-square warp can produce a non-positive mapping Jacobian");
    }
    const int local_i = cells_i / 2;
    const int local_j = cells_j / 2;
    const int ni = local_i + 1;
    const int nj = local_j + 1;
    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &file), "cg_open periodic square");
    try {
        int base = 0;
        check_cgns(
            cg_base_write(file, "WCNSPeriodicSquare", 2, 2, &base),
            "cg_base_write periodic square");
        const auto square_zone_name = [](int i, int j) {
            return "Zone" + std::to_string(j * 2 + i + 1);
        };
        for (int zone_j = 0; zone_j < 2; ++zone_j) {
            for (int zone_i = 0; zone_i < 2; ++zone_i) {
                const auto name = square_zone_name(zone_i, zone_j);
                cgsize_t size[6] = {ni, nj, local_i, local_j, 0, 0};
                int zone = 0;
                check_cgns(
                    cg_zone_write(
                        file, base, name.c_str(), size, Structured, &zone),
                    "cg_zone_write periodic square");
                std::vector<double> x(static_cast<std::size_t>(ni * nj));
                std::vector<double> y(x.size());
                for (int j = 0; j < nj; ++j) {
                    for (int i = 0; i < ni; ++i) {
                        const auto index = static_cast<std::size_t>(j * ni + i);
                        const double xi = length
                            * static_cast<double>(zone_i * local_i + i)
                            / static_cast<double>(cells_i);
                        const double eta = length
                            * static_cast<double>(zone_j * local_j + j)
                            / static_cast<double>(cells_j);
                        // This is the x/y mapping used by the case03 reference
                        // Isentropic_curl.cpp.  Its periodicity also makes the
                        // existing translational CGNS 1-to-1 links exact.
                        x[index] = xi + x_warp_amplitude
                            * std::sin(2.0 * pi * eta / length);
                        y[index] = eta + y_warp_amplitude
                            * std::sin(4.0 * pi * xi / length);
                    }
                }
                int coordinate = 0;
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateX", x.data(), &coordinate),
                    "cg_coord_write periodic square X");
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateY", y.data(), &coordinate),
                    "cg_coord_write periodic square Y");
            }
        }

        for (int zone_j = 0; zone_j < 2; ++zone_j) {
            for (int zone_i = 0; zone_i < 2; ++zone_i) {
                const int zone = zone_j * 2 + zone_i + 1;
                const int other_i = 1 - zone_i;
                const int other_j = 1 - zone_j;
                const auto add = [&](
                    const std::string& name,
                    const std::string& donor,
                    const std::vector<cgsize_t>& range,
                    const std::vector<cgsize_t>& donor_range,
                    std::array<float, 3> translation) {
                    const int connection = write_connection(
                        file, base, zone, name, donor,
                        range, donor_range, 2);
                    if (translation[0] != 0.0F || translation[1] != 0.0F) {
                        std::array<float, 3> center {{0.0F, 0.0F, 0.0F}};
                        std::array<float, 3> angle {{0.0F, 0.0F, 0.0F}};
                        check_cgns(
                            cg_1to1_periodic_write(
                                file, base, zone, connection,
                                center.data(), angle.data(), translation.data()),
                            "cg_1to1_periodic_write periodic square");
                    }
                };
                add(
                    "imin", square_zone_name(other_i, zone_j),
                    {1, 1, 1, nj},
                    {ni, 1, ni, nj},
                    {{zone_i == 0 ? static_cast<float>(length) : 0.0F, 0.0F, 0.0F}});
                add(
                    "imax", square_zone_name(other_i, zone_j),
                    {ni, 1, ni, nj},
                    {1, 1, 1, nj},
                    {{zone_i == 1 ? -static_cast<float>(length) : 0.0F, 0.0F, 0.0F}});
                add(
                    "jmin", square_zone_name(zone_i, other_j),
                    {1, 1, ni, 1},
                    {1, nj, ni, nj},
                    {{0.0F, zone_j == 0 ? static_cast<float>(length) : 0.0F, 0.0F}});
                add(
                    "jmax", square_zone_name(zone_i, other_j),
                    {1, nj, ni, nj},
                    {1, 1, ni, 1},
                    {{0.0F, zone_j == 1 ? -static_cast<float>(length) : 0.0F, 0.0F}});
            }
        }
        check_cgns(cg_close(file), "cg_close periodic square");
        file = 0;
    } catch (...) {
        if (file != 0) cg_close(file);
        throw;
    }
}

void generate_periodic_channel(
    const std::string& path,
    int cells_i,
    int cells_j,
    int cells_k,
    int zones_i,
    int zones_k,
    double length_x,
    double length_y,
    double length_z,
    double wall_cluster_strength,
    double origin_y = 0.0)
{
    if (cells_i % zones_i != 0 || cells_k % zones_k != 0
        || zones_i < 2 || zones_k < 2
        || !(length_x > 0.0) || !(length_y > 0.0) || !(length_z > 0.0)
        || !std::isfinite(wall_cluster_strength) || wall_cluster_strength < 0.0
        || !std::isfinite(origin_y)) {
        throw std::invalid_argument(
            "periodic-channel requires divisible cell counts, at least two zones "
            "in each periodic direction, positive lengths and nonnegative clustering");
    }
    const int local_i = cells_i / zones_i;
    const int local_k = cells_k / zones_k;
    const int ni = local_i + 1;
    const int nj = cells_j + 1;
    const int nk = local_k + 1;
    const auto channel_zone_name = [zones_i](int zone_i, int zone_k) {
        return "Zone" + std::to_string(zone_k * zones_i + zone_i + 1);
    };

    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &file), "cg_open periodic channel");
    try {
        int base = 0;
        check_cgns(
            cg_base_write(file, "WCNSPeriodicChannel", 3, 3, &base),
            "cg_base_write periodic channel");
        for (int zone_k = 0; zone_k < zones_k; ++zone_k) {
            for (int zone_i = 0; zone_i < zones_i; ++zone_i) {
                const auto name = channel_zone_name(zone_i, zone_k);
                cgsize_t size[9] = {
                    ni, nj, nk, local_i, cells_j, local_k, 0, 0, 0,
                };
                int zone = 0;
                check_cgns(
                    cg_zone_write(
                        file, base, name.c_str(), size, Structured, &zone),
                    "cg_zone_write periodic channel");
                const std::size_t count = static_cast<std::size_t>(ni)
                    * static_cast<std::size_t>(nj) * static_cast<std::size_t>(nk);
                std::vector<double> x(count);
                std::vector<double> y(count);
                std::vector<double> z(count);
                for (int k = 0; k < nk; ++k) {
                    const int global_k = zone_k * local_k + k;
                    const double zeta
                        = static_cast<double>(global_k) / static_cast<double>(cells_k);
                    for (int j = 0; j < nj; ++j) {
                        const double eta
                            = static_cast<double>(j) / static_cast<double>(cells_j);
                        for (int i = 0; i < ni; ++i) {
                            const int global_i = zone_i * local_i + i;
                            const double xi
                                = static_cast<double>(global_i) / static_cast<double>(cells_i);
                            const auto index
                                = static_cast<std::size_t>((k * nj + j) * ni + i);
                            x[index] = length_x * xi;
                            y[index] = origin_y + length_y
                                * wall_clustered_unit_coordinate(eta, wall_cluster_strength);
                            z[index] = length_z * zeta;
                        }
                    }
                }
                int coordinate = 0;
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateX", x.data(), &coordinate),
                    "cg_coord_write periodic channel X");
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateY", y.data(), &coordinate),
                    "cg_coord_write periodic channel Y");
                check_cgns(
                    cg_coord_write(
                        file, base, zone, RealDouble,
                        "CoordinateZ", z.data(), &coordinate),
                    "cg_coord_write periodic channel Z");
                write_boundary(
                    file, base, zone, "bottom", {1, 1, 1, ni, 1, nk});
                write_boundary(
                    file, base, zone, "top", {1, nj, 1, ni, nj, nk});
            }
        }

        for (int zone_k = 0; zone_k < zones_k; ++zone_k) {
            for (int zone_i = 0; zone_i < zones_i; ++zone_i) {
                const int zone = zone_k * zones_i + zone_i + 1;
                const int previous_i = (zone_i + zones_i - 1) % zones_i;
                const int next_i = (zone_i + 1) % zones_i;
                const int previous_k = (zone_k + zones_k - 1) % zones_k;
                const int next_k = (zone_k + 1) % zones_k;
                const auto add = [&] (
                    const std::string& name,
                    const std::string& donor,
                    const std::vector<cgsize_t>& range,
                    const std::vector<cgsize_t>& donor_range,
                    std::array<float, 3> translation) {
                    const int connection = write_connection(
                        file, base, zone, name, donor,
                        range, donor_range, 3);
                    if (translation[0] != 0.0F || translation[2] != 0.0F) {
                        std::array<float, 3> center {{0.0F, 0.0F, 0.0F}};
                        std::array<float, 3> angle {{0.0F, 0.0F, 0.0F}};
                        check_cgns(
                            cg_1to1_periodic_write(
                                file, base, zone, connection,
                                center.data(), angle.data(), translation.data()),
                            "cg_1to1_periodic_write periodic channel");
                    }
                };
                add(
                    "imin", channel_zone_name(previous_i, zone_k),
                    {1, 1, 1, 1, nj, nk},
                    {ni, 1, 1, ni, nj, nk},
                    {{zone_i == 0 ? static_cast<float>(length_x) : 0.0F, 0.0F, 0.0F}});
                add(
                    "imax", channel_zone_name(next_i, zone_k),
                    {ni, 1, 1, ni, nj, nk},
                    {1, 1, 1, 1, nj, nk},
                    {{zone_i + 1 == zones_i ? -static_cast<float>(length_x) : 0.0F,
                      0.0F, 0.0F}});
                add(
                    "kmin", channel_zone_name(zone_i, previous_k),
                    {1, 1, 1, ni, nj, 1},
                    {1, 1, nk, ni, nj, nk},
                    {{0.0F, 0.0F,
                      zone_k == 0 ? static_cast<float>(length_z) : 0.0F}});
                add(
                    "kmax", channel_zone_name(zone_i, next_k),
                    {1, 1, nk, ni, nj, nk},
                    {1, 1, 1, ni, nj, 1},
                    {{0.0F, 0.0F,
                      zone_k + 1 == zones_k ? -static_cast<float>(length_z) : 0.0F}});
            }
        }
        check_cgns(cg_close(file), "cg_close periodic channel");
        file = 0;
    } catch (...) {
        if (file != 0) cg_close(file);
        throw;
    }
}

double cylinder_radial_coordinate(double logical, double strength)
{
    if (strength == 0.0) return logical;
    if (!std::isfinite(strength) || strength < 0.0) {
        throw std::invalid_argument(
            "cylinder radial-cluster strength must be finite and nonnegative");
    }
    return std::expm1(strength * logical) / std::expm1(strength);
}

void generate_cylinder_o_grid(
    const std::string& path,
    int cells_theta,
    int cells_radial,
    int zones_theta,
    double diameter,
    double outer_radius,
    double radial_cluster_strength)
{
    if (zones_theta < 2 || cells_theta % zones_theta != 0
        || cells_theta / zones_theta < 1 || cells_radial < 1
        || !std::isfinite(diameter) || diameter <= 0.0
        || !std::isfinite(outer_radius) || outer_radius <= 0.5 * diameter
        || !std::isfinite(radial_cluster_strength)
        || radial_cluster_strength < 0.0) {
        throw std::invalid_argument(
            "cylinder-o requires at least two equal angular zones, positive cell "
            "counts and diameter, outer_radius > diameter/2, and nonnegative "
            "radial clustering");
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    const double cylinder_radius = 0.5 * diameter;
    const int local_cells_theta = cells_theta / zones_theta;
    const int ni = local_cells_theta + 1;
    const int nj = cells_radial + 1;
    const auto cylinder_zone_name = [](int zone_index) {
        return "CylinderZone" + std::to_string(zone_index + 1);
    };

    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &file), "cg_open cylinder O-grid");
    try {
        int base = 0;
        check_cgns(
            cg_base_write(file, "WCNSCylinderOGrid", 2, 2, &base),
            "cg_base_write cylinder O-grid");
        for (int zone_index = 0; zone_index < zones_theta; ++zone_index) {
            cgsize_t size[6] = {
                static_cast<cgsize_t>(ni), static_cast<cgsize_t>(nj),
                static_cast<cgsize_t>(local_cells_theta),
                static_cast<cgsize_t>(cells_radial), 0, 0,
            };
            int zone = 0;
            const auto name = cylinder_zone_name(zone_index);
            check_cgns(
                cg_zone_write(
                    file, base, name.c_str(), size, Structured, &zone),
                "cg_zone_write cylinder O-grid");

            const auto count = static_cast<std::size_t>(ni)
                * static_cast<std::size_t>(nj);
            std::vector<double> x(count);
            std::vector<double> y(count);
            for (int j = 0; j < nj; ++j) {
                const double eta
                    = static_cast<double>(j) / static_cast<double>(cells_radial);
                const double radial_fraction = cylinder_radial_coordinate(
                    eta, radial_cluster_strength);
                const double radius = cylinder_radius
                    + (outer_radius - cylinder_radius) * radial_fraction;
                for (int i = 0; i < ni; ++i) {
                    const int global_i = zone_index * local_cells_theta + i;
                    // Clockwise theta makes J=partial(x,y)/partial(i,j) positive
                    // while the radial index grows from the cylinder to farfield.
                    const double theta = -2.0 * pi
                        * static_cast<double>(global_i)
                        / static_cast<double>(cells_theta);
                    const auto index = static_cast<std::size_t>(j * ni + i);
                    x[index] = radius * std::cos(theta);
                    y[index] = radius * std::sin(theta);
                }
            }
            int coordinate = 0;
            check_cgns(
                cg_coord_write(
                    file, base, zone, RealDouble,
                    "CoordinateX", x.data(), &coordinate),
                "cg_coord_write cylinder O-grid X");
            check_cgns(
                cg_coord_write(
                    file, base, zone, RealDouble,
                    "CoordinateY", y.data(), &coordinate),
                "cg_coord_write cylinder O-grid Y");

            write_boundary(
                file, base, zone, "cylinder", {1, 1, ni, 1}, BCWall);
            write_boundary(
                file, base, zone, "farfield", {1, nj, ni, nj}, BCFarfield);
        }

        for (int zone_index = 0; zone_index < zones_theta; ++zone_index) {
            const int zone = zone_index + 1;
            const int previous = (zone_index + zones_theta - 1) % zones_theta;
            const int next = (zone_index + 1) % zones_theta;
            write_connection(
                file, base, zone, "theta-min", cylinder_zone_name(previous),
                {1, 1, 1, nj}, {ni, 1, ni, nj}, 2);
            write_connection(
                file, base, zone, "theta-max", cylinder_zone_name(next),
                {ni, 1, ni, nj}, {1, 1, 1, nj}, 2);
        }
        check_cgns(cg_close(file), "cg_close cylinder O-grid");
        file = 0;
    } catch (...) {
        if (file != 0) cg_close(file);
        throw;
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 13 && argc != 12 && argc != 9 && argc != 8 && argc != 6 && argc != 5) {
        std::cerr
            << "usage: wcns_generate_release_cgns <output.cgns> <dimension> "
               "<cells_i> <cells_j> <cells_k> <zones_i> <warp> <periodic_x>\n"
               "   or: wcns_generate_release_cgns periodic-square <output.cgns> "
               "<cells_i> <cells_j> <length>\n"
               "   or: wcns_generate_release_cgns warped-periodic-square "
               "<output.cgns> <cells_i> <cells_j> <length> "
               "<x_warp_amplitude> <y_warp_amplitude>\n"
               "   or: wcns_generate_release_cgns rectangle <output.cgns> "
               "<cells_i> <cells_j> <zones_i> <length_x> <length_y> "
               "<periodic_x>\n"
               "   or: wcns_generate_release_cgns clustered-rectangle "
               "<output.cgns> <cells_i> <cells_j> <zones_i> "
               "<length_x> <length_y> <cluster_x> <cluster_y> "
               "<strength> <periodic_x>\n"
               "   or: wcns_generate_release_cgns periodic-channel "
               "<output.cgns> <cells_i> <cells_j> <cells_k> "
               "<zones_i> <zones_k> <length_x> <length_y> <length_z> "
               "<wall_cluster_strength> [origin_y]\n"
               "   or: wcns_generate_release_cgns cylinder-o <output.cgns> "
               "<cells_theta> <cells_radial> <zones_theta> <diameter> "
               "<outer_radius> <radial_cluster_strength>\n"
               "   or: wcns_generate_release_cgns invalid-one-sided "
               "<output.cgns> <cells_i> <cells_j>\n";
        return EXIT_FAILURE;
    }
    try {
        check_cgns(cg_set_file_type(CG_FILE_ADF), "cg_set_file_type release grid");
        std::string output;
        if (argc == 5 && std::string(argv[1]) == "invalid-one-sided") {
            output = argv[2];
            generate(
                output, 2,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"), 1, 2, 0.0, false,
                1.0, 1.0, 1.0, true);
        } else if (argc == 6 && std::string(argv[1]) == "periodic-square") {
            output = argv[2];
            generate_periodic_square(
                output,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"),
                parse_positive_real(argv[5], "length"));
        } else if (argc == 8
                   && std::string(argv[1]) == "warped-periodic-square") {
            output = argv[2];
            generate_periodic_square(
                output,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"),
                parse_positive_real(argv[5], "length"),
                parse_nonnegative_real(argv[6], "x_warp_amplitude"),
                parse_nonnegative_real(argv[7], "y_warp_amplitude"));
        } else if ((argc == 12 || argc == 13)
                   && std::string(argv[1]) == "periodic-channel") {
            output = argv[2];
            generate_periodic_channel(
                output,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"),
                parse_positive(argv[5], "cells_k"),
                parse_positive(argv[6], "zones_i"),
                parse_positive(argv[7], "zones_k"),
                parse_positive_real(argv[8], "length_x"),
                parse_positive_real(argv[9], "length_y"),
                parse_positive_real(argv[10], "length_z"),
                parse_nonnegative_real(argv[11], "wall_cluster_strength"),
                argc == 13 ? parse_finite_real(argv[12], "origin_y") : 0.0);
        } else if (argc == 12 && std::string(argv[1]) == "clustered-rectangle") {
            output = argv[2];
            const double length_x = parse_positive_real(argv[6], "length_x");
            const double length_y = parse_positive_real(argv[7], "length_y");
            generate(
                output, 2,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"), 1,
                parse_positive(argv[5], "zones_i"), 0.0,
                parse_bool(argv[11]), length_x, length_y, 1.0, false,
                parse_positive_real(argv[8], "cluster_x"),
                parse_positive_real(argv[9], "cluster_y"),
                parse_positive_real(argv[10], "strength"));
        } else if (argc == 9 && std::string(argv[1]) == "cylinder-o") {
            output = argv[2];
            generate_cylinder_o_grid(
                output,
                parse_positive(argv[3], "cells_theta"),
                parse_positive(argv[4], "cells_radial"),
                parse_positive(argv[5], "zones_theta"),
                parse_positive_real(argv[6], "diameter"),
                parse_positive_real(argv[7], "outer_radius"),
                parse_nonnegative_real(argv[8], "radial_cluster_strength"));
        } else if (argc == 9 && std::string(argv[1]) == "rectangle") {
            output = argv[2];
            generate(
                output, 2,
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"), 1,
                parse_positive(argv[5], "zones_i"), 0.0,
                parse_bool(argv[8]),
                parse_positive_real(argv[6], "length_x"),
                parse_positive_real(argv[7], "length_y"), 1.0);
        } else if (argc == 9) {
            output = argv[1];
            generate(
                output, parse_positive(argv[2], "dimension"),
                parse_positive(argv[3], "cells_i"),
                parse_positive(argv[4], "cells_j"),
                parse_positive(argv[5], "cells_k"),
                parse_positive(argv[6], "zones_i"),
                parse_warp(argv[7]), parse_bool(argv[8]),
                1.0, 1.0, 1.0);
        } else {
            throw std::invalid_argument("invalid release grid generator mode");
        }
        std::cout << "generated release CGNS grid: " << output << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "release grid generation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
