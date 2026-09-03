#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
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

double parse_real(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
    return value;
}

struct ZoneFields {
    int dimension = 0;
    std::vector<cgsize_t> cells;
    std::vector<std::array<double, 3>> cell_centers;
    std::map<std::string, std::vector<double>> fields;
};

using FieldFile = std::map<std::string, ZoneFields>;

FieldFile read_fields(const std::string& path)
{
    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_READ, &file), "cg_open release field");
    try {
        int bases = 0;
        check_cgns(cg_nbases(file, &bases), "cg_nbases release field");
        if (bases != 1) throw std::runtime_error("release field requires one CGNS base");
        char base_name[33] = {};
        int dimension = 0;
        int physical_dimension = 0;
        check_cgns(
            cg_base_read(file, 1, base_name, &dimension, &physical_dimension),
            "cg_base_read release field");
        if (dimension != 2 && dimension != 3) {
            throw std::runtime_error("release field dimension must be two or three");
        }
        int zone_count = 0;
        check_cgns(cg_nzones(file, 1, &zone_count), "cg_nzones release field");
        FieldFile result;
        for (int zone = 1; zone <= zone_count; ++zone) {
            char zone_name[33] = {};
            cgsize_t size[9] = {};
            check_cgns(
                cg_zone_read(file, 1, zone, zone_name, size),
                "cg_zone_read release field");
            ZoneFields values;
            values.dimension = dimension;
            std::size_t count = 1;
            std::vector<cgsize_t> vertices;
            for (int axis = 0; axis < dimension; ++axis) {
                vertices.push_back(size[axis]);
                values.cells.push_back(size[dimension + axis]);
                count *= static_cast<std::size_t>(size[dimension + axis]);
            }
            std::size_t vertex_count = 1;
            for (const auto extent : vertices) {
                vertex_count *= static_cast<std::size_t>(extent);
            }
            std::array<std::vector<double>, 3> coordinates;
            for (auto& coordinate : coordinates) coordinate.assign(vertex_count, 0.0);
            const std::array<const char*, 3> coordinate_names {{
                "CoordinateX", "CoordinateY", "CoordinateZ",
            }};
            std::vector<cgsize_t> vertex_lower(static_cast<std::size_t>(dimension), 1);
            for (int axis = 0; axis < dimension; ++axis) {
                check_cgns(
                    cg_coord_read(
                        file, 1, zone, coordinate_names[static_cast<std::size_t>(axis)],
                        RealDouble, vertex_lower.data(), vertices.data(),
                        coordinates[static_cast<std::size_t>(axis)].data()),
                    "cg_coord_read release field");
            }
            values.cell_centers.resize(count);
            const auto vertex_index = [&](int i, int j, int k) {
                const std::size_t ni = static_cast<std::size_t>(vertices[0]);
                const std::size_t nj = static_cast<std::size_t>(vertices[1]);
                return (static_cast<std::size_t>(k) * nj
                    + static_cast<std::size_t>(j)) * ni
                    + static_cast<std::size_t>(i);
            };
            const int cell_nk = dimension == 3 ? static_cast<int>(values.cells[2]) : 1;
            std::size_t cell = 0;
            for (int k = 0; k < cell_nk; ++k) {
                for (int j = 0; j < static_cast<int>(values.cells[1]); ++j) {
                    for (int i = 0; i < static_cast<int>(values.cells[0]); ++i) {
                        std::array<double, 3> center {{0.0, 0.0, 0.0}};
                        const int corner_count = dimension == 3 ? 8 : 4;
                        for (int corner = 0; corner < corner_count; ++corner) {
                            const int di = corner & 1;
                            const int dj = (corner >> 1) & 1;
                            const int dk = dimension == 3 ? ((corner >> 2) & 1) : 0;
                            const auto vertex = vertex_index(i + di, j + dj, k + dk);
                            for (int axis = 0; axis < dimension; ++axis) {
                                center[static_cast<std::size_t>(axis)]
                                    += coordinates[static_cast<std::size_t>(axis)][vertex]
                                    / static_cast<double>(corner_count);
                            }
                        }
                        values.cell_centers[cell++] = center;
                    }
                }
            }
            int solutions = 0;
            check_cgns(
                cg_nsols(file, 1, zone, &solutions),
                "cg_nsols release field");
            if (solutions != 1) {
                throw std::runtime_error("release field zone requires one solution");
            }
            int field_count = 0;
            check_cgns(
                cg_nfields(file, 1, zone, 1, &field_count),
                "cg_nfields release field");
            std::vector<cgsize_t> lower(static_cast<std::size_t>(dimension), 1);
            std::vector<cgsize_t> upper = values.cells;
            for (int field = 1; field <= field_count; ++field) {
                DataType_t data_type = DataTypeNull;
                char field_name[33] = {};
                check_cgns(
                    cg_field_info(file, 1, zone, 1, field, &data_type, field_name),
                    "cg_field_info release field");
                std::vector<double> field_values(count);
                check_cgns(
                    cg_field_read(
                        file, 1, zone, 1, field_name, RealDouble,
                        lower.data(), upper.data(), field_values.data()),
                    "cg_field_read release field");
                values.fields.emplace(field_name, std::move(field_values));
            }
            if (!result.emplace(zone_name, std::move(values)).second) {
                throw std::runtime_error("duplicate zone name in release field");
            }
        }
        check_cgns(cg_close(file), "cg_close release field");
        return result;
    } catch (...) {
        cg_close(file);
        throw;
    }
}

const std::vector<double>& require_field(
    const ZoneFields& zone,
    const std::string& name)
{
    const auto iterator = zone.fields.find(name);
    if (iterator == zone.fields.end()) {
        throw std::runtime_error("release output is missing field " + name);
    }
    return iterator->second;
}

void require_tolerance(double tolerance)
{
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance)) {
        throw std::invalid_argument("tolerance must be finite and nonnegative");
    }
}

double validate_uniform(
    const std::string& path,
    const std::map<std::string, double>& expected,
    double tolerance)
{
    require_tolerance(tolerance);
    const auto file = read_fields(path);
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const auto& [field_name, reference] : expected) {
            for (const double value : require_field(zone, field_name)) {
                if (!std::isfinite(value)) {
                    throw std::runtime_error("uniform field contains a non-finite value");
                }
                maximum = std::max(maximum, std::abs(value - reference));
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error(
            "uniform maximum error " + std::to_string(maximum)
            + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17)
              << "check=uniform samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
    return maximum;
}

void validate_finite(const std::string& path)
{
    const auto file = read_fields(path);
    std::size_t samples = 0;
    double minimum_density = std::numeric_limits<double>::infinity();
    double minimum_pressure = std::numeric_limits<double>::infinity();
    double minimum_temperature = std::numeric_limits<double>::infinity();
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const auto& [field_name, values] : zone.fields) {
            for (const double value : values) {
                if (!std::isfinite(value)) {
                    throw std::runtime_error("release field contains a non-finite value");
                }
                ++samples;
                if (field_name == "Density") minimum_density = std::min(minimum_density, value);
                if (field_name == "Pressure") minimum_pressure = std::min(minimum_pressure, value);
                if (field_name == "Temperature") {
                    minimum_temperature = std::min(minimum_temperature, value);
                }
            }
        }
    }
    for (const auto& value : {
             std::pair<const char*, double> {"density", minimum_density},
             {"pressure", minimum_pressure}, {"temperature", minimum_temperature}}) {
        if (std::isfinite(value.second) && !(value.second > 0.0)) {
            throw std::runtime_error(std::string(value.first) + " is not positive");
        }
    }
    std::cout << std::setprecision(17)
              << "check=finite samples=" << samples
              << " min_rho=" << minimum_density
              << " min_p=" << minimum_pressure
              << " min_T=" << minimum_temperature << '\n';
}

void compare_fields(
    const std::string& lhs_path,
    const std::string& rhs_path,
    double tolerance)
{
    require_tolerance(tolerance);
    const auto lhs = read_fields(lhs_path);
    const auto rhs = read_fields(rhs_path);
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("release field comparison zone counts differ");
    }
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [name, lhs_zone] : lhs) {
        const auto rhs_iterator = rhs.find(name);
        if (rhs_iterator == rhs.end()) {
            throw std::runtime_error("release field comparison zone names differ");
        }
        const auto& rhs_zone = rhs_iterator->second;
        if (lhs_zone.dimension != rhs_zone.dimension
            || lhs_zone.cells != rhs_zone.cells
            || lhs_zone.fields.size() != rhs_zone.fields.size()) {
            throw std::runtime_error("release field comparison metadata differs");
        }
        for (const auto& [field_name, lhs_values] : lhs_zone.fields) {
            const auto& rhs_values = require_field(rhs_zone, field_name);
            if (lhs_values.size() != rhs_values.size()) {
                throw std::runtime_error("release field comparison array sizes differ");
            }
            for (std::size_t index = 0; index < lhs_values.size(); ++index) {
                if (!std::isfinite(lhs_values[index])
                    || !std::isfinite(rhs_values[index])) {
                    throw std::runtime_error("release field comparison is non-finite");
                }
                maximum = std::max(
                    maximum, std::abs(lhs_values[index] - rhs_values[index]));
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error(
            "field comparison maximum difference " + std::to_string(maximum)
            + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17)
              << "check=compare samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
}

struct SpatialValue {
    std::array<double, 3> coordinates;
    std::vector<double> fields;
};

std::vector<SpatialValue> spatial_values(
    const FieldFile& file,
    const std::vector<std::string>& field_names)
{
    std::vector<SpatialValue> result;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            SpatialValue value;
            value.coordinates = zone.cell_centers[cell];
            value.fields.reserve(field_names.size());
            for (const auto& name : field_names) {
                value.fields.push_back(require_field(zone, name).at(cell));
            }
            result.push_back(std::move(value));
        }
    }
    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.coordinates < rhs.coordinates;
    });
    return result;
}

void compare_spatial_fields(
    const std::string& lhs_path,
    const std::string& rhs_path,
    double field_tolerance,
    double coordinate_tolerance)
{
    require_tolerance(field_tolerance);
    require_tolerance(coordinate_tolerance);
    const auto lhs = read_fields(lhs_path);
    const auto rhs = read_fields(rhs_path);
    if (lhs.empty() || rhs.empty()) {
        throw std::runtime_error("spatial comparison requires nonempty field files");
    }
    std::vector<std::string> fields;
    for (const auto& [name, values] : lhs.begin()->second.fields) {
        static_cast<void>(values);
        fields.push_back(name);
    }
    const auto check_field_set = [&](const FieldFile& file) {
        for (const auto& [zone_name, zone] : file) {
            static_cast<void>(zone_name);
            std::vector<std::string> names;
            for (const auto& [name, values] : zone.fields) {
                static_cast<void>(values);
                names.push_back(name);
            }
            if (names != fields) {
                throw std::runtime_error("spatial comparison field sets differ");
            }
        }
    };
    check_field_set(lhs);
    check_field_set(rhs);
    const auto lhs_values = spatial_values(lhs, fields);
    const auto rhs_values = spatial_values(rhs, fields);
    if (lhs_values.size() != rhs_values.size()) {
        throw std::runtime_error("spatial comparison cell counts differ");
    }
    double maximum_coordinate = 0.0;
    double maximum_field = 0.0;
    for (std::size_t cell = 0; cell < lhs_values.size(); ++cell) {
        for (std::size_t component = 0; component < 3; ++component) {
            maximum_coordinate = std::max(
                maximum_coordinate,
                std::abs(lhs_values[cell].coordinates[component]
                    - rhs_values[cell].coordinates[component]));
        }
        for (std::size_t field = 0; field < fields.size(); ++field) {
            if (!std::isfinite(lhs_values[cell].fields[field])
                || !std::isfinite(rhs_values[cell].fields[field])) {
                throw std::runtime_error("spatial comparison field is non-finite");
            }
            maximum_field = std::max(
                maximum_field,
                std::abs(lhs_values[cell].fields[field]
                    - rhs_values[cell].fields[field]));
        }
    }
    if (maximum_coordinate > coordinate_tolerance
        || maximum_field > field_tolerance) {
        throw std::runtime_error("spatial field comparison exceeds tolerance");
    }
    std::cout << std::setprecision(17)
              << "check=compare_spatial cells=" << lhs_values.size()
              << " max_coordinate=" << maximum_coordinate
              << " coordinate_tolerance=" << coordinate_tolerance
              << " max_field=" << maximum_field
              << " field_tolerance=" << field_tolerance << '\n';
}

void validate_constant_series(const std::string& path, double tolerance)
{
    require_tolerance(tolerance);
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open release series: " + path);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::istringstream values(line);
        std::vector<double> row;
        double value = 0.0;
        while (values >> value) row.push_back(value);
        if (!values.eof() || row.size() < 3) {
            throw std::runtime_error("release series contains an invalid data row");
        }
        rows.push_back(std::move(row));
    }
    if (rows.size() < 2) {
        throw std::runtime_error("release series requires at least two samples");
    }
    double maximum = 0.0;
    for (const auto& row : rows) {
        if (row.size() != rows.front().size()) {
            throw std::runtime_error("release series row widths differ");
        }
        for (std::size_t column = 2; column < row.size(); ++column) {
            if (!std::isfinite(row[column])) {
                throw std::runtime_error("release series contains a non-finite value");
            }
            maximum = std::max(
                maximum,
                std::abs(row[column] - rows.front()[column])
                    / std::max(1.0, std::abs(rows.front()[column])));
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error(
            "series relative drift " + std::to_string(maximum)
            + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17)
              << "check=series_constant rows=" << rows.size()
              << " max_relative_drift=" << maximum
              << " tolerance=" << tolerance << '\n';
}

void validate_isentropic_vortex(
    const std::string& path,
    double time,
    double length,
    double x0,
    double y0,
    double beta,
    double background_u,
    double background_v,
    double gamma,
    double mach,
    double l1_tolerance)
{
    if (!(length > 0.0) || !(gamma > 1.0) || !(mach > 0.0)) {
        throw std::invalid_argument("vortex length, gamma and Mach are invalid");
    }
    require_tolerance(l1_tolerance);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double center_x = std::fmod(x0 + background_u * time, length);
    const double center_y = std::fmod(y0 + background_v * time, length);
    const auto file = read_fields(path);
    double density_l1 = 0.0;
    double density_l2 = 0.0;
    double density_linf = 0.0;
    double maximum_primitive = 0.0;
    std::size_t cells = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& density = require_field(zone, "Density");
        const auto& velocity_x = require_field(zone, "VelocityX");
        const auto& velocity_y = require_field(zone, "VelocityY");
        const auto& temperature = require_field(zone, "Temperature");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            const double dx = std::remainder(
                zone.cell_centers[cell][0] - center_x, length);
            const double dy = std::remainder(
                zone.cell_centers[cell][1] - center_y, length);
            const double radius_squared = dx * dx + dy * dy;
            const double exponential = std::exp(0.5 * (1.0 - radius_squared));
            const double exact_u = background_u
                - beta * exponential * dy / (2.0 * pi);
            const double exact_v = background_v
                + beta * exponential * dx / (2.0 * pi);
            const double exact_temperature = 1.0
                - (gamma - 1.0) * mach * mach * beta * beta
                    * std::exp(1.0 - radius_squared)
                    / (8.0 * pi * pi);
            const double exact_density
                = std::pow(exact_temperature, 1.0 / (gamma - 1.0));
            const double density_error = std::abs(density[cell] - exact_density);
            density_l1 += density_error;
            density_l2 += density_error * density_error;
            density_linf = std::max(density_linf, density_error);
            maximum_primitive = std::max({
                maximum_primitive,
                density_error,
                std::abs(velocity_x[cell] - exact_u),
                std::abs(velocity_y[cell] - exact_v),
                std::abs(temperature[cell] - exact_temperature),
            });
            ++cells;
        }
    }
    if (cells == 0) throw std::runtime_error("vortex output has no cells");
    density_l1 /= static_cast<double>(cells);
    density_l2 = std::sqrt(density_l2 / static_cast<double>(cells));
    if (density_l1 > l1_tolerance) {
        throw std::runtime_error(
            "vortex density L1 error " + std::to_string(density_l1)
            + " exceeds tolerance " + std::to_string(l1_tolerance));
    }
    std::cout << std::setprecision(17)
              << "check=isentropic_vortex cells=" << cells
              << " rho_l1=" << density_l1
              << " rho_l2=" << density_l2
              << " rho_linf=" << density_linf
              << " primitive_linf=" << maximum_primitive
              << " tolerance=" << l1_tolerance << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 3 && std::string(argv[1]) == "finite") {
            validate_finite(argv[2]);
        } else if (argc == 4 && std::string(argv[1]) == "series-constant") {
            validate_constant_series(
                argv[2], parse_real(argv[3], "tolerance"));
        } else if (argc == 5 && std::string(argv[1]) == "compare") {
            compare_fields(argv[2], argv[3], parse_real(argv[4], "tolerance"));
        } else if (argc == 6 && std::string(argv[1]) == "compare-spatial") {
            compare_spatial_fields(
                argv[2], argv[3],
                parse_real(argv[4], "field tolerance"),
                parse_real(argv[5], "coordinate tolerance"));
        } else if (argc == 9 && std::string(argv[1]) == "uniform") {
            validate_uniform(
                argv[2],
                {{"Density", parse_real(argv[3], "rho")},
                 {"VelocityX", parse_real(argv[4], "u")},
                 {"VelocityY", parse_real(argv[5], "v")},
                 {"VelocityZ", parse_real(argv[6], "w")},
                 {"Temperature", parse_real(argv[7], "temperature")}},
                parse_real(argv[8], "tolerance"));
        } else if (argc == 13 && std::string(argv[1]) == "vortex") {
            validate_isentropic_vortex(
                argv[2],
                parse_real(argv[3], "time"),
                parse_real(argv[4], "length"),
                parse_real(argv[5], "x0"),
                parse_real(argv[6], "y0"),
                parse_real(argv[7], "beta"),
                parse_real(argv[8], "background_u"),
                parse_real(argv[9], "background_v"),
                parse_real(argv[10], "gamma"),
                parse_real(argv[11], "Mach"),
                parse_real(argv[12], "tolerance"));
        } else {
            std::cerr
                << "usage:\n"
                   "  wcns_validate_release_case finite <field.cgns>\n"
                   "  wcns_validate_release_case series-constant <series.txt> <tol>\n"
                   "  wcns_validate_release_case compare <lhs.cgns> <rhs.cgns> <tol>\n"
                   "  wcns_validate_release_case compare-spatial <lhs.cgns> "
                   "<rhs.cgns> <field-tol> <coordinate-tol>\n"
                   "  wcns_validate_release_case uniform <field.cgns> "
                   "<rho> <u> <v> <w> <T> <tol>\n"
                   "  wcns_validate_release_case vortex <field.cgns> <time> "
                   "<length> <x0> <y0> <beta> <u0> <v0> <gamma> <Mach> "
                   "<L1-tol>\n";
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "release validation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
