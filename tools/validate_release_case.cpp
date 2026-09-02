#include <cgnslib.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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
            for (int axis = 0; axis < dimension; ++axis) {
                values.cells.push_back(size[dimension + axis]);
                count *= static_cast<std::size_t>(size[dimension + axis]);
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

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 3 && std::string(argv[1]) == "finite") {
            validate_finite(argv[2]);
        } else if (argc == 5 && std::string(argv[1]) == "compare") {
            compare_fields(argv[2], argv[3], parse_real(argv[4], "tolerance"));
        } else if (argc == 9 && std::string(argv[1]) == "uniform") {
            validate_uniform(
                argv[2],
                {{"Density", parse_real(argv[3], "rho")},
                 {"VelocityX", parse_real(argv[4], "u")},
                 {"VelocityY", parse_real(argv[5], "v")},
                 {"VelocityZ", parse_real(argv[6], "w")},
                 {"Temperature", parse_real(argv[7], "temperature")}},
                parse_real(argv[8], "tolerance"));
        } else {
            std::cerr
                << "usage:\n"
                   "  wcns_validate_release_case finite <field.cgns>\n"
                   "  wcns_validate_release_case compare <lhs.cgns> <rhs.cgns> <tol>\n"
                   "  wcns_validate_release_case uniform <field.cgns> "
                   "<rho> <u> <v> <w> <T> <tol>\n";
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "release validation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
