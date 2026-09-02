#include <wcns/runtime/case_config.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace wcns {
namespace {

using EntryMap = std::unordered_map<std::string, std::string>;

std::string trim(std::string_view text)
{
    std::size_t first = 0;
    while (first < text.size()
        && std::isspace(static_cast<unsigned char>(text[first])) != 0) {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first
        && std::isspace(static_cast<unsigned char>(text[last - 1])) != 0) {
        --last;
    }
    return std::string(text.substr(first, last - first));
}

EntryMap parse_entries(const std::string& text)
{
    EntryMap entries;
    std::istringstream stream(text);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const auto clean = trim(line);
        if (clean.empty() || clean.front() == '#') {
            continue;
        }
        const auto equal = clean.find('=');
        if (equal == std::string::npos || clean.find('=', equal + 1) != std::string::npos) {
            throw CaseConfigurationError(
                "configuration line " + std::to_string(line_number)
                + " must contain exactly one '='");
        }
        const auto key = trim(std::string_view(clean).substr(0, equal));
        const auto value = trim(std::string_view(clean).substr(equal + 1));
        if (key.empty() || value.empty()) {
            throw CaseConfigurationError(
                "configuration line " + std::to_string(line_number)
                + " has an empty key or value");
        }
        if (!entries.emplace(key, value).second) {
            throw CaseConfigurationError("duplicate configuration key: " + key);
        }
    }
    return entries;
}

const std::string& require(const EntryMap& entries, const std::string& key)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) {
        throw CaseConfigurationError("missing required configuration key: " + key);
    }
    return iterator->second;
}

Real parse_real(const std::string& value, const std::string& key)
{
    std::size_t consumed = 0;
    Real result = 0.0;
    try {
        result = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw CaseConfigurationError("configuration key is not a real number: " + key);
    }
    if (consumed != value.size() || !std::isfinite(result)) {
        throw CaseConfigurationError(
            "configuration key requires a finite real number: " + key);
    }
    return result;
}

long long parse_integer(const std::string& value, const std::string& key)
{
    std::size_t consumed = 0;
    long long result = 0;
    try {
        result = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        throw CaseConfigurationError("configuration key is not an integer: " + key);
    }
    if (consumed != value.size()) {
        throw CaseConfigurationError("configuration key is not an integer: " + key);
    }
    return result;
}

bool parse_bool(const std::string& value, const std::string& key)
{
    if (value == "true") return true;
    if (value == "false") return false;
    throw CaseConfigurationError(
        "configuration key requires true or false: " + key);
}

PartitionMode parse_partition_mode(const std::string& value)
{
    if (value == "zones_only") return PartitionMode::ZonesOnly;
    if (value == "auto_split") return PartitionMode::AutoSplit;
    if (value == "force_split") return PartitionMode::ForceSplit;
    throw CaseConfigurationError("unknown partition mode: " + value);
}

RunMode parse_run_mode(const std::string& value)
{
    if (value == "steady") return RunMode::Steady;
    if (value == "unsteady") return RunMode::Unsteady;
    throw CaseConfigurationError("unknown run mode: " + value);
}

FieldOutputFormat parse_field_output_format(const std::string& value)
{
    if (value == "cgns") return FieldOutputFormat::Cgns;
    if (value == "tecplot") return FieldOutputFormat::Tecplot;
    if (value == "both") return FieldOutputFormat::Both;
    throw CaseConfigurationError("unknown field output format: " + value);
}

SeriesOutputFormat parse_series_output_format(const std::string& value)
{
    if (value == "txt") return SeriesOutputFormat::Text;
    if (value == "tecplot") return SeriesOutputFormat::Tecplot;
    throw CaseConfigurationError("unknown series output format: " + value);
}

ReconstructionVariables parse_reconstruction_variables(const std::string& value)
{
    if (value == "conservative") return ReconstructionVariables::Conservative;
    if (value == "primitive") return ReconstructionVariables::Primitive;
    if (value == "characteristic") return ReconstructionVariables::Characteristic;
    throw CaseConfigurationError("unknown reconstruction variable space: " + value);
}

BoundaryType parse_boundary_type(const std::string& value)
{
    if (value == "farfield") return BoundaryType::Farfield;
    if (value == "inflow") return BoundaryType::Inflow;
    if (value == "outflow") return BoundaryType::Outflow;
    if (value == "slip_wall") return BoundaryType::SlipWall;
    if (value == "no_slip_adiabatic_wall") return BoundaryType::NoSlipAdiabaticWall;
    if (value == "no_slip_isothermal_wall") return BoundaryType::NoSlipIsothermalWall;
    if (value == "symmetry") return BoundaryType::Symmetry;
    if (value == "periodic") return BoundaryType::Periodic;
    throw CaseConfigurationError("unknown boundary type: " + value);
}

std::vector<std::string> split_list(const std::string& value)
{
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find(',', begin);
        const auto item = trim(std::string_view(value).substr(
            begin, end == std::string::npos ? value.size() - begin : end - begin));
        if (item.empty()) {
            throw CaseConfigurationError("list contains an empty item");
        }
        result.push_back(item);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

SourceModelKind parse_source_model(const std::string& value)
{
    if (value == "uniform_conservative") return SourceModelKind::UniformConservative;
    if (value == "body_force") return SourceModelKind::BodyForce;
    if (value == "manufactured") return SourceModelKind::ManufacturedSolution;
    throw CaseConfigurationError("unknown source model: " + value);
}

std::uint64_t fnv1a(const std::string& text)
{
    std::uint64_t result = 14695981039346656037ull;
    for (const unsigned char byte : text) {
        result ^= static_cast<std::uint64_t>(byte);
        result *= 1099511628211ull;
    }
    return result;
}

std::string canonical_entries(const EntryMap& entries)
{
    std::vector<std::pair<std::string, std::string>> sorted(
        entries.begin(), entries.end());
    std::sort(sorted.begin(), sorted.end());
    std::ostringstream result;
    for (const auto& [key, value] : sorted) {
        result << key << '=' << value << '\n';
    }
    return result.str();
}

const std::set<std::string>& fixed_keys()
{
    static const std::set<std::string> keys {
        "schema_version", "case.name", "mesh.path",
        "algorithm.profile", "algorithm.reconstruction",
        "algorithm.reconstruction_variables", "algorithm.riemann",
        "gas.gamma", "gas.molar_mass", "gas.specific_gas_constant",
        "reference.velocity", "reference.density", "reference.temperature",
        "reference.length", "reference.viscosity",
        "partition.mode", "partition.allow_idle_ranks",
        "partition.max_load_ratio", "partition.min_cells_per_active_direction",
        "initial.type", "initial.rho", "initial.u", "initial.v", "initial.w",
        "initial.pressure", "initial.temperature", "initial.x0", "initial.y0",
        "initial.beta", "initial.background_u", "initial.background_v",
        "initial.left_rho", "initial.left_u", "initial.left_v", "initial.left_p",
        "initial.right_rho", "initial.right_u", "initial.right_v", "initial.right_p",
        "initial.ne_rho", "initial.ne_u", "initial.ne_v", "initial.ne_p",
        "initial.nw_rho", "initial.nw_u", "initial.nw_v", "initial.nw_p",
        "initial.sw_rho", "initial.sw_u", "initial.sw_v", "initial.sw_p",
        "initial.se_rho", "initial.se_u", "initial.se_v", "initial.se_p",
        "boundary.default", "source.enabled", "source.models",
        "source.uniform.rho", "source.uniform.momentum_x",
        "source.uniform.momentum_y", "source.uniform.momentum_z",
        "source.uniform.energy", "source.body.ax", "source.body.ay",
        "source.body.az", "source.manufactured.rho",
        "source.manufactured.momentum_x", "source.manufactured.momentum_y",
        "source.manufactured.momentum_z", "source.manufactured.energy",
        "run.mode", "run.viscous", "run.cfl", "run.max_steps",
        "run.t_end", "run.max_wall_time",
        "steady.min_steps", "steady.check_interval_steps",
        "steady.consecutive_checks", "steady.reference_floor",
        "steady.l2_absolute", "steady.l2_relative",
        "steady.linf_enabled", "steady.linf_absolute",
        "steady.linf_relative", "output.directory",
        "output.allow_existing", "output.dimensional",
        "output.field.enabled", "output.field.format",
        "output.field.every_steps", "output.field.every_time",
        "output.field.explicit_times", "output.field.write_initial",
        "output.field.write_final", "output.field.quantities",
        "output.history.enabled", "output.history.format",
        "output.history.every_steps", "output.history.every_time",
        "output.history.explicit_times", "output.history.write_initial",
        "output.history.write_final", "output.history.quantities",
        "output.statistics.enabled", "output.statistics.format",
        "output.statistics.every_steps", "output.statistics.every_time",
        "output.statistics.explicit_times", "output.statistics.write_initial",
        "output.statistics.write_final", "output.statistics.quantities",
        "output.checkpoint.enabled", "output.checkpoint.every_steps",
        "output.checkpoint.every_time", "output.checkpoint.explicit_times",
        "output.checkpoint.write_initial", "output.checkpoint.write_final",
        "restart.path",
    };
    return keys;
}

bool dynamic_boundary_key(const std::string& key)
{
    constexpr std::string_view prefix = "boundary.";
    constexpr std::string_view suffix = ".type";
    return key.size() > prefix.size() + suffix.size()
        && key.compare(0, prefix.size(), prefix.data(), prefix.size()) == 0
        && key.compare(
            key.size() - suffix.size(),
            suffix.size(),
            suffix.data(),
            suffix.size()) == 0;
}

void reject_unknown_keys(const EntryMap& entries)
{
    for (const auto& [key, value] : entries) {
        static_cast<void>(value);
        if (key == "Re" || key == "Ma" || key == "reference.reynolds"
            || key == "reference.mach") {
            throw CaseConfigurationError(
                "Re and Ma are derived values and cannot be configured: " + key);
        }
        if (fixed_keys().find(key) == fixed_keys().end()
            && !dynamic_boundary_key(key)) {
            throw CaseConfigurationError("unknown configuration key: " + key);
        }
    }
}

Real optional_real(
    const EntryMap& entries,
    const std::string& key,
    Real default_value)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end()
        ? default_value : parse_real(iterator->second, key);
}

bool optional_bool(
    const EntryMap& entries,
    const std::string& key,
    bool default_value)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end()
        ? default_value : parse_bool(iterator->second, key);
}

std::size_t optional_size(
    const EntryMap& entries,
    const std::string& key,
    std::size_t default_value)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) return default_value;
    const auto value = parse_integer(iterator->second, key);
    if (value < 0
        || static_cast<unsigned long long>(value)
            > std::numeric_limits<std::size_t>::max()) {
        throw CaseConfigurationError(
            "configuration key is outside size_t range: " + key);
    }
    return static_cast<std::size_t>(value);
}

std::vector<Real> optional_real_list(
    const EntryMap& entries,
    const std::string& key)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) return {};
    std::vector<Real> result;
    for (const auto& item : split_list(iterator->second)) {
        result.push_back(parse_real(item, key));
    }
    return result;
}

std::vector<std::string> optional_string_list(
    const EntryMap& entries,
    const std::string& key)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end()
        ? std::vector<std::string> {} : split_list(iterator->second);
}

OutputScheduleConfig parse_schedule(
    const EntryMap& entries,
    const std::string& prefix,
    bool default_final)
{
    OutputScheduleConfig result;
    result.every_steps = optional_size(entries, prefix + ".every_steps", 0);
    result.every_time = optional_real(entries, prefix + ".every_time", 0.0);
    result.explicit_times = optional_real_list(entries, prefix + ".explicit_times");
    result.write_initial = optional_bool(
        entries, prefix + ".write_initial", false);
    result.write_final = optional_bool(
        entries, prefix + ".write_final", default_final);
    return result;
}

} // namespace

const char* partition_mode_name(PartitionMode mode)
{
    switch (mode) {
    case PartitionMode::ZonesOnly: return "zones_only";
    case PartitionMode::AutoSplit: return "auto_split";
    case PartitionMode::ForceSplit: return "force_split";
    }
    throw CaseConfigurationError("invalid partition mode");
}

const char* boundary_type_name(BoundaryType type)
{
    switch (type) {
    case BoundaryType::Farfield: return "farfield";
    case BoundaryType::Inflow: return "inflow";
    case BoundaryType::Outflow: return "outflow";
    case BoundaryType::SlipWall: return "slip_wall";
    case BoundaryType::NoSlipAdiabaticWall: return "no_slip_adiabatic_wall";
    case BoundaryType::NoSlipIsothermalWall: return "no_slip_isothermal_wall";
    case BoundaryType::Symmetry: return "symmetry";
    case BoundaryType::Periodic: return "periodic";
    case BoundaryType::Undefined: return "undefined";
    }
    throw CaseConfigurationError("invalid boundary type");
}

const char* run_mode_name(RunMode mode)
{
    switch (mode) {
    case RunMode::Steady: return "steady";
    case RunMode::Unsteady: return "unsteady";
    }
    throw CaseConfigurationError("invalid run mode");
}

const char* field_output_format_name(FieldOutputFormat format)
{
    switch (format) {
    case FieldOutputFormat::Cgns: return "cgns";
    case FieldOutputFormat::Tecplot: return "tecplot";
    case FieldOutputFormat::Both: return "both";
    }
    throw CaseConfigurationError("invalid field output format");
}

const char* series_output_format_name(SeriesOutputFormat format)
{
    switch (format) {
    case SeriesOutputFormat::Text: return "txt";
    case SeriesOutputFormat::Tecplot: return "tecplot";
    }
    throw CaseConfigurationError("invalid series output format");
}

void PartitionConfig::validate(AlgorithmProfileKind profile) const
{
    if (!std::isfinite(max_load_ratio) || max_load_ratio < 1.0) {
        throw CaseConfigurationError("partition max load ratio must be finite and >= 1");
    }
    const int strict_minimum
        = profile == AlgorithmProfileKind::PhengleiWcns ? 4 : 5;
    if (min_cells_per_active_direction < strict_minimum) {
        throw CaseConfigurationError(
            "partition minimum cells are below the selected profile limit");
    }
}

std::string PartitionConfig::summary() const
{
    std::ostringstream result;
    result << "partition(mode=" << partition_mode_name(mode)
           << ",allow_idle=" << (allow_idle_ranks ? "true" : "false")
           << ",max_load_ratio=" << std::setprecision(17) << max_load_ratio
           << ",min_cells=" << min_cells_per_active_direction << ')';
    return result.str();
}

Real InitialConditionConfig::parameter(
    const std::string& name,
    Real default_value) const
{
    const auto iterator = parameters.find(name);
    return iterator == parameters.end() ? default_value : iterator->second;
}

void InitialConditionConfig::validate(int dimension) const
{
    static const std::set<std::string> valid_types {
        "uniform", "quadrant_riemann", "sod_x", "isentropic_vortex",
        "couette", "manufactured_periodic",
    };
    if (valid_types.find(type) == valid_types.end()) {
        throw CaseConfigurationError("unknown initial condition type: " + type);
    }
    if (dimension != 2 && dimension != 3) {
        throw CaseConfigurationError("initial condition dimension must be 2 or 3");
    }
    for (const auto& [name, value] : parameters) {
        if (!std::isfinite(value)) {
            throw CaseConfigurationError(
                "initial condition parameter is not finite: " + name);
        }
    }
    if (parameter("rho", 1.0) <= 0.0
        || (parameters.find("temperature") != parameters.end()
            && parameter("temperature", 0.0) <= 0.0)
        || (parameters.find("pressure") != parameters.end()
            && parameter("pressure", 0.0) <= 0.0)) {
        throw CaseConfigurationError(
            "initial density, temperature and pressure must be positive");
    }
}

std::string InitialConditionConfig::summary() const
{
    std::vector<std::pair<std::string, Real>> sorted(
        parameters.begin(), parameters.end());
    std::sort(sorted.begin(), sorted.end());
    std::ostringstream result;
    result << "initial(type=" << type;
    for (const auto& [name, value] : sorted) {
        result << ',' << name << '=' << std::setprecision(17) << value;
    }
    result << ')';
    return result.str();
}

void SteadyConvergenceConfig::validate() const
{
    if (min_steps == 0 || check_interval_steps == 0
        || consecutive_checks == 0) {
        throw CaseConfigurationError(
            "steady step/check counts must be positive");
    }
    const std::array<std::pair<const char*, Real>, 5> values {{
        {"reference_floor", reference_floor},
        {"l2_absolute", l2_absolute},
        {"l2_relative", l2_relative},
        {"linf_absolute", linf_absolute},
        {"linf_relative", linf_relative},
    }};
    for (const auto& value : values) {
        if (!std::isfinite(value.second) || value.second <= 0.0) {
            throw CaseConfigurationError(
                std::string("steady ") + value.first
                + " must be finite and positive");
        }
    }
}

std::string SteadyConvergenceConfig::summary() const
{
    std::ostringstream result;
    result << "steady(min_steps=" << min_steps
           << ",check_interval=" << check_interval_steps
           << ",consecutive=" << consecutive_checks
           << ",reference_floor=" << std::setprecision(17) << reference_floor
           << ",l2_abs=" << l2_absolute << ",l2_rel=" << l2_relative
           << ",linf_enabled=" << (linf_enabled ? "true" : "false")
           << ",linf_abs=" << linf_absolute
           << ",linf_rel=" << linf_relative << ')';
    return result.str();
}

void OutputScheduleConfig::validate() const
{
    if (!std::isfinite(every_time) || every_time < 0.0) {
        throw CaseConfigurationError(
            "output schedule every_time must be finite and non-negative");
    }
    Real previous = -1.0;
    for (const Real time : explicit_times) {
        if (!std::isfinite(time) || time < 0.0 || time <= previous) {
            throw CaseConfigurationError(
                "output explicit times must be finite, non-negative and strictly increasing");
        }
        previous = time;
    }
}

std::string OutputScheduleConfig::summary() const
{
    std::ostringstream result;
    result << "schedule(every_steps=" << every_steps
           << ",every_time=" << std::setprecision(17) << every_time
           << ",explicit_times=";
    for (std::size_t index = 0; index < explicit_times.size(); ++index) {
        if (index != 0) result << ':';
        result << explicit_times[index];
    }
    result << ",initial=" << (write_initial ? "true" : "false")
           << ",final=" << (write_final ? "true" : "false") << ')';
    return result.str();
}

void FieldOutputConfig::validate() const
{
    schedule.validate();
    if (enabled && quantities.empty()) {
        throw CaseConfigurationError(
            "enabled field output requires at least one quantity");
    }
}

std::string FieldOutputConfig::summary() const
{
    std::ostringstream result;
    result << "field(enabled=" << (enabled ? "true" : "false")
           << ",format=" << field_output_format_name(format)
           << ',' << schedule.summary() << ",quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    result << ')';
    return result.str();
}

void SeriesOutputConfig::validate(const char* label) const
{
    schedule.validate();
    if (enabled && quantities.empty() && std::string(label) == "statistics") {
        throw CaseConfigurationError(
            "enabled statistics output requires at least one quantity");
    }
}

std::string SeriesOutputConfig::summary(const char* label) const
{
    std::ostringstream result;
    result << label << "(enabled=" << (enabled ? "true" : "false")
           << ",format=" << series_output_format_name(format)
           << ',' << schedule.summary() << ",quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    result << ')';
    return result.str();
}

void CheckpointOutputConfig::validate() const
{
    schedule.validate();
}

std::string CheckpointOutputConfig::summary() const
{
    return std::string("checkpoint(enabled=")
        + (enabled ? "true," : "false,") + schedule.summary() + ')';
}

void OutputConfig::validate() const
{
    if (directory.empty()) {
        throw CaseConfigurationError("output directory must not be empty");
    }
    field.validate();
    history.validate("history");
    statistics.validate("statistics");
    checkpoint.validate();
}

std::string OutputConfig::summary() const
{
    std::ostringstream result;
    result << "output(directory=" << directory
           << ",allow_existing=" << (allow_existing ? "true" : "false")
           << ",dimensional=" << (dimensional ? "true" : "false")
           << ',' << field.summary() << ',' << history.summary("history")
           << ',' << statistics.summary("statistics") << ','
           << checkpoint.summary() << ')';
    return result.str();
}

void CaseRunConfig::validate() const
{
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw CaseConfigurationError("run CFL must be finite and positive");
    }
    if (max_steps == 0) {
        throw CaseConfigurationError("run max_steps must be positive");
    }
    if (!std::isfinite(end_time) || end_time < 0.0
        || !std::isfinite(max_wall_time) || max_wall_time < 0.0) {
        throw CaseConfigurationError(
            "run t_end and max_wall_time must be finite and non-negative");
    }
    if (mode == RunMode::Unsteady && end_time <= 0.0) {
        throw CaseConfigurationError(
            "unsteady run requires a positive t_end");
    }
    steady.validate();
}

std::string CaseRunConfig::summary() const
{
    std::ostringstream result;
    result << "run(mode=" << run_mode_name(mode)
           << ",viscous=" << (viscous ? "true" : "false")
           << ",cfl=" << std::setprecision(17) << cfl
           << ",max_steps=" << max_steps << ",t_end=" << end_time
           << ",max_wall_time=" << max_wall_time << ','
           << steady.summary() << ')';
    return result.str();
}

CaseConfig CaseConfig::from_text(const std::string& text)
{
    const auto entries = parse_entries(text);
    reject_unknown_keys(entries);

    CaseConfig result;
    const auto version = parse_integer(require(entries, "schema_version"), "schema_version");
    if (version != supported_schema_version) {
        throw CaseConfigurationError("unsupported configuration schema version");
    }
    result.schema_version = static_cast<int>(version);
    result.case_name = require(entries, "case.name");
    result.mesh_path = require(entries, "mesh.path");
    result.profile = ProfileFactory::from_string(
        require(entries, "algorithm.profile")).kind();
    result.reconstruction.scheme = require(entries, "algorithm.reconstruction");
    result.reconstruction.variables = parse_reconstruction_variables(
        require(entries, "algorithm.reconstruction_variables"));
    result.riemann.scheme = require(entries, "algorithm.riemann");

    result.gas.gamma = parse_real(require(entries, "gas.gamma"), "gas.gamma");
    if (const auto iterator = entries.find("gas.molar_mass"); iterator != entries.end()) {
        result.gas.molar_mass = parse_real(iterator->second, iterator->first);
    }
    if (const auto iterator = entries.find("gas.specific_gas_constant");
        iterator != entries.end()) {
        result.gas.specific_gas_constant = parse_real(iterator->second, iterator->first);
    }

    result.reference.velocity = parse_real(
        require(entries, "reference.velocity"), "reference.velocity");
    result.reference.density = parse_real(
        require(entries, "reference.density"), "reference.density");
    result.reference.temperature = parse_real(
        require(entries, "reference.temperature"), "reference.temperature");
    result.reference.length = parse_real(
        require(entries, "reference.length"), "reference.length");
    result.reference.viscosity = parse_real(
        require(entries, "reference.viscosity"), "reference.viscosity");

    result.partition.mode = parse_partition_mode(require(entries, "partition.mode"));
    result.partition.allow_idle_ranks = parse_bool(
        require(entries, "partition.allow_idle_ranks"),
        "partition.allow_idle_ranks");
    result.partition.max_load_ratio = parse_real(
        require(entries, "partition.max_load_ratio"),
        "partition.max_load_ratio");
    const auto min_cells = parse_integer(
        require(entries, "partition.min_cells_per_active_direction"),
        "partition.min_cells_per_active_direction");
    if (min_cells < std::numeric_limits<int>::min()
        || min_cells > std::numeric_limits<int>::max()) {
        throw CaseConfigurationError("partition minimum cells exceed int range");
    }
    result.partition.min_cells_per_active_direction = static_cast<int>(min_cells);

    result.initial.type = require(entries, "initial.type");
    constexpr std::string_view initial_prefix = "initial.";
    for (const auto& [key, value] : entries) {
        if (key.size() >= initial_prefix.size()
            && key.compare(
                0,
                initial_prefix.size(),
                initial_prefix.data(),
                initial_prefix.size()) == 0
            && key != "initial.type") {
            result.initial.parameters.emplace(
                key.substr(initial_prefix.size()), parse_real(value, key));
        }
    }

    result.default_boundary = parse_boundary_type(require(entries, "boundary.default"));
    for (const auto& [key, value] : entries) {
        if (!dynamic_boundary_key(key)) continue;
        constexpr std::size_t prefix_size = std::string_view("boundary.").size();
        constexpr std::size_t suffix_size = std::string_view(".type").size();
        const auto name = key.substr(
            prefix_size, key.size() - prefix_size - suffix_size);
        result.boundary_overrides.emplace(name, parse_boundary_type(value));
    }

    result.source_terms.enable_source_terms = parse_bool(
        require(entries, "source.enabled"), "source.enabled");
    if (const auto iterator = entries.find("source.models"); iterator != entries.end()) {
        for (const auto& name : split_list(iterator->second)) {
            result.source_terms.models.push_back(parse_source_model(name));
        }
    }
    const std::array<std::string, 5> uniform_keys {{
        "source.uniform.rho", "source.uniform.momentum_x",
        "source.uniform.momentum_y", "source.uniform.momentum_z",
        "source.uniform.energy",
    }};
    const std::array<std::string, 5> manufactured_keys {{
        "source.manufactured.rho", "source.manufactured.momentum_x",
        "source.manufactured.momentum_y", "source.manufactured.momentum_z",
        "source.manufactured.energy",
    }};
    for (std::size_t component = 0; component < 5; ++component) {
        result.source_terms.uniform_conservative[component]
            = optional_real(entries, uniform_keys[component], 0.0);
        result.source_terms.manufactured_amplitude[component]
            = optional_real(entries, manufactured_keys[component], 0.0);
    }
    result.source_terms.body_acceleration = {{
        optional_real(entries, "source.body.ax", 0.0),
        optional_real(entries, "source.body.ay", 0.0),
        optional_real(entries, "source.body.az", 0.0),
    }};

    result.run.mode = parse_run_mode(require(entries, "run.mode"));
    result.run.viscous = parse_bool(require(entries, "run.viscous"), "run.viscous");
    result.run.cfl = parse_real(require(entries, "run.cfl"), "run.cfl");
    const auto max_steps = parse_integer(
        require(entries, "run.max_steps"), "run.max_steps");
    if (max_steps <= 0
        || static_cast<unsigned long long>(max_steps)
            > std::numeric_limits<std::size_t>::max()) {
        throw CaseConfigurationError("run max_steps is outside size_t range");
    }
    result.run.max_steps = static_cast<std::size_t>(max_steps);
    result.run.end_time = optional_real(entries, "run.t_end", 0.0);
    result.run.max_wall_time = optional_real(entries, "run.max_wall_time", 0.0);
    result.run.steady.min_steps = optional_size(
        entries, "steady.min_steps", 1);
    result.run.steady.check_interval_steps = optional_size(
        entries, "steady.check_interval_steps", 1);
    result.run.steady.consecutive_checks = optional_size(
        entries, "steady.consecutive_checks", 1);
    result.run.steady.reference_floor = optional_real(
        entries, "steady.reference_floor", 1.0e-30);
    result.run.steady.l2_absolute = optional_real(
        entries, "steady.l2_absolute", 1.0e-12);
    result.run.steady.l2_relative = optional_real(
        entries, "steady.l2_relative", 1.0e-8);
    result.run.steady.linf_enabled = optional_bool(
        entries, "steady.linf_enabled", true);
    result.run.steady.linf_absolute = optional_real(
        entries, "steady.linf_absolute", 1.0e-11);
    result.run.steady.linf_relative = optional_real(
        entries, "steady.linf_relative", 1.0e-8);

    result.output.directory = require(entries, "output.directory");
    result.output.allow_existing = parse_bool(
        require(entries, "output.allow_existing"), "output.allow_existing");
    result.output.dimensional = parse_bool(
        require(entries, "output.dimensional"), "output.dimensional");

    result.output.field.enabled = parse_bool(
        require(entries, "output.field.enabled"), "output.field.enabled");
    if (const auto iterator = entries.find("output.field.format");
        iterator != entries.end()) {
        result.output.field.format = parse_field_output_format(iterator->second);
    }
    result.output.field.schedule = parse_schedule(entries, "output.field", true);
    result.output.field.quantities = optional_string_list(
        entries, "output.field.quantities");

    result.output.history.enabled = parse_bool(
        require(entries, "output.history.enabled"), "output.history.enabled");
    if (const auto iterator = entries.find("output.history.format");
        iterator != entries.end()) {
        result.output.history.format = parse_series_output_format(iterator->second);
    }
    result.output.history.schedule = parse_schedule(entries, "output.history", true);
    result.output.history.quantities = optional_string_list(
        entries, "output.history.quantities");

    result.output.statistics.enabled = parse_bool(
        require(entries, "output.statistics.enabled"),
        "output.statistics.enabled");
    if (const auto iterator = entries.find("output.statistics.format");
        iterator != entries.end()) {
        result.output.statistics.format = parse_series_output_format(iterator->second);
    }
    result.output.statistics.schedule = parse_schedule(
        entries, "output.statistics", true);
    result.output.statistics.quantities = optional_string_list(
        entries, "output.statistics.quantities");

    result.output.checkpoint.enabled = parse_bool(
        require(entries, "output.checkpoint.enabled"),
        "output.checkpoint.enabled");
    result.output.checkpoint.schedule = parse_schedule(
        entries, "output.checkpoint", true);
    if (const auto iterator = entries.find("restart.path");
        iterator != entries.end()) {
        result.restart_path = iterator->second;
    }
    result.digest_ = fnv1a(canonical_entries(entries));
    result.validate();
    return result;
}

CaseConfig CaseConfig::from_file(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw CaseConfigurationError("cannot open case configuration: " + path);
    }
    std::ostringstream text;
    text << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw CaseConfigurationError("failed to read case configuration: " + path);
    }
    return from_text(text.str());
}

void CaseConfig::validate() const
{
    if (schema_version != supported_schema_version) {
        throw CaseConfigurationError("unsupported configuration schema version");
    }
    if (case_name.empty() || mesh_path.empty()) {
        throw CaseConfigurationError("case name and mesh path must not be empty");
    }
    const auto gas_model = make_gas_model();
    static_cast<void>(make_reference_scales(gas_model));
    static_cast<void>(make_profile());
    partition.validate(profile);
    initial.validate();
    run.validate();
    output.validate();
    if (run.max_wall_time > 0.0 && !output.checkpoint.enabled) {
        throw CaseConfigurationError(
            "positive max_wall_time requires checkpoint output");
    }
    if (default_boundary == BoundaryType::Undefined) {
        throw CaseConfigurationError("default boundary type must be defined");
    }
    for (const auto& [name, type] : boundary_overrides) {
        if (name.empty() || type == BoundaryType::Undefined) {
            throw CaseConfigurationError("boundary override is invalid");
        }
    }
    auto inviscid = make_inviscid_config();
    inviscid.validate();
}

GasModel CaseConfig::make_gas_model() const
{
    return GasModel::from_input(gas);
}

ReferenceScales CaseConfig::make_reference_scales(
    const GasModel& gas_model) const
{
    return ReferenceScales::derive(reference, gas_model);
}

AlgorithmProfile CaseConfig::make_profile() const
{
    return ProfileFactory::create(profile);
}

InviscidWcnsConfig CaseConfig::make_inviscid_config() const
{
    InviscidWcnsConfig result;
    result.reconstruction = reconstruction;
    result.riemann = riemann;
    result.source_terms = source_terms;
    return result;
}

std::string CaseConfig::summary() const
{
    const auto gas_model = make_gas_model();
    const auto reference_scales = make_reference_scales(gas_model);
    std::vector<std::pair<std::string, BoundaryType>> boundaries(
        boundary_overrides.begin(), boundary_overrides.end());
    std::sort(boundaries.begin(), boundaries.end());
    std::ostringstream result;
    result << "case(schema=" << schema_version << ",name=" << case_name
           << ",mesh=" << mesh_path << ",profile=" << make_profile().name()
           << "," << reconstruction.summary() << ',' << riemann.summary()
           << ',' << gas_model.summary() << ',' << reference_scales.summary()
           << ',' << partition.summary() << ',' << initial.summary()
           << ",boundary.default=" << boundary_type_name(default_boundary);
    for (const auto& [name, type] : boundaries) {
        result << ",boundary." << name << '=' << boundary_type_name(type);
    }
    result << ',' << source_terms.summary() << ',' << run.summary()
           << ',' << output.summary()
           << ",restart.path=" << (restart_path.empty() ? "<none>" : restart_path)
           << ",digest=0x" << std::hex << digest_ << ')';
    return result.str();
}

std::string CaseConfig::restart_signature() const
{
    std::vector<std::pair<std::string, BoundaryType>> boundaries(
        boundary_overrides.begin(), boundary_overrides.end());
    std::sort(boundaries.begin(), boundaries.end());
    std::ostringstream result;
    result << "schema=" << schema_version << ";profile="
           << make_profile().restart_signature() << ";reconstruction="
           << reconstruction.restart_signature() << ";riemann="
           << riemann.restart_signature() << ";gas="
           << make_gas_model().restart_signature() << ";reference="
           << make_reference_scales(make_gas_model()).restart_signature()
           << ";boundary.default=" << boundary_type_name(default_boundary);
    for (const auto& boundary : boundaries) {
        result << ";boundary." << boundary.first << '='
               << boundary_type_name(boundary.second);
    }
    result << ";source=" << source_terms.restart_signature()
           << ";viscous=" << (run.viscous ? "true" : "false");
    return result.str();
}

} // namespace wcns
