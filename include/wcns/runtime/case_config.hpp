#pragma once

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/physics/source_terms.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/transport_model.hpp>

#include <cstddef>
#include <cstdint>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

class CaseConfigurationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class PartitionMode {
    ZonesOnly,
    AutoSplit,
    ForceSplit,
};

struct PartitionConfig {
    PartitionMode mode = PartitionMode::AutoSplit;
    bool allow_idle_ranks = false;
    Real max_load_ratio = 1.2;
    int min_cells_per_active_direction = 8;

    void validate(AlgorithmProfileKind profile) const;
    [[nodiscard]] std::string summary() const;
};

struct InitialConditionConfig {
    std::string type = "uniform";
    std::unordered_map<std::string, Real> parameters;

    void validate(int dimension = 3) const;
    [[nodiscard]] Real parameter(
        const std::string& name,
        Real default_value) const;
    [[nodiscard]] std::string summary() const;
};

struct BoundaryPhysicalDataConfig {
    std::array<std::optional<Real>, 3> wall_velocity;
    std::optional<Real> wall_temperature;
    std::optional<Real> rho;
    std::optional<Real> u;
    std::optional<Real> v;
    std::optional<Real> w;
    std::optional<Real> temperature;
    std::optional<Real> pressure;

    [[nodiscard]] bool has_wall_velocity() const noexcept;
    [[nodiscard]] bool has_target_state() const noexcept;
    void validate() const;
    [[nodiscard]] std::string summary() const;
};

enum class RunMode {
    Steady,
    Unsteady,
};

struct SteadyConvergenceConfig {
    std::size_t min_steps = 1;
    std::size_t check_interval_steps = 1;
    std::size_t consecutive_checks = 1;
    Real reference_floor = 1.0e-30;
    Real l2_absolute = 1.0e-12;
    Real l2_relative = 1.0e-8;
    bool linf_enabled = true;
    Real linf_absolute = 1.0e-11;
    Real linf_relative = 1.0e-8;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct OutputScheduleConfig {
    std::size_t every_steps = 0;
    Real every_time = 0.0;
    std::vector<Real> explicit_times;
    bool write_initial = false;
    bool write_final = true;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

enum class FieldOutputFormat {
    Cgns,
    Tecplot,
    Both,
};

enum class SeriesOutputFormat {
    Text,
    Tecplot,
};

struct FieldOutputConfig {
    bool enabled = false;
    FieldOutputFormat format = FieldOutputFormat::Cgns;
    OutputScheduleConfig schedule;
    std::vector<std::string> quantities;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct SeriesOutputConfig {
    bool enabled = false;
    SeriesOutputFormat format = SeriesOutputFormat::Text;
    OutputScheduleConfig schedule;
    std::vector<std::string> quantities;

    void validate(const char* label) const;
    [[nodiscard]] std::string summary(const char* label) const;
};

struct CheckpointOutputConfig {
    bool enabled = false;
    OutputScheduleConfig schedule;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct OutputConfig {
    std::string directory = "output";
    bool allow_existing = false;
    bool dimensional = false;
    FieldOutputConfig field;
    SeriesOutputConfig history;
    SeriesOutputConfig statistics;
    CheckpointOutputConfig checkpoint;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct CaseRunConfig {
    RunMode mode = RunMode::Steady;
    bool viscous = false;
    Real cfl = 0.2;
    std::size_t max_steps = 1;
    Real end_time = 0.0;
    Real max_wall_time = 0.0;
    SteadyConvergenceConfig steady;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct CaseConfig {
    static constexpr int supported_schema_version = 1;

    int schema_version = supported_schema_version;
    std::string case_name;
    std::string mesh_path;
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    ReconstructionConfig reconstruction {};
    RiemannConfig riemann {};
    GasModelInput gas;
    ReferenceInput reference;
    PartitionConfig partition;
    InitialConditionConfig initial;
    BoundaryType default_boundary = BoundaryType::Farfield;
    std::unordered_map<std::string, BoundaryType> boundary_overrides;
    std::unordered_map<std::string, BoundaryPhysicalDataConfig> boundary_data;
    SourceTermConfig source_terms;
    CaseRunConfig run;
    OutputConfig output;
    std::string restart_path;

    [[nodiscard]] static CaseConfig from_text(const std::string& text);
    [[nodiscard]] static CaseConfig from_file(const std::string& path);

    void validate() const;
    [[nodiscard]] GasModel make_gas_model() const;
    [[nodiscard]] ReferenceScales make_reference_scales(
        const GasModel& gas_model) const;
    [[nodiscard]] AlgorithmProfile make_profile() const;
    [[nodiscard]] InviscidWcnsConfig make_inviscid_config() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
    [[nodiscard]] std::uint64_t digest() const noexcept { return digest_; }

private:
    std::uint64_t digest_ = 0;
};

[[nodiscard]] const char* partition_mode_name(PartitionMode mode);
[[nodiscard]] const char* boundary_type_name(BoundaryType type);
[[nodiscard]] const char* run_mode_name(RunMode mode);
[[nodiscard]] const char* field_output_format_name(FieldOutputFormat format);
[[nodiscard]] const char* series_output_format_name(SeriesOutputFormat format);

} // namespace wcns
