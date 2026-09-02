#pragma once

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/physics/source_terms.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/transport_model.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

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

struct CaseRunConfig {
    bool viscous = false;
    Real cfl = 0.2;
    std::size_t max_steps = 1;

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
    SourceTermConfig source_terms;
    CaseRunConfig run;

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

} // namespace wcns
