#pragma once

#include <wcns/runtime/output_manager.hpp>
#include <wcns/runtime/quantity_registry.hpp>

#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

struct OriginalZoneField {
    BlockId source_zone = invalid_block_id;
    std::string name;
    int cell_dimension = 0;
    Extent3 cell_extent {};
    std::vector<Real> x;
    std::vector<Real> y;
    std::vector<Real> z;
    std::unordered_map<std::string, std::vector<Real>> quantities;
};

struct OriginalFieldSnapshot {
    std::vector<OriginalZoneField> zones;
};

[[nodiscard]] OriginalFieldSnapshot gather_original_zone_fields(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    const FieldQuantityRegistry& registry,
    const std::vector<std::string>& quantities,
    const QuantityContext& context,
    RankId root = 0);

class ProductionFieldWriter {
public:
    ProductionFieldWriter(
        const MpiRuntime& mpi,
        const CaseConfig& config,
        const StructuredPartitionPlan& partition,
        const LocalBlockSet& local_blocks,
        const BlockMetricMap& metrics,
        QuantityContext quantity_context,
        std::string mesh_path,
        FieldQuantityRegistry registry = FieldQuantityRegistry::create_builtin());

    [[nodiscard]] std::vector<std::string> write(
        const SimulationState& state) const;

private:
    void write_cgns(
        const OriginalFieldSnapshot& snapshot,
        const SimulationState& state,
        const std::string& path) const;
    void write_tecplot(
        const OriginalFieldSnapshot& snapshot,
        const SimulationState& state,
        const std::string& path) const;

    const MpiRuntime& mpi_;
    const CaseConfig& config_;
    const StructuredPartitionPlan& partition_;
    const LocalBlockSet& local_blocks_;
    const BlockMetricMap& metrics_;
    QuantityContext quantity_context_;
    std::string mesh_path_;
    FieldQuantityRegistry registry_;
};

} // namespace wcns
