#pragma once

#include <wcns/core/topology_field.hpp>
#include <wcns/runtime/structured_partition.hpp>
#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/transport_model.hpp>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class QuantityScale {
    Dimensionless,
    Density,
    Velocity,
    Pressure,
    Temperature,
    Momentum,
    Energy,
    SpecificEnergy,
    Viscosity,
    LengthPower,
};

struct QuantityDescriptor {
    std::string name;
    TopologyLocation location = TopologyLocation::Cell;
    std::vector<std::string> dependencies;
    std::string nondimensional_unit = "1";
    std::string dimensional_unit = "1";
    QuantityScale scale = QuantityScale::Dimensionless;
    int length_power = 0;

    void validate() const;
};

struct QuantityContext {
    GasModel gas;
    ReferenceScales reference;
    NumericalFloors floors;
    TransportModel transport;
    bool dimensional = false;
};

struct QuantityField {
    QuantityDescriptor descriptor;
    Extent3 extent {};
    std::vector<Real> values;
};

class IFieldQuantity {
public:
    virtual ~IFieldQuantity() = default;
    [[nodiscard]] virtual const QuantityDescriptor& descriptor() const = 0;
    [[nodiscard]] virtual Real evaluate_cell(
        const StructuredBlock& block,
        const MetricField& metric,
        Index3 index,
        const QuantityContext& context) const = 0;
};

class FieldQuantityRegistry {
public:
    FieldQuantityRegistry() = default;
    [[nodiscard]] static FieldQuantityRegistry create_builtin();

    void register_quantity(std::shared_ptr<const IFieldQuantity> quantity);
    [[nodiscard]] bool contains(const std::string& name) const noexcept;
    [[nodiscard]] const QuantityDescriptor& descriptor(
        const std::string& name) const;
    [[nodiscard]] QuantityField evaluate(
        const std::string& name,
        const StructuredBlock& block,
        const MetricField& metric,
        const QuantityContext& context) const;
    void validate_selection(const std::vector<std::string>& names) const;

private:
    void validate_dependencies(const std::string& root) const;

    std::unordered_map<std::string, std::shared_ptr<const IFieldQuantity>> quantities_;
};

struct StatisticContext {
    const MpiRuntime& mpi;
    const LocalBlockSet& local_blocks;
    const BlockMetricMap& metrics;
    const StructuredPartitionPlan& partition;
    const GlobalConservationWeights& conservation_weights;
    AlgorithmProfile profile;
    QuantityContext quantities;
};

class IStatisticQuantity {
public:
    virtual ~IStatisticQuantity() = default;
    [[nodiscard]] virtual const QuantityDescriptor& descriptor() const = 0;
    [[nodiscard]] virtual Real evaluate(const StatisticContext& context) const = 0;
};

class StatisticRegistry {
public:
    StatisticRegistry() = default;
    [[nodiscard]] static StatisticRegistry create_builtin();

    void register_quantity(std::shared_ptr<const IStatisticQuantity> quantity);
    [[nodiscard]] bool contains(const std::string& name) const noexcept;
    [[nodiscard]] Real evaluate(
        const std::string& name,
        const StatisticContext& context) const;
    void validate_selection(const std::vector<std::string>& names) const;

private:
    std::unordered_map<std::string, std::shared_ptr<const IStatisticQuantity>> quantities_;
};

[[nodiscard]] Real quantity_scale_factor(
    const QuantityDescriptor& descriptor,
    const QuantityContext& context,
    int dimension);

} // namespace wcns
