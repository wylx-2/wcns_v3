#pragma once

#include <wcns/runtime/case_config.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>

namespace wcns {

class FlowInitializationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class FlowInitializer {
public:
    static void initialize_block(
        StructuredBlock& block,
        const MetricField& metrics,
        const InitialConditionConfig& config,
        const GasModel& gas,
        const ReferenceScales& reference,
        const NumericalFloors& floors = {});

    static void initialize_local_blocks(
        LocalBlockSet& local_blocks,
        const BlockMetricMap& metrics,
        const InitialConditionConfig& config,
        const GasModel& gas,
        const ReferenceScales& reference,
        const NumericalFloors& floors = {});

    [[nodiscard]] static TemperaturePrimitiveState evaluate(
        const InitialConditionConfig& config,
        std::array<Real, 3> coordinates,
        const GasModel& gas,
        const ReferenceScales& reference,
        const NumericalFloors& floors,
        int dimension);
};

} // namespace wcns
