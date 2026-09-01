#pragma once

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/physics/source_terms.hpp>
#include <wcns/solver/euler.hpp>

namespace wcns {

void add_source_terms(
    StructuredBlock& block,
    const MetricField& metric,
    const SourceTermRegistry& registry,
    Real stage_time);

[[nodiscard]] ConservativeState volume_weighted_source(
    const StructuredBlock& block,
    const MetricField& metric,
    const SourceTermRegistry& registry,
    Real stage_time);

} // namespace wcns
