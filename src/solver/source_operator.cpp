#include <wcns/solver/source_operator.hpp>

#include <cmath>
#include <stdexcept>

namespace wcns {

void add_source_terms(
    StructuredBlock& block,
    const MetricField& metric,
    const SourceTermRegistry& registry,
    Real stage_time)
{
    if (registry.empty()) return;
    if (metric.dimension() != block.cell_dimension()) {
        throw std::invalid_argument("source metric and block dimensions differ");
    }
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto source = registry.evaluate(
                    load_conservative(block.flow.conservative, {i, j, k}),
                    {{metric.cell_coordinates().x(i, j, k),
                        metric.cell_coordinates().y(i, j, k),
                        metric.cell_coordinates().z(i, j, k)}},
                    stage_time, block.cell_dimension());
                for (int component = 0; component < euler_components; ++component) {
                    auto& residual = block.flow.residual(i, j, k, component);
                    residual += source[static_cast<std::size_t>(component)];
                    if (!std::isfinite(residual)) {
                        throw PhysicsError("source accumulation produced non-finite residual");
                    }
                }
            }
        }
    }
}

ConservativeState volume_weighted_source(
    const StructuredBlock& block,
    const MetricField& metric,
    const SourceTermRegistry& registry,
    Real stage_time)
{
    ConservativeState result {};
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto source = registry.evaluate(
                    load_conservative(block.flow.conservative, {i, j, k}),
                    {{metric.cell_coordinates().x(i, j, k),
                        metric.cell_coordinates().y(i, j, k),
                        metric.cell_coordinates().z(i, j, k)}},
                    stage_time, block.cell_dimension());
                const Real jacobian = metric.jacobian()(i, j, k);
                for (int component = 0; component < euler_components; ++component) {
                    result[static_cast<std::size_t>(component)]
                        += jacobian * source[static_cast<std::size_t>(component)];
                }
            }
        }
    }
    return result;
}

} // namespace wcns
