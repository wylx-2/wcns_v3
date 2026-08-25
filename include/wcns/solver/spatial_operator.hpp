#pragma once

#include <wcns/mesh/structured_block.hpp>
#include <wcns/solver/euler.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

namespace wcns {

struct SpatialParameters {
    IdealGas gas {};
    WcnsParameters wcns {};
    Real cfl = 0.4;

    void validate() const;
};

// Computes dU/dt in block.flow.residual. Primitive ghost cells must already
// have been populated by physical BC and connectivity exchange operations.
void compute_euler_residual(
    StructuredBlock& block,
    const IdealGas& gas = {},
    const WcnsParameters& parameters = {});

[[nodiscard]] Real stable_time_step(
    const StructuredBlock& block,
    Real cfl,
    const IdealGas& gas = {});

[[nodiscard]] Real residual_l2(const StructuredBlock& block);

} // namespace wcns
