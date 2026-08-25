#pragma once

#include <wcns/mesh/structured_block.hpp>
#include <wcns/solver/euler.hpp>

namespace wcns {

void update_primitive_cell(
    StructuredBlock& block,
    Index3 index,
    const IdealGas& gas = {});

void update_primitive_interior(
    StructuredBlock& block,
    const IdealGas& gas = {});

// Populates every physical-boundary ghost layer. Connectivity and periodic
// ghosts are intentionally owned by the halo exchanger, not this routine.
void fill_physical_boundaries(
    StructuredBlock& block,
    const PrimitiveState& prescribed,
    const IdealGas& gas = {});

} // namespace wcns
