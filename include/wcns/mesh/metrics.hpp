#pragma once

#include <wcns/mesh/structured_block.hpp>

#include <stdexcept>

namespace wcns {

class GeometryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Computes interior cell centers, cell volumes/Jacobians, and positive
// computational-direction face normals and areas from node coordinates.
void compute_metrics(StructuredBlock& block);

} // namespace wcns

