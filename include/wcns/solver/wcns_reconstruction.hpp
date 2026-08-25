#pragma once

#include <wcns/core/field.hpp>
#include <wcns/mesh/topology.hpp>
#include <wcns/solver/euler.hpp>

#include <array>

namespace wcns {

struct WcnsParameters {
    Real epsilon = 1.0e-6;
    int nonlinear_power = 2;

    void validate() const;
};

struct ScalarFaceStates {
    Real left = 0.0;
    Real right = 0.0;
};

struct EulerFaceStates {
    PrimitiveState left {};
    PrimitiveState right {};
};

[[nodiscard]] Real wcns5_left_interpolation(
    const std::array<Real, 5>& stencil,
    const WcnsParameters& parameters = {});

[[nodiscard]] ScalarFaceStates wcns5_reconstruct(
    const std::array<Real, 6>& stencil,
    const WcnsParameters& parameters = {});

[[nodiscard]] EulerFaceStates reconstruct_euler_face(
    const Field<Real>& primitive,
    Axis axis,
    Index3 face,
    const WcnsParameters& parameters = {});

} // namespace wcns

