#pragma once

#include <wcns/core/field.hpp>
#include <wcns/mesh/topology.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <array>
#include <cstddef>
#include <string>

namespace wcns {

struct WcnsParameters {
    Real epsilon = 1.0e-6;
    int nonlinear_power = 2;

    void validate() const;
};

enum class ReconstructionKind {
    Linear5,
    WcnsJs,
};

enum class ReconstructionVariables {
    Conservative,
    Primitive,
};

struct ReconstructionConfig {
    ReconstructionKind kind = ReconstructionKind::WcnsJs;
    ReconstructionVariables variables = ReconstructionVariables::Primitive;
    WcnsParameters nonlinear {};
    ReconstructionScaling scaling {};
    NumericalFloors floors {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct ReconstructionDiagnostics {
    std::size_t nonlinear_faces = 0;
    std::size_t linear_faces = 0;
    std::size_t linear_fallbacks = 0;
    std::size_t first_order_fallbacks = 0;
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

[[nodiscard]] ScalarFaceStates linear5_reconstruct(
    const std::array<Real, 6>& stencil);

[[nodiscard]] ScalarFaceStates wcns5_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters = {});

[[nodiscard]] EulerFaceStates reconstruct_euler_face(
    const Field<Real>& primitive,
    Axis axis,
    Index3 face,
    const WcnsParameters& parameters = {});

[[nodiscard]] EulerFaceStates reconstruct_thermodynamic_face(
    const Field<Real>& conservative,
    const Field<Real>& pressure_primitive_field,
    Axis axis,
    Index3 face,
    const ReconstructionConfig& config,
    const GasModel& gas,
    const ReferenceScales& reference,
    ReconstructionDiagnostics& diagnostics,
    int dimension);

} // namespace wcns
