#pragma once

#include <wcns/core/field.hpp>
#include <wcns/mesh/topology.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>
#include <wcns/solver/numerical_diagnostics.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace wcns {

struct WcnsParameters {
    Real epsilon = 1.0e-6;
    int nonlinear_power = 2;
    int weno_z_power = 2;
    Real mdcd_dispersion = 0.0463783;
    Real mdcd_dissipation = 0.01;
    Real mdcd_sensor_epsilon = 5.625e-5;
    Real mdcd_sensor_threshold = 0.4;
    Real mdcd_weight_bias = 20.0;
    Real mdcd_weight_epsilon = 1.0e-40;

    void validate() const;
};

enum class ReconstructionKind {
    Linear5,
    WenoJs,
    WcnsJs = WenoJs,
    WenoZ,
    MdcdLinear,
    MdcdHybrid,
};

enum class ReconstructionVariables {
    Conservative,
    Primitive,
    Characteristic,
};

enum class TraceSide {
    Left,
    Right,
};

class ScalarStencilView {
public:
    ScalarStencilView(const Real* values, std::size_t size)
        : values_(values), size_(size)
    {
    }

    template<std::size_t Size>
    ScalarStencilView(const std::array<Real, Size>& values)
        : ScalarStencilView(values.data(), values.size())
    {
    }

    ScalarStencilView(const std::vector<Real>& values)
        : ScalarStencilView(values.data(), values.size())
    {
    }

    [[nodiscard]] const Real* begin() const noexcept { return values_; }
    [[nodiscard]] const Real* end() const noexcept { return values_ + size_; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] const Real& operator[](std::size_t index) const
    {
        if (index >= size_) throw std::out_of_range("scalar stencil index is out of range");
        return values_[index];
    }

private:
    const Real* values_ = nullptr;
    std::size_t size_ = 0;
};

struct StencilRequirement {
    int first_offset = -2;
    int point_count = 6;
    int halo_width = 3;
};

struct ReconstructionContext {
    Real scale = 1.0;
    WcnsParameters parameters {};
};

class IReconstructionScheme {
public:
    virtual ~IReconstructionScheme() = default;
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;
    [[nodiscard]] virtual StencilRequirement stencil_requirement() const noexcept = 0;
    [[nodiscard]] virtual Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const = 0;
};

class ReconstructionRegistry {
public:
    using Factory = std::function<std::unique_ptr<IReconstructionScheme>()>;

    void register_scheme(std::string name, Factory factory);
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] std::unique_ptr<IReconstructionScheme> create(
        std::string_view name) const;
    [[nodiscard]] std::vector<std::string> names() const;

    [[nodiscard]] static ReconstructionRegistry with_builtins();

private:
    std::unordered_map<std::string, Factory> factories_;
};

[[nodiscard]] std::string_view reconstruction_name(ReconstructionKind kind);

struct ReconstructionConfig {
    std::string scheme = "weno_js";
    ReconstructionVariables variables = ReconstructionVariables::Primitive;
    WcnsParameters nonlinear {};
    ReconstructionScaling scaling {};
    NumericalFloors floors {};

    void validate() const;
    void validate(const ReconstructionRegistry& registry) const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

enum class ReconstructionFallbackReason {
    InvalidCharacteristicState,
    InvalidReconstructedState,
};

struct ReconstructionFallbackEvent {
    FaceDiagnosticLocation location {};
    std::string requested_scheme;
    std::string from_strategy;
    std::string to_strategy;
    ReconstructionFallbackReason reason
        = ReconstructionFallbackReason::InvalidReconstructedState;
};

struct ReconstructionDiagnostics {
    std::size_t nonlinear_faces = 0;
    std::size_t linear_faces = 0;
    std::size_t characteristic_faces = 0;
    std::size_t characteristic_fallbacks = 0;
    std::size_t primitive_fallbacks = 0;
    std::size_t linear_fallbacks = 0;
    std::size_t first_order_fallbacks = 0;
    std::vector<ReconstructionFallbackEvent> fallback_events;

    void record_fallback(
        FaceDiagnosticLocation location,
        std::string requested_scheme,
        std::string from_strategy,
        std::string to_strategy,
        ReconstructionFallbackReason reason);
};

struct ScalarFaceStates {
    Real left = 0.0;
    Real right = 0.0;
};

struct EulerFaceStates {
    PrimitiveState left {};
    PrimitiveState right {};
};

using CharacteristicMatrix = std::array<
    std::array<Real, euler_components>, euler_components>;

struct EulerCharacteristicBasis {
    Normal3 normal {};
    Normal3 tangent_1 {};
    Normal3 tangent_2 {};
    CharacteristicMatrix left {};
    CharacteristicMatrix right {};
};

[[nodiscard]] EulerCharacteristicBasis make_roe_characteristic_basis(
    const PressurePrimitiveState& left,
    const PressurePrimitiveState& right,
    Normal3 unit_normal,
    const GasModel& gas,
    const NumericalFloors& floors,
    int dimension);

[[nodiscard]] ConservativeState project_characteristic(
    const ConservativeState& conservative,
    const EulerCharacteristicBasis& basis);

[[nodiscard]] ConservativeState restore_characteristic(
    const ConservativeState& characteristic,
    const EulerCharacteristicBasis& basis);

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

[[nodiscard]] ScalarFaceStates weno_z_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters = {});

[[nodiscard]] ScalarFaceStates mdcd_linear_reconstruct(
    const std::array<Real, 6>& stencil,
    const WcnsParameters& parameters = {});

[[nodiscard]] ScalarFaceStates mdcd_hybrid_reconstruct_scaled(
    const std::array<Real, 6>& stencil,
    Real scale,
    const WcnsParameters& parameters = {});

[[nodiscard]] Real mdcd_six_point_smoothness(
    const std::array<Real, 6>& stencil,
    Real scale);

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
    int dimension,
    Normal3 unit_normal = {1.0, 0.0, 0.0},
    FaceDiagnosticLocation location = {});

} // namespace wcns
