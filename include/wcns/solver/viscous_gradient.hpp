#pragma once

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/solver/viscous_flux.hpp>

#include <cstdint>

namespace wcns {

inline constexpr int gradient_operand_components
    = viscous_primitive_components * 3;

class GradientOperandFaceField {
public:
    GradientOperandFaceField(
        Extent3 cell_extent,
        int dimension,
        AlgorithmProfileKind profile,
        std::uint64_t version);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] std::uint64_t version() const noexcept { return version_; }
    [[nodiscard]] int halo_layers() const noexcept { return halo_layers_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] Field<Real>& field(Axis axis);
    [[nodiscard]] const Field<Real>& field(Axis axis) const;

private:
    AlgorithmProfileKind profile_;
    std::uint64_t version_;
    int halo_layers_;
    int dimension_;
    Field<Real> i_;
    Field<Real> j_;
    Field<Real> k_;
};

class PrimitiveGradientField {
public:
    PrimitiveGradientField(
        Extent3 cell_extent,
        int dimension,
        AlgorithmProfileKind profile,
        std::uint64_t version);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] std::uint64_t version() const noexcept { return version_; }
    [[nodiscard]] int halo_layers() const noexcept { return halo_layers_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] Field<Real>& values() noexcept { return values_; }
    [[nodiscard]] const Field<Real>& values() const noexcept { return values_; }

    [[nodiscard]] Real& operator()(
        Index3 cell, ViscousPrimitive variable, int cartesian_direction);
    [[nodiscard]] Real operator()(
        Index3 cell, ViscousPrimitive variable, int cartesian_direction) const;

private:
    AlgorithmProfileKind profile_;
    std::uint64_t version_;
    int halo_layers_;
    int dimension_;
    Field<Real> values_;
};

[[nodiscard]] TemperaturePrimitiveState interpolate_temperature_face(
    const StructuredBlock& block,
    const AlgorithmProfile& profile,
    Axis axis,
    Index3 face);

[[nodiscard]] GradientOperandFaceField compute_gradient_face_operands(
    const StructuredBlock& block,
    const MetricField& metric,
    const AlgorithmProfile& profile,
    std::uint64_t version);

[[nodiscard]] PrimitiveGradientField compute_primitive_gradients(
    const StructuredBlock& block,
    const MetricField& metric,
    const GradientOperandFaceField& operands,
    const AlgorithmProfile& profile);

} // namespace wcns
