#pragma once

#include <wcns/core/array3d.hpp>
#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/metrics.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/physics/thermodynamics.hpp>

#include <cstddef>

namespace wcns {

enum class MetricFallback {
    Strict,
    PhengleiFiniteVolume,
};

struct MetricBuildOptions {
    NumericalFloors floors {};
    MetricFallback fallback = MetricFallback::Strict;
    Real maximum_reference_relative_difference = 0.2;

    void validate(const AlgorithmProfile& profile) const;
};

struct CellCoordinates {
    explicit CellCoordinates(Extent3 extent)
        : x(extent)
        , y(extent)
        , z(extent)
    {
    }

    Array3D<Real> x;
    Array3D<Real> y;
    Array3D<Real> z;
};

struct FaceAreaVectors {
    explicit FaceAreaVectors(Extent3 extent)
        : x(extent)
        , y(extent)
        , z(extent)
    {
    }

    [[nodiscard]] Real area(int i, int j, int k) const;

    Array3D<Real> x;
    Array3D<Real> y;
    Array3D<Real> z;
};

class MetricField {
public:
    MetricField(AlgorithmProfileKind profile, Extent3 cell_extent, int dimension);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] const CellCoordinates& cell_coordinates() const noexcept
    {
        return cell_coordinates_;
    }
    [[nodiscard]] const Array3D<Real>& jacobian() const noexcept { return jacobian_; }
    [[nodiscard]] const FaceAreaVectors& i_faces() const noexcept { return i_faces_; }
    [[nodiscard]] const FaceAreaVectors& j_faces() const noexcept { return j_faces_; }
    [[nodiscard]] const FaceAreaVectors& k_faces() const;

private:
    friend struct MetricFieldBuilderAccess;
    friend class SharedMetricSynchronizer;

    AlgorithmProfileKind profile_;
    int dimension_;
    CellCoordinates cell_coordinates_;
    Array3D<Real> jacobian_;
    FaceAreaVectors i_faces_;
    FaceAreaVectors j_faces_;
    FaceAreaVectors k_faces_;
};

struct GeometryDiagnostics {
    explicit GeometryDiagnostics(Extent3 cell_extent)
        : reference_volume(cell_extent)
        , reference_i_area({cell_extent.ni + 1, cell_extent.nj, cell_extent.nk})
        , reference_j_area({cell_extent.ni, cell_extent.nj + 1, cell_extent.nk})
        , reference_k_area({cell_extent.ni, cell_extent.nj, cell_extent.nk + 1})
    {
    }

    Array3D<Real> reference_volume;
    Array3D<Real> reference_i_area;
    Array3D<Real> reference_j_area;
    Array3D<Real> reference_k_area;
    Real maximum_jacobian_relative_difference = 0.0;
    std::size_t fallback_cell_count = 0;
};

struct MetricInitializationResult {
    MetricField metric;
    GeometryDiagnostics diagnostics;
};

[[nodiscard]] MetricInitializationResult initialize_metric_field(
    StructuredBlock& block,
    const AlgorithmProfile& profile,
    const MetricBuildOptions& options = {});

} // namespace wcns
