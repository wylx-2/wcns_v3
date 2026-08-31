#include <wcns/mesh/high_order_metrics.hpp>

#include <wcns/mesh/linear_operators.hpp>
#include <wcns/mesh/metrics.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace wcns {

struct MetricFieldBuilderAccess {
    static CellCoordinates& cell_coordinates(MetricField& metric)
    {
        return metric.cell_coordinates_;
    }

    static Array3D<Real>& jacobian(MetricField& metric)
    {
        return metric.jacobian_;
    }

    static FaceAreaVectors& i_faces(MetricField& metric) { return metric.i_faces_; }
    static FaceAreaVectors& j_faces(MetricField& metric) { return metric.j_faces_; }
    static FaceAreaVectors& k_faces(MetricField& metric) { return metric.k_faces_; }
};

namespace {

using ScalarField = Array3D<Real>;

struct VectorFields {
    explicit VectorFields(Extent3 extent)
        : x(extent)
        , y(extent)
        , z(extent)
    {
    }

    ScalarField x;
    ScalarField y;
    ScalarField z;
};

Extent3 face_extent(Extent3 cells, int axis)
{
    ++cells[static_cast<std::size_t>(axis)];
    return cells;
}

Real finite_relative_difference(Real value, Real reference)
{
    return std::abs(value - reference) / reference;
}

ScalarField multiply(const ScalarField& lhs, const ScalarField& rhs)
{
    if (lhs.interior_extent() != rhs.interior_extent()) {
        throw GeometryError("cannot multiply geometry fields with different extents");
    }
    const auto extent = lhs.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k) = lhs(i, j, k) * rhs(i, j, k);
            }
        }
    }
    return result;
}

ScalarField combine(
    const ScalarField& a,
    const ScalarField& b,
    const ScalarField& c,
    const ScalarField& d,
    Real scale)
{
    const auto extent = a.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k)
                    = scale * (a(i, j, k) + b(i, j, k) - c(i, j, k) - d(i, j, k));
            }
        }
    }
    return result;
}

template<class LineOperation>
ScalarField apply_lines(
    const ScalarField& input,
    int axis,
    Extent3 output_extent,
    LineOperation operation)
{
    const auto input_extent = input.interior_extent();
    ScalarField output(output_extent);
    const int line_count = input_extent[static_cast<std::size_t>(axis)];

    const int outer_k = axis == 2 ? 1 : input_extent.nk;
    const int outer_j = axis == 1 ? 1 : input_extent.nj;
    const int outer_i = axis == 0 ? 1 : input_extent.ni;
    for (int k0 = 0; k0 < outer_k; ++k0) {
        for (int j0 = 0; j0 < outer_j; ++j0) {
            for (int i0 = 0; i0 < outer_i; ++i0) {
                std::vector<Real> line(static_cast<std::size_t>(line_count));
                for (int coordinate = 0; coordinate < line_count; ++coordinate) {
                    Index3 index {i0, j0, k0};
                    index[static_cast<std::size_t>(axis)] = coordinate;
                    line[static_cast<std::size_t>(coordinate)]
                        = input(index.i, index.j, index.k);
                }
                const auto transformed = operation(line);
                const int output_count
                    = output_extent[static_cast<std::size_t>(axis)];
                if (static_cast<int>(transformed.size()) != output_count) {
                    throw GeometryError("line operator returned an unexpected extent");
                }
                for (int coordinate = 0; coordinate < output_count; ++coordinate) {
                    Index3 index {i0, j0, k0};
                    index[static_cast<std::size_t>(axis)] = coordinate;
                    output(index.i, index.j, index.k)
                        = transformed[static_cast<std::size_t>(coordinate)];
                }
            }
        }
    }
    return output;
}

ScalarField refined_derivative(const ScalarField& field, int axis)
{
    return apply_lines(
        field,
        axis,
        field.interior_extent(),
        [](const std::vector<Real>& line) { return grid_delta(line); });
}

ScalarField product_derivative(
    const ScalarField& first,
    const ScalarField& second_derivative,
    int outer_axis)
{
    return refined_derivative(multiply(first, second_derivative), outer_axis);
}

ScalarField symmetric_component(
    const ScalarField& first,
    const ScalarField& second,
    int first_tangent,
    int second_tangent)
{
    const auto d_first_first = refined_derivative(first, first_tangent);
    const auto d_first_second = refined_derivative(first, second_tangent);
    const auto d_second_first = refined_derivative(second, first_tangent);
    const auto d_second_second = refined_derivative(second, second_tangent);
    const auto a = product_derivative(first, d_second_first, second_tangent);
    const auto b = product_derivative(second, d_first_second, first_tangent);
    const auto c = product_derivative(first, d_second_second, first_tangent);
    const auto d = product_derivative(second, d_first_first, second_tangent);
    return combine(a, b, c, d, 0.5);
}

ScalarField dot_product(
    const VectorFields& coordinates,
    const VectorFields& vectors)
{
    const auto extent = coordinates.x.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k) = coordinates.x(i, j, k) * vectors.x(i, j, k)
                    + coordinates.y(i, j, k) * vectors.y(i, j, k)
                    + coordinates.z(i, j, k) * vectors.z(i, j, k);
            }
        }
    }
    return result;
}

VectorFields refined_coordinates(const StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    const Extent3 refined {
        2 * vertices.ni - 1,
        2 * vertices.nj - 1,
        block.cell_dimension() == 3 ? 2 * vertices.nk - 1 : 1,
    };
    VectorFields result(refined);
    for (int rk = 0; rk < refined.nk; ++rk) {
        for (int rj = 0; rj < refined.nj; ++rj) {
            for (int ri = 0; ri < refined.ni; ++ri) {
                const int i0 = ri / 2;
                const int j0 = rj / 2;
                const int k0 = rk / 2;
                const int i_count = ri % 2 == 0 ? 1 : 2;
                const int j_count = rj % 2 == 0 ? 1 : 2;
                const int k_count
                    = block.cell_dimension() == 2 || rk % 2 == 0 ? 1 : 2;
                Real x = 0.0;
                Real y = 0.0;
                Real z = 0.0;
                const Real weight = 1.0 / static_cast<Real>(i_count * j_count * k_count);
                for (int dk = 0; dk < k_count; ++dk) {
                    for (int dj = 0; dj < j_count; ++dj) {
                        for (int di = 0; di < i_count; ++di) {
                            x += weight * block.coordinates.x(i0 + di, j0 + dj, k0 + dk);
                            y += weight * block.coordinates.y(i0 + di, j0 + dj, k0 + dk);
                            z += weight * block.coordinates.z(i0 + di, j0 + dj, k0 + dk);
                        }
                    }
                }
                result.x(ri, rj, rk) = x;
                result.y(ri, rj, rk) = y;
                result.z(ri, rj, rk) = z;
            }
        }
    }
    return result;
}

void compute_refined_symmetric_metrics(
    const VectorFields& coordinates,
    int dimension,
    VectorFields& s_i,
    VectorFields& s_j,
    VectorFields& s_k,
    ScalarField& jacobian)
{
    if (dimension == 2) {
        s_i.x = refined_derivative(coordinates.y, 1);
        const auto dx_eta = refined_derivative(coordinates.x, 1);
        const auto dy_xi = refined_derivative(coordinates.y, 0);
        s_j.y = refined_derivative(coordinates.x, 0);
        const auto extent = coordinates.x.interior_extent();
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                s_i.y(i, j, 0) = -dx_eta(i, j, 0);
                s_i.z(i, j, 0) = 0.0;
                s_j.x(i, j, 0) = -dy_xi(i, j, 0);
                s_j.z(i, j, 0) = 0.0;
            }
        }
        const auto di = refined_derivative(dot_product(coordinates, s_i), 0);
        const auto dj = refined_derivative(dot_product(coordinates, s_j), 1);
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, 0) = 0.5 * (di(i, j, 0) + dj(i, j, 0));
            }
        }
        return;
    }

    s_i.x = symmetric_component(coordinates.z, coordinates.y, 1, 2);
    s_i.y = symmetric_component(coordinates.x, coordinates.z, 1, 2);
    s_i.z = symmetric_component(coordinates.y, coordinates.x, 1, 2);
    s_j.x = symmetric_component(coordinates.z, coordinates.y, 2, 0);
    s_j.y = symmetric_component(coordinates.x, coordinates.z, 2, 0);
    s_j.z = symmetric_component(coordinates.y, coordinates.x, 2, 0);
    s_k.x = symmetric_component(coordinates.z, coordinates.y, 0, 1);
    s_k.y = symmetric_component(coordinates.x, coordinates.z, 0, 1);
    s_k.z = symmetric_component(coordinates.y, coordinates.x, 0, 1);

    const auto di = refined_derivative(dot_product(coordinates, s_i), 0);
    const auto dj = refined_derivative(dot_product(coordinates, s_j), 1);
    const auto dk = refined_derivative(dot_product(coordinates, s_k), 2);
    const auto extent = coordinates.x.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, k)
                    = (di(i, j, k) + dj(i, j, k) + dk(i, j, k)) / 3.0;
            }
        }
    }
}

void copy_refined_to_metric(
    const VectorFields& coordinates,
    const VectorFields& s_i,
    const VectorFields& s_j,
    const VectorFields& s_k,
    const ScalarField& jacobian,
    MetricField& metric)
{
    const auto cells = metric.jacobian().interior_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Index3 refined {2 * i + 1, 2 * j + 1,
                    metric.dimension() == 3 ? 2 * k + 1 : 0};
                MetricFieldBuilderAccess::cell_coordinates(metric).x(i, j, k)
                    = coordinates.x(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::cell_coordinates(metric).y(i, j, k)
                    = coordinates.y(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::cell_coordinates(metric).z(i, j, k)
                    = coordinates.z(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::jacobian(metric)(i, j, k)
                    = jacobian(refined.i, refined.j, refined.k);
            }
        }
    }
    const auto copy_faces = [&](const VectorFields& source, FaceAreaVectors& target,
                                int normal_axis) {
        const auto extent = target.x.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    Index3 refined {2 * i + 1, 2 * j + 1,
                        metric.dimension() == 3 ? 2 * k + 1 : 0};
                    refined[static_cast<std::size_t>(normal_axis)]
                        = 2 * Index3 {i, j, k}[static_cast<std::size_t>(normal_axis)];
                    target.x(i, j, k) = source.x(refined.i, refined.j, refined.k);
                    target.y(i, j, k) = source.y(refined.i, refined.j, refined.k);
                    target.z(i, j, k) = source.z(refined.i, refined.j, refined.k);
                }
            }
        }
    };
    copy_faces(s_i, MetricFieldBuilderAccess::i_faces(metric), 0);
    copy_faces(s_j, MetricFieldBuilderAccess::j_faces(metric), 1);
    if (metric.dimension() == 3) {
        copy_faces(s_k, MetricFieldBuilderAccess::k_faces(metric), 2);
    }
}

MetricField build_phenglei_metric(
    const StructuredBlock& block,
    const AlgorithmProfile& profile)
{
    const auto vertices = block.vertex_extent();
    for (int axis = 0; axis < block.cell_dimension(); ++axis) {
        if (vertices[static_cast<std::size_t>(axis)] < 3) {
            throw GeometryError("phenglei_wcns metrics require at least three vertices per active direction");
        }
    }
    auto coordinates = refined_coordinates(block);
    const auto refined_extent = coordinates.x.interior_extent();
    VectorFields s_i(refined_extent);
    VectorFields s_j(refined_extent);
    VectorFields s_k(refined_extent);
    ScalarField jacobian(refined_extent);
    compute_refined_symmetric_metrics(
        coordinates, block.cell_dimension(), s_i, s_j, s_k, jacobian);
    MetricField metric(profile.kind(), block.cell_extent(), block.cell_dimension());
    copy_refined_to_metric(coordinates, s_i, s_j, s_k, jacobian, metric);
    return metric;
}

GeometryDiagnostics capture_reference_geometry(const StructuredBlock& block)
{
    const auto cells = block.cell_extent();
    GeometryDiagnostics diagnostics(cells);
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                diagnostics.reference_volume(i, j, k)
                    = block.cell_metrics.volume(i, j, k);
            }
        }
    }
    const auto copy_area = [](const FaceMetric& source, ScalarField& target) {
        const auto extent = target.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    target(i, j, k) = source.area(i, j, k);
                }
            }
        }
    };
    copy_area(block.face_metrics.i_faces, diagnostics.reference_i_area);
    copy_area(block.face_metrics.j_faces, diagnostics.reference_j_area);
    if (block.cell_dimension() == 3) {
        copy_area(block.face_metrics.k_faces, diagnostics.reference_k_area);
    }
    return diagnostics;
}

void validate_metric(
    MetricField& metric,
    GeometryDiagnostics& diagnostics,
    const MetricBuildOptions& options)
{
    const auto cells = metric.jacobian().interior_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Real reference = diagnostics.reference_volume(i, j, k);
                Real& value = MetricFieldBuilderAccess::jacobian(metric)(i, j, k);
                const bool finite_positive = std::isfinite(value)
                    && value > options.floors.jacobian_floor(reference);
                const Real difference = finite_positive
                    ? finite_relative_difference(value, reference)
                    : std::numeric_limits<Real>::infinity();
                diagnostics.maximum_jacobian_relative_difference = std::max(
                    diagnostics.maximum_jacobian_relative_difference, difference);
                if (!finite_positive
                    || difference > options.maximum_reference_relative_difference) {
                    if (metric.profile() == AlgorithmProfileKind::PhengleiWcns
                        && options.fallback == MetricFallback::PhengleiFiniteVolume) {
                        value = reference;
                        ++diagnostics.fallback_cell_count;
                    } else {
                        throw GeometryError(
                            "high-order Jacobian failed strict reference validation at cell ("
                            + std::to_string(i) + ',' + std::to_string(j) + ','
                            + std::to_string(k) + ')');
                    }
                }
            }
        }
    }

    const auto validate_faces = [&](const FaceAreaVectors& faces,
                                    const ScalarField& reference) {
        const auto extent = reference.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real area = faces.area(i, j, k);
                    if (!std::isfinite(area)
                        || area <= options.floors.face_area_floor(reference(i, j, k))) {
                        throw GeometryError("high-order metric has a degenerate real face");
                    }
                }
            }
        }
    };
    validate_faces(metric.i_faces(), diagnostics.reference_i_area);
    validate_faces(metric.j_faces(), diagnostics.reference_j_area);
    if (metric.dimension() == 3) {
        validate_faces(metric.k_faces(), diagnostics.reference_k_area);
    }
}

} // namespace

void MetricBuildOptions::validate(const AlgorithmProfile& profile) const
{
    floors.validate();
    if (!std::isfinite(maximum_reference_relative_difference)
        || maximum_reference_relative_difference < 0.0) {
        throw std::invalid_argument(
            "maximum metric/reference relative difference must be finite and non-negative");
    }
    if (fallback == MetricFallback::PhengleiFiniteVolume
        && profile.kind() != AlgorithmProfileKind::PhengleiWcns) {
        throw ProfileError("finite-volume metric fallback is only valid for phenglei_wcns");
    }
}

Real FaceAreaVectors::area(int i, int j, int k) const
{
    return std::sqrt(
        x(i, j, k) * x(i, j, k) + y(i, j, k) * y(i, j, k)
        + z(i, j, k) * z(i, j, k));
}

MetricField::MetricField(
    AlgorithmProfileKind profile,
    Extent3 cell_extent,
    int dimension)
    : profile_(profile)
    , dimension_(dimension)
    , cell_coordinates_(cell_extent)
    , jacobian_(cell_extent)
    , i_faces_(face_extent(cell_extent, 0))
    , j_faces_(face_extent(cell_extent, 1))
    , k_faces_(face_extent(cell_extent, 2))
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("metric field dimension must be 2 or 3");
    }
}

const FaceAreaVectors& MetricField::k_faces() const
{
    if (dimension_ != 3) {
        throw std::logic_error("a 2D metric field has no K-face family");
    }
    return k_faces_;
}

MetricInitializationResult initialize_metric_field(
    StructuredBlock& block,
    const AlgorithmProfile& profile,
    const MetricBuildOptions& options)
{
    options.validate(profile);
    compute_metrics(block);
    auto diagnostics = capture_reference_geometry(block);
    MetricField metric = [&] {
        if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
            return build_phenglei_metric(block, profile);
        }
        throw ProfileError("scmm6_wcns metric method is not available before Stage I3");
    }();
    validate_metric(metric, diagnostics, options);
    return {std::move(metric), std::move(diagnostics)};
}

} // namespace wcns
