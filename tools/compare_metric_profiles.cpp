#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/high_order_metrics.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct Difference {
    double sum_absolute = 0.0;
    double sum_squared = 0.0;
    double maximum_absolute = 0.0;
    double maximum_scale = 0.0;
    std::size_t samples = 0;

    void add(double lhs, double rhs)
    {
        if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
            throw std::runtime_error("metric comparison encountered non-finite data");
        }
        const double delta = std::abs(lhs - rhs);
        sum_absolute += delta;
        sum_squared += delta * delta;
        maximum_absolute = std::max(maximum_absolute, delta);
        maximum_scale = std::max({maximum_scale, std::abs(lhs), std::abs(rhs)});
        ++samples;
    }

    void print(const char* quantity) const
    {
        if (samples == 0) throw std::runtime_error("empty metric comparison");
        const double count = static_cast<double>(samples);
        const double relative = maximum_absolute
            / std::max(maximum_scale, std::numeric_limits<double>::min());
        std::cout << "comparison=" << quantity
                  << " samples=" << samples
                  << " l1=" << sum_absolute / count
                  << " l2=" << std::sqrt(sum_squared / count)
                  << " linf=" << maximum_absolute
                  << " relative_linf=" << relative << '\n';
    }
};

struct ProfileStatistics {
    double jacobian_minimum = std::numeric_limits<double>::infinity();
    double jacobian_maximum = -std::numeric_limits<double>::infinity();
    double volume_sum = 0.0;
    double maximum_reference_relative_difference = 0.0;
    double gcl_closure_linf = 0.0;
    std::size_t fallback_cells = 0;
    std::size_t cells = 0;

    void add(
        const wcns::MetricField& metric,
        const wcns::GeometryDiagnostics& diagnostics)
    {
        const auto extent = metric.jacobian().interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const double jacobian = metric.jacobian()(i, j, k);
                    if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
                        throw std::runtime_error("metric Jacobian is not finite and positive");
                    }
                    jacobian_minimum = std::min(jacobian_minimum, jacobian);
                    jacobian_maximum = std::max(jacobian_maximum, jacobian);
                    volume_sum += jacobian;

                    const double closure_x
                        = metric.i_faces().x(i + 1, j, k) - metric.i_faces().x(i, j, k)
                        + metric.j_faces().x(i, j + 1, k) - metric.j_faces().x(i, j, k);
                    const double closure_y
                        = metric.i_faces().y(i + 1, j, k) - metric.i_faces().y(i, j, k)
                        + metric.j_faces().y(i, j + 1, k) - metric.j_faces().y(i, j, k);
                    const double closure_z
                        = metric.i_faces().z(i + 1, j, k) - metric.i_faces().z(i, j, k)
                        + metric.j_faces().z(i, j + 1, k) - metric.j_faces().z(i, j, k);
                    gcl_closure_linf = std::max(
                        gcl_closure_linf,
                        std::sqrt(closure_x * closure_x
                            + closure_y * closure_y + closure_z * closure_z));
                    ++cells;
                }
            }
        }
        maximum_reference_relative_difference = std::max(
            maximum_reference_relative_difference,
            diagnostics.maximum_jacobian_relative_difference);
        fallback_cells += diagnostics.fallback_cell_count;
    }

    void print(const char* profile) const
    {
        std::cout << "profile=" << profile
                  << " cells=" << cells
                  << " jacobian_min=" << jacobian_minimum
                  << " jacobian_max=" << jacobian_maximum
                  << " volume_sum=" << volume_sum
                  << " max_reference_relative_difference="
                  << maximum_reference_relative_difference
                  << " fallback_cells=" << fallback_cells
                  << " gcl_closure_linf=" << gcl_closure_linf << '\n';
    }
};

void compare_array(
    const wcns::Array3D<wcns::Real>& lhs,
    const wcns::Array3D<wcns::Real>& rhs,
    Difference& difference)
{
    if (lhs.interior_extent() != rhs.interior_extent()) {
        throw std::runtime_error("metric arrays have different extents");
    }
    const auto extent = lhs.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                difference.add(lhs(i, j, k), rhs(i, j, k));
            }
        }
    }
}

void compare_faces(
    const wcns::FaceAreaVectors& lhs,
    const wcns::FaceAreaVectors& rhs,
    Difference& x,
    Difference& y,
    Difference& z)
{
    compare_array(lhs.x, rhs.x, x);
    compare_array(lhs.y, rhs.y, y);
    compare_array(lhs.z, rhs.z, z);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: wcns_compare_metric_profiles <mesh.cgns>\n";
        return 1;
    }
    try {
        const std::string mesh_path = argv[1];
        const auto phenglei_profile
            = wcns::ProfileFactory::from_string("phenglei_wcns");
        const auto scmm_profile
            = wcns::ProfileFactory::from_string("scmm6_wcns");
        wcns::CgnsReader reader;
        const auto metadata = reader.read_metadata(mesh_path);
        if (metadata.zones.empty()) {
            throw std::runtime_error("metric comparison mesh has no zones");
        }

        ProfileStatistics phenglei_statistics;
        ProfileStatistics scmm_statistics;
        Difference coordinate_x;
        Difference coordinate_y;
        Difference coordinate_z;
        Difference jacobian;
        Difference i_face_x;
        Difference i_face_y;
        Difference i_face_z;
        Difference j_face_x;
        Difference j_face_y;
        Difference j_face_z;

        for (const auto& zone : metadata.zones) {
            auto phenglei_block = reader.read_block(mesh_path, zone, 0, 0);
            auto scmm_block = reader.read_block(mesh_path, zone, 0, 0);
            const auto phenglei
                = wcns::initialize_metric_field(phenglei_block, phenglei_profile);
            const auto scmm
                = wcns::initialize_metric_field(scmm_block, scmm_profile);

            phenglei_statistics.add(phenglei.metric, phenglei.diagnostics);
            scmm_statistics.add(scmm.metric, scmm.diagnostics);
            compare_array(
                phenglei.metric.cell_coordinates().x,
                scmm.metric.cell_coordinates().x, coordinate_x);
            compare_array(
                phenglei.metric.cell_coordinates().y,
                scmm.metric.cell_coordinates().y, coordinate_y);
            compare_array(
                phenglei.metric.cell_coordinates().z,
                scmm.metric.cell_coordinates().z, coordinate_z);
            compare_array(phenglei.metric.jacobian(), scmm.metric.jacobian(), jacobian);
            compare_faces(
                phenglei.metric.i_faces(), scmm.metric.i_faces(),
                i_face_x, i_face_y, i_face_z);
            compare_faces(
                phenglei.metric.j_faces(), scmm.metric.j_faces(),
                j_face_x, j_face_y, j_face_z);
        }

        std::cout << std::setprecision(17)
                  << "check=metric_profile_comparison zones=" << metadata.zones.size()
                  << " dimension=" << metadata.zones.front().cell_dimension << '\n';
        phenglei_statistics.print("phenglei_wcns");
        scmm_statistics.print("scmm6_wcns");
        coordinate_x.print("cell_coordinate_x");
        coordinate_y.print("cell_coordinate_y");
        coordinate_z.print("cell_coordinate_z");
        jacobian.print("jacobian");
        i_face_x.print("i_face_x");
        i_face_y.print("i_face_y");
        i_face_z.print("i_face_z");
        j_face_x.print("j_face_x");
        j_face_y.print("j_face_y");
        j_face_z.print("j_face_z");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "metric profile comparison failed: " << error.what() << '\n';
        return 1;
    }
}
