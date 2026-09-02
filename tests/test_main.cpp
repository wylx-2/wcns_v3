#include <wcns/core/types.hpp>

#include <cstdlib>
#include <iostream>
#include <type_traits>

void test_index();
void test_array3d();
void test_field();
void test_topology_field();
void test_structured_block();
void test_cgns_link();
void test_metrics();
void test_topology();
void test_distribution();
void test_euler();
void test_thermodynamics();
void test_wcns();
void test_spatial_operator();
void test_source_terms();
void test_algorithm_profile();
void test_geometry_line_operators();
void test_phenglei_high_order_metrics();
void test_phenglei_high_order_metrics_3d();
void test_scmm6_high_order_metrics();
void test_scmm6_high_order_metrics_3d();
void test_warped_metric_convergence();
void test_geometry_halo_plan();
void test_shared_metric_synchronization();
void test_periodic_shared_metric_synchronization();
void test_line_conservation_weights();
void test_global_conservation_weights();
void test_physical_ghost_state();
void test_inviscid_boundary_face_state();
void test_stage_j_scalar_reconstruction();
void test_reconstruction_positivity_fallback();
void test_rusanov_riemann_solver();
void test_wcns_inviscid_freestream();
void test_wcns_strong_wall_flux();
void test_face_flux_halo_plan();
void test_wcns_flux_divergence_convergence();
void test_stage_j_source_models();
void test_source_operator_balance();
void test_timed_ssprk3_sources();
void test_transport_model();
void test_viscous_cartesian_flux();
void test_viscous_linear_gradients();
void test_viscous_gradient_periodic_transform();
void test_wall_dirichlet_derivative();
void test_viscous_boundary_trace();

// 顺序运行不依赖独立输入文件的全部单元验收入口。
int main()
{
    try {
        static_assert(std::is_same_v<wcns::Real, double>);
        static_assert(std::is_signed_v<wcns::BlockId>);
        static_assert(std::is_signed_v<wcns::RankId>);

        if (wcns::invalid_block_id >= 0 || wcns::invalid_rank_id >= 0) {
            std::cerr << "invalid identifiers must be negative\n";
            return EXIT_FAILURE;
        }

        test_index();
        test_array3d();
        test_field();
        test_topology_field();
        test_structured_block();
        test_cgns_link();
        test_metrics();
        test_topology();
        test_distribution();
        test_euler();
        test_thermodynamics();
        test_wcns();
        test_spatial_operator();
        test_source_terms();
        test_algorithm_profile();
        test_geometry_line_operators();
        test_phenglei_high_order_metrics();
        test_phenglei_high_order_metrics_3d();
        test_scmm6_high_order_metrics();
        test_scmm6_high_order_metrics_3d();
        test_warped_metric_convergence();
        test_geometry_halo_plan();
        test_shared_metric_synchronization();
        test_periodic_shared_metric_synchronization();
        test_line_conservation_weights();
        test_global_conservation_weights();
        test_physical_ghost_state();
        test_inviscid_boundary_face_state();
        test_stage_j_scalar_reconstruction();
        test_reconstruction_positivity_fallback();
        test_rusanov_riemann_solver();
        test_wcns_inviscid_freestream();
        test_wcns_strong_wall_flux();
        test_face_flux_halo_plan();
        test_wcns_flux_divergence_convergence();
        test_stage_j_source_models();
        test_source_operator_balance();
        test_timed_ssprk3_sources();
        test_transport_model();
        test_viscous_cartesian_flux();
        test_viscous_linear_gradients();
        test_viscous_gradient_periodic_transform();
        test_wall_dirichlet_derivative();
        test_viscous_boundary_trace();

        std::cout << "WCNS unit tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
