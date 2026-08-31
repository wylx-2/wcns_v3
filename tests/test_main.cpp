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

        std::cout << "WCNS unit tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
