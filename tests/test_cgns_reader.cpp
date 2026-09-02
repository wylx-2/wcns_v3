#include "test_support.hpp"

#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/metrics.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>

#include <cstdlib>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

// 验收二维 CGNS 坐标、边界范围和参考度量读取。
void test_2d(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.bases.size() == 1);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    WCNS_REQUIRE(metadata.bases.front().name == "Base2D");
    const auto& zone = metadata.zones.front();
    WCNS_REQUIRE(zone.block_id == 0);
    WCNS_REQUIRE(zone.name == "Zone2D");
    WCNS_REQUIRE(zone.cell_dimension == 2);
    WCNS_REQUIRE(zone.physical_dimension == 2);
    WCNS_REQUIRE(zone.vertex_extent == (wcns::Extent3 {5, 4, 1}));
    WCNS_REQUIRE(zone.cell_extent == (wcns::Extent3 {4, 3, 1}));

    auto block = reader.read_block(path, zone, 0, 3);
    WCNS_REQUIRE(block.id() == 0);
    WCNS_REQUIRE(block.owner_rank() == 0);
    WCNS_REQUIRE(block.ghost_width() == 3);
    WCNS_REQUIRE(block.coordinates.x(0, 0, 0) == 0.0);
    WCNS_REQUIRE(block.coordinates.x(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.y(4, 3, 0) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(4, 3, 0) == 0.0);
    WCNS_REQUIRE(block.boundaries.size() == 4);

    const auto& imin = block.boundaries[0];
    WCNS_REQUIRE(imin.name == "imin");
    WCNS_REQUIRE(imin.type == wcns::BoundaryType::Inflow);
    WCNS_REQUIRE(imin.face == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Lower}));
    WCNS_REQUIRE(imin.vertex_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 3, 0}}));
    WCNS_REQUIRE(imin.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 2, 0}}));
    WCNS_REQUIRE(imin.boundary_face_range.untyped() == (wcns::IndexRange3 {{0, 0, 0}, {0, 2, 0}}));

    const auto& imax = block.boundaries[1];
    WCNS_REQUIRE(imax.type == wcns::BoundaryType::Outflow);
    WCNS_REQUIRE(imax.face == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(imax.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{3, 0, 0}, {3, 2, 0}}));
    WCNS_REQUIRE(imax.boundary_face_range.untyped() == (wcns::IndexRange3 {{4, 0, 0}, {4, 2, 0}}));

    const auto& jmax = block.boundaries[3];
    WCNS_REQUIRE(jmax.type == wcns::BoundaryType::SlipWall);
    WCNS_REQUIRE(jmax.face == (wcns::FaceLocation {wcns::Axis::J, wcns::Side::Upper}));
    WCNS_REQUIRE(jmax.vertex_range.untyped() == (wcns::IndexRange3 {{4, 3, 0}, {0, 3, 0}}));
    WCNS_REQUIRE(jmax.adjacent_cell_range.untyped() == (wcns::IndexRange3 {{3, 2, 0}, {0, 2, 0}}));
    WCNS_REQUIRE(jmax.boundary_face_range.untyped() == (wcns::IndexRange3 {{3, 3, 0}, {0, 3, 0}}));

    wcns::compute_metrics(block);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_x(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_y(0, 0, 0), 1.0 / 6.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.volume(0, 0, 0), 1.0 / 12.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.jacobian(3, 2, 0), 1.0 / 12.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.area(0, 0, 0), 1.0 / 3.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_x(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_y(0, 0, 0), 0.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.area(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.normal_y(0, 0, 0), 1.0, 1.0e-14);

    double total_area = 0.0;
    for (int j = 0; j < block.cell_extent().nj; ++j) {
        for (int i = 0; i < block.cell_extent().ni; ++i) {
            total_area += block.cell_metrics.volume(i, j, 0);
        }
    }
    WCNS_REQUIRE_NEAR(total_area, 1.0, 1.0e-13);
}

// 验收三维 CGNS 坐标、边界范围和参考度量读取。
void test_3d(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.bases.size() == 1);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    const auto& zone = metadata.zones.front();
    WCNS_REQUIRE(zone.name == "Zone3D");
    WCNS_REQUIRE(zone.cell_dimension == 3);
    WCNS_REQUIRE(zone.physical_dimension == 3);
    WCNS_REQUIRE(zone.vertex_extent == (wcns::Extent3 {4, 3, 3}));
    WCNS_REQUIRE(zone.cell_extent == (wcns::Extent3 {3, 2, 2}));

    auto block = reader.read_block(path, zone, 2, 2);
    WCNS_REQUIRE(block.owner_rank() == 2);
    WCNS_REQUIRE(block.coordinates.x(3, 2, 2) == 3.0);
    WCNS_REQUIRE(block.coordinates.y(3, 2, 2) == 1.0);
    WCNS_REQUIRE(block.coordinates.z(3, 2, 2) == 0.5);
    WCNS_REQUIRE(block.boundaries.size() == 6);
    WCNS_REQUIRE(block.boundaries[0].type == wcns::BoundaryType::Farfield);
    WCNS_REQUIRE(
        block.boundaries[3].adjacent_cell_range.untyped()
        == (wcns::IndexRange3 {{2, 1, 0}, {0, 1, 1}}));
    WCNS_REQUIRE(
        block.boundaries[5].face
        == (wcns::FaceLocation {wcns::Axis::K, wcns::Side::Upper}));

    wcns::compute_metrics(block);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_x(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_y(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.center_z(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.cell_metrics.volume(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.area(0, 0, 0), 0.125, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.i_faces.normal_x(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.area(0, 0, 0), 0.25, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.j_faces.normal_y(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.k_faces.area(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(block.face_metrics.k_faces.normal_z(0, 0, 0), 1.0, 1.0e-14);

    double total_volume = 0.0;
    for (int k = 0; k < block.cell_extent().nk; ++k) {
        for (int j = 0; j < block.cell_extent().nj; ++j) {
            for (int i = 0; i < block.cell_extent().ni; ++i) {
                total_volume += block.cell_metrics.volume(i, j, k);
            }
        }
    }
    WCNS_REQUIRE_NEAR(total_volume, 1.5, 1.0e-13);
}

// 验收越出顶点范围的 CGNS 物理边界被确定性拒绝。
void test_out_of_extent_boundary(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    WCNS_REQUIRE(metadata.zones.size() == 1);
    try {
        static_cast<void>(reader.read_block(path, metadata.zones.front(), 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find(
                "boundary vertex PointRange end J index 4 is outside [0, 3]")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("out-of-extent CGNS boundary was accepted");
}

// 验收二维多块连接的强类型范围、变换和共享面。
void test_multiblock_2d(const char* path)
{
    wcns::CgnsReader reader;
    auto mesh = reader.read_mesh(path, 2, 3);
    WCNS_REQUIRE(mesh.block_count() == 2);
    const auto& left = mesh.block(0);
    const auto& right = mesh.block(1);
    WCNS_REQUIRE(left.name() == "Left2D");
    WCNS_REQUIRE(right.name() == "Right2D");
    WCNS_REQUIRE(left.connectivities.size() == 1);
    WCNS_REQUIRE(right.connectivities.size() == 1);

    const auto& connection = left.connectivities.front();
    WCNS_REQUIRE(connection.receiver_block == 0);
    WCNS_REQUIRE(connection.donor_block == 1);
    WCNS_REQUIRE(connection.donor_rank == 2);
    WCNS_REQUIRE(
        connection.receiver_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(
        connection.donor_face
        == (wcns::FaceLocation {wcns::Axis::J, wcns::Side::Lower}));
    WCNS_REQUIRE(
        connection.receiver_vertex_range.untyped()
        == (wcns::IndexRange3 {{4, 0, 0}, {4, 3, 0}}));
    WCNS_REQUIRE(
        connection.donor_vertex_range.untyped()
        == (wcns::IndexRange3 {{3, 0, 0}, {0, 0, 0}}));
    WCNS_REQUIRE(connection.transform == (wcns::IndexTransform {{{2, -1, 3}}}));
    WCNS_REQUIRE(
        connection.shared_face_range.untyped()
        == (wcns::IndexRange3 {{4, 0, 0}, {4, 2, 0}}));
    WCNS_REQUIRE(connection.ghost_width == 3);
}

// 验收三维多块连接的轴置换和 donor 范围映射。
void test_multiblock_3d(const char* path)
{
    wcns::CgnsReader reader;
    auto mesh = reader.read_mesh(path, 0, 2);
    WCNS_REQUIRE(mesh.block_count() == 2);
    const auto& connection = mesh.block(0).connectivities.front();
    WCNS_REQUIRE(connection.receiver_block == 0);
    WCNS_REQUIRE(connection.donor_block == 1);
    WCNS_REQUIRE(
        connection.receiver_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Upper}));
    WCNS_REQUIRE(
        connection.donor_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Lower}));
    WCNS_REQUIRE(connection.transform == (wcns::IndexTransform {{{1, -2, 3}}}));
    WCNS_REQUIRE(
        connection.donor_vertex_range.untyped()
        == (wcns::IndexRange3 {{0, 2, 0}, {0, 0, 2}}));
    WCNS_REQUIRE(
        connection.donor_adjacent_cell_range.untyped()
        == (wcns::IndexRange3 {{0, 1, 0}, {0, 0, 1}}));
}

// 验收缺少互反记录的单边连接不能进入求解拓扑。
void test_one_sided_connectivity(const char* path)
{
    wcns::CgnsReader reader;
    try {
        static_cast<void>(reader.read_mesh(path, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find("has no reciprocal donor record")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("one-sided CGNS connectivity was accepted");
}

// 验收引用未知 donor zone 的 CGNS 连接被拒绝。
void test_unknown_donor(const char* path)
{
    wcns::CgnsReader reader;
    try {
        static_cast<void>(reader.read_mesh(path, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find("unknown donor zone: Missing2D")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("unknown CGNS donor zone was accepted");
}

// 验收越出 zone 范围的连接 PointRange 被拒绝。
void test_out_of_extent_connectivity(const char* path)
{
    wcns::CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    auto reduced_zone = metadata.zones.front();
    reduced_zone.vertex_extent.nj = 3;
    reduced_zone.cell_extent.nj = 2;
    try {
        static_cast<void>(reader.read_block(path, reduced_zone, 0, 3));
    } catch (const wcns::CgnsError& error) {
        WCNS_REQUIRE(
            std::string(error.what()).find(
                "connectivity receiver vertex range end J index 3 is outside [0, 2]")
            != std::string::npos);
        return;
    }
    throw std::runtime_error("out-of-extent CGNS connectivity was accepted");
}

// 验收单 CGNS zone 的坐标 hyperslab、物理边界切片和兄弟块互反连接。
void test_partitioned_zone(const char* path)
{
    wcns::CgnsReader reader;
    const std::vector<wcns::CgnsPartitionLeaf> leaves {
        {0, 0, {0, 0, 0}, {2, 3, 1}, 0},
        {1, 0, {2, 0, 0}, {4, 3, 1}, 1},
    };
    auto partitioned = reader.read_partitioned_mesh(path, leaves, 1, 3);
    WCNS_REQUIRE(partitioned.global_mesh.block_count() == 2);
    WCNS_REQUIRE(partitioned.local_blocks.size() == 1);
    const auto& local = partitioned.local_blocks.front();
    WCNS_REQUIRE(local.id() == 1);
    WCNS_REQUIRE(local.owner_rank() == 1);
    WCNS_REQUIRE(local.cell_extent() == (wcns::Extent3 {2, 3, 1}));
    WCNS_REQUIRE_NEAR(local.coordinates.x(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(local.coordinates.x(2, 3, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(local.coordinates.y(2, 3, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE(local.boundaries.size() == 3);
    WCNS_REQUIRE(local.connectivities.size() == 1);
    WCNS_REQUIRE(local.connectivities.front().donor_block == 0);
    WCNS_REQUIRE(local.connectivities.front().donor_rank == 0);
    WCNS_REQUIRE(
        local.connectivities.front().receiver_face
        == (wcns::FaceLocation {wcns::Axis::I, wcns::Side::Lower}));
    WCNS_REQUIRE(
        local.connectivities.front().receiver_vertex_range.untyped()
        == (wcns::IndexRange3 {{0, 0, 0}, {0, 3, 0}}));
    WCNS_REQUIRE(std::isnan(
        partitioned.global_mesh.block(0).coordinates.x(0, 0, 0)));
    partitioned.global_mesh.validate_connectivities(false);
}

// 验收轴置换原连接在 receiver/donor 两侧二次切片后仍成对且可建立通信计划。
void test_partitioned_multiblock(const char* path)
{
    using namespace wcns;
    CgnsReader reader;
    const std::vector<CgnsPartitionLeaf> leaves {
        {0, 0, {0, 0, 0}, {4, 1, 1}, 0},
        {1, 0, {0, 1, 0}, {4, 3, 1}, 0},
        {2, 1, {0, 0, 0}, {2, 4, 1}, 0},
        {3, 1, {2, 0, 0}, {3, 4, 1}, 0},
    };
    auto partitioned = reader.read_partitioned_mesh(path, leaves, 0, 1);
    WCNS_REQUIRE(partitioned.global_mesh.block_count() == 4);
    WCNS_REQUIRE(partitioned.local_blocks.size() == 4);
    partitioned.global_mesh.validate_connectivities(false);

    const auto& first_interface = partitioned.global_mesh.block(0).connectivities;
    WCNS_REQUIRE(first_interface.size() == 2);
    bool found_axis_swapped_original = false;
    for (const auto& connection : first_interface) {
        if (connection.donor_block == 3) {
            found_axis_swapped_original = true;
            WCNS_REQUIRE(
                connection.transform == (IndexTransform {{{2, -1, 3}}}));
        }
    }
    WCNS_REQUIRE(found_axis_swapped_original);

    const auto distribution = BlockDistribution::balanced(
        {{0, 4}, {1, 8}, {2, 8}, {3, 4}},
        1);
    const auto topology = DistributedTopology::build(
        partitioned.global_mesh,
        distribution,
        false);
    WCNS_REQUIRE(topology.exchanges().size() == 8);
}

// 验收 CGNS Periodic_t 平移由当前面映射到 donor 面，并在结构二次切分后保持互逆。
void test_periodic_translation_2d(const char* path)
{
    using namespace wcns;
    CgnsReader reader;
    const auto mesh = reader.read_mesh(path, 0, 3);
    WCNS_REQUIRE(mesh.block_count() == 2);
    WCNS_REQUIRE(mesh.block(0).connectivities.size() == 2);
    WCNS_REQUIRE(mesh.block(1).connectivities.size() == 2);
    const auto& forward = mesh.block(0).connectivities[1].periodic;
    const auto& reverse = mesh.block(1).connectivities[1].periodic;
    WCNS_REQUIRE_NEAR(forward.translation[0], 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(reverse.translation[0], -1.0, 1.0e-14);
    WCNS_REQUIRE(reverse == forward.inverse());

    const std::vector<CgnsPartitionLeaf> leaves {
        {0, 0, {0, 0, 0}, {4, 8, 1}, 0},
        {1, 0, {4, 0, 0}, {8, 8, 1}, 1},
        {2, 1, {0, 0, 0}, {4, 8, 1}, 0},
        {3, 1, {4, 0, 0}, {8, 8, 1}, 1},
    };
    auto partitioned = reader.read_partitioned_mesh(path, leaves, 0, 3);
    partitioned.global_mesh.validate_connectivities(false);
    WCNS_REQUIRE(partitioned.global_mesh.block_count() == 4);
    int periodic_connections = 0;
    for (const auto& block : partitioned.global_mesh.blocks()) {
        WCNS_REQUIRE(block.connectivities.size() == 2);
        for (const auto& connection : block.connectivities) {
            if (std::abs(connection.periodic.translation[0]) > 0.5) {
                ++periodic_connections;
            }
        }
    }
    WCNS_REQUIRE(periodic_connections == 2);
}

// 验收 CGNS x-y-z 欧拉角按 current->donor 方向构造旋转，并对反向记录取严格逆变换。
void test_periodic_rotation_3d(const char* path)
{
    using namespace wcns;
    CgnsReader reader;
    const auto mesh = reader.read_mesh(path, 0, 2);
    WCNS_REQUIRE(mesh.block_count() == 2);
    const auto& forward = mesh.block(0).connectivities.front().periodic;
    const auto& reverse = mesh.block(1).connectivities.front().periodic;
    const auto rotated = forward.apply_vector({{1.0, 0.0, 0.0}});
    WCNS_REQUIRE(std::abs(rotated[0]) < 1.0e-6);
    WCNS_REQUIRE_NEAR(rotated[1], 1.0, 1.0e-12);
    WCNS_REQUIRE_NEAR(rotated[2], 0.0, 1.0e-12);
    const auto recovered = reverse.apply_vector(rotated);
    WCNS_REQUIRE_NEAR(recovered[0], 1.0, 1.0e-12);
    WCNS_REQUIRE_NEAR(recovered[1], 0.0, 1.0e-12);
    WCNS_REQUIRE_NEAR(recovered[2], 0.0, 1.0e-12);
    const auto inverse = forward.inverse();
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            WCNS_REQUIRE_NEAR(
                reverse.rotation[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)],
                inverse.rotation[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)],
                1.0e-7);
        }
    }
}

} // namespace

// 运行全部 CGNS 读取验收子项并汇总进程退出状态。
int main(int argc, char** argv)
{
    if (argc != 10) {
        std::cerr
            << "usage: wcns_cgns_reader_tests <2d.cgns> <3d.cgns> "
               "<invalid-2d.cgns> <multi-2d.cgns> <multi-3d.cgns> "
               "<one-sided-2d.cgns> <unknown-donor-2d.cgns> "
               "<periodic-2d.cgns> <periodic-3d.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        test_2d(argv[1]);
        test_3d(argv[2]);
        test_out_of_extent_boundary(argv[3]);
        test_multiblock_2d(argv[4]);
        test_multiblock_3d(argv[5]);
        test_one_sided_connectivity(argv[6]);
        test_unknown_donor(argv[7]);
        test_out_of_extent_connectivity(argv[4]);
        test_partitioned_zone(argv[1]);
        test_partitioned_multiblock(argv[4]);
        test_periodic_translation_2d(argv[8]);
        test_periodic_rotation_3d(argv[9]);
        std::cout << "CGNS reader tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
