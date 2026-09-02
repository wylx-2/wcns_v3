#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/viscous_wcns_solver.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <utility>
#include <vector>

namespace {

void add_farfield(
    wcns::StructuredBlock& block,
    const char* name,
    wcns::FaceLocation face)
{
    using namespace wcns;
    const auto vertices = block.vertex_extent();
    const auto cells = block.cell_extent();
    BoundaryPatch patch;
    patch.name = name;
    patch.type = BoundaryType::Farfield;
    patch.face = face;
    if (face.axis == Axis::I) {
        const int vi = face.side == Side::Lower ? 0 : vertices.ni - 1;
        const int ci = face.side == Side::Lower ? 0 : cells.ni - 1;
        const int fi = face.side == Side::Lower ? 0 : cells.ni;
        patch.vertex_range = {{vi, 0, 0}, {vi, vertices.nj - 1, 0}};
        patch.adjacent_cell_range = {{ci, 0, 0}, {ci, cells.nj - 1, 0}};
        patch.boundary_face_range = {{fi, 0, 0}, {fi, cells.nj - 1, 0}};
    } else {
        const int vj = face.side == Side::Lower ? 0 : vertices.nj - 1;
        const int cj = face.side == Side::Lower ? 0 : cells.nj - 1;
        const int fj = face.side == Side::Lower ? 0 : cells.nj;
        patch.vertex_range = {{0, vj, 0}, {vertices.ni - 1, vj, 0}};
        patch.adjacent_cell_range = {{0, cj, 0}, {cells.ni - 1, cj, 0}};
        patch.boundary_face_range = {{0, fj, 0}, {cells.ni - 1, fj, 0}};
    }
    block.boundaries.push_back(std::move(patch));
}

wcns::StructuredMesh make_mesh()
{
    using namespace wcns;
    constexpr int vertices = 9;
    StructuredBlock left(0, "left", 0, 2, 2, {vertices, vertices, 1}, 3);
    StructuredBlock right(1, "right", 0, 2, 2, {vertices, vertices, 1}, 3);
    for (int j = 0; j < vertices; ++j) {
        for (int i = 0; i < vertices; ++i) {
            left.coordinates.x(i, j, 0) = static_cast<Real>(i);
            left.coordinates.y(i, j, 0) = static_cast<Real>(j);
            left.coordinates.z(i, j, 0) = 0.0;
            right.coordinates.x(i, j, 0) = static_cast<Real>(i + vertices - 1);
            right.coordinates.y(i, j, 0) = static_cast<Real>(j);
            right.coordinates.z(i, j, 0) = 0.0;
        }
    }
    left.connectivities.push_back({
        "left-right", 0, 1, 0,
        {Axis::I, Side::Upper}, {Axis::I, Side::Lower},
        {{vertices - 1, 0, 0}, {vertices - 1, vertices - 1, 0}},
        {{0, 0, 0}, {0, vertices - 1, 0}},
        {{vertices - 2, 0, 0}, {vertices - 2, vertices - 2, 0}},
        {{0, 0, 0}, {0, vertices - 2, 0}},
        {{vertices - 1, 0, 0}, {vertices - 1, vertices - 2, 0}},
        {{{1, 2, 3}}}, 3});
    right.connectivities.push_back({
        "right-left", 1, 0, 0,
        {Axis::I, Side::Lower}, {Axis::I, Side::Upper},
        {{0, 0, 0}, {0, vertices - 1, 0}},
        {{vertices - 1, 0, 0}, {vertices - 1, vertices - 1, 0}},
        {{0, 0, 0}, {0, vertices - 2, 0}},
        {{vertices - 2, 0, 0}, {vertices - 2, vertices - 2, 0}},
        {{0, 0, 0}, {0, vertices - 2, 0}},
        {{{1, 2, 3}}}, 3});
    add_farfield(left, "left-i-lower", {Axis::I, Side::Lower});
    add_farfield(left, "left-j-lower", {Axis::J, Side::Lower});
    add_farfield(left, "left-j-upper", {Axis::J, Side::Upper});
    add_farfield(right, "right-i-upper", {Axis::I, Side::Upper});
    add_farfield(right, "right-j-lower", {Axis::J, Side::Lower});
    add_farfield(right, "right-j-upper", {Axis::J, Side::Upper});
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

wcns::GasModel make_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

void run_profile(
    const wcns::MpiRuntime& mpi,
    wcns::AlgorithmProfileKind kind)
{
    using namespace wcns;
    auto mesh = make_mesh();
    std::vector<BlockLoad> loads;
    for (const auto& block : mesh.blocks()) {
        loads.push_back({block.id(), block.cell_extent().size()});
    }
    const auto distribution = BlockDistribution::balanced(
        std::move(loads), mpi.size());
    distribution.apply(mesh);
    const auto topology = DistributedTopology::build(mesh, distribution);
    std::vector<StructuredBlock> local_storage;
    for (const auto& block : mesh.blocks()) {
        if (block.owner_rank() == mpi.rank()) local_storage.push_back(block);
    }
    LocalBlockSet local(mpi.rank(), std::move(local_storage), distribution);
    const auto gas = make_gas();
    const auto reference = ReferenceScales::derive(
        {340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
    const NumericalFloors floors;
    const TemperaturePrimitiveState freestream {{1.0, 0.2, -0.1, 0.0, 1.0}};
    const auto conservative = thermodynamic_conservative(
        freestream, gas, reference, floors, 2);
    const auto profile = ProfileFactory::create(kind);
    BlockMetricMap metrics;
    BlockBoundaryDataMap boundary_data;
    for (auto& block : local.blocks()) {
        const auto cells = block.cell_extent();
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                store_state(block.flow.conservative, {i, j, 0}, conservative);
            }
        }
        metrics.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(
                initialize_metric_field(block, profile).metric));
        BoundaryDataMap data;
        for (const auto& patch : block.boundaries) {
            BoundaryData patch_data;
            patch_data.target_state = freestream;
            data.emplace(patch.name, patch_data);
        }
        boundary_data.emplace(block.id(), std::move(data));
    }
    ViscousWcnsConfig config;
    config.inviscid.reconstruction.kind = ReconstructionKind::Linear5;
    ViscousWcnsSolver solver(
        mpi, local, mesh, topology, distribution.rank_count(), metrics,
        boundary_data, profile, gas, reference, floors, config);
    const Real time_step = solver.global_time_step(0.2);
    WCNS_REQUIRE(std::isfinite(time_step));
    WCNS_REQUIRE(time_step > 0.0);
    WCNS_REQUIRE_NEAR(mpi.max(time_step), mpi.min(time_step), 0.0);
    solver.compute_residuals(0.0);
    WCNS_REQUIRE(solver.global_residual_l2() < 8.0e-11);
    solver.advance(std::min(time_step, 1.0e-3), 0.0);
    Real local_error = 0.0;
    for (const auto& block : local.blocks()) {
        const auto cells = block.cell_extent();
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto state = load_conservative(
                    block.flow.conservative, {i, j, 0});
                for (int component = 0; component < euler_components; ++component) {
                    local_error = std::max(local_error,
                        std::abs(state[static_cast<std::size_t>(component)]
                            - conservative[static_cast<std::size_t>(component)]));
                }
            }
        }
    }
    WCNS_REQUIRE(mpi.max(local_error) < 8.0e-11);
}

} // namespace

// 验收两套粘性 WCNS 驱动在单 rank/双 rank 多块网格上的自由流、halo 和一步状态一致性。
int main(int argc, char** argv)
{
    try {
        wcns::MpiRuntime mpi(argc, argv);
        run_profile(mpi, wcns::AlgorithmProfileKind::PhengleiWcns);
        run_profile(mpi, wcns::AlgorithmProfileKind::Scmm6Wcns);
        if (mpi.rank() == 0) {
            std::cout << "viscous WCNS solver tests passed with "
                      << mpi.size() << " ranks\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
