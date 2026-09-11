#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/viscous_wcns_solver.hpp>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <utility>
#include <vector>

namespace {

bool allocation_probe_enabled = false;
std::size_t allocation_count = 0;
std::size_t allocation_bytes = 0;

void* allocate(std::size_t bytes)
{
    const auto actual = bytes == 0 ? std::size_t {1} : bytes;
    if (void* memory = std::malloc(actual)) {
        if (allocation_probe_enabled) {
            ++allocation_count;
            allocation_bytes += bytes;
        }
        return memory;
    }
    throw std::bad_alloc();
}

struct AllocationSample {
    std::size_t count = 0;
    std::size_t bytes = 0;
};

template <class Function> AllocationSample measure(Function&& function)
{
    allocation_count = 0;
    allocation_bytes = 0;
    allocation_probe_enabled = true;
    try {
        function();
    } catch (...) {
        allocation_probe_enabled = false;
        throw;
    }
    allocation_probe_enabled = false;
    return {allocation_count, allocation_bytes};
}

void add_boundary(wcns::StructuredBlock& block, const char* name, wcns::FaceLocation face)
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
    constexpr int vertices = 25;
    StructuredBlock block(0, "allocation-probe", 0, 2, 2, {vertices, vertices, 1}, 3);
    for (int j = 0; j < vertices; ++j) {
        for (int i = 0; i < vertices; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i) / (vertices - 1);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j) / (vertices - 1);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    add_boundary(block, "i-lower", {Axis::I, Side::Lower});
    add_boundary(block, "i-upper", {Axis::I, Side::Upper});
    add_boundary(block, "j-lower", {Axis::J, Side::Lower});
    add_boundary(block, "j-upper", {Axis::J, Side::Upper});
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(block));
    return StructuredMesh(std::move(blocks));
}

} // namespace

void* operator new(std::size_t bytes)
{
    return allocate(bytes);
}

void* operator new[](std::size_t bytes)
{
    return allocate(bytes);
}

void operator delete(void* memory) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory) noexcept
{
    std::free(memory);
}

void operator delete(void* memory, std::size_t) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory, std::size_t) noexcept
{
    std::free(memory);
}

int main(int argc, char** argv)
{
    try {
        using namespace wcns;
        MpiRuntime mpi(argc, argv);
        if (mpi.size() != 1) {
            throw std::runtime_error("stage T allocation probe requires one rank");
        }
        auto mesh = make_mesh();
        const auto distribution
            = BlockDistribution::balanced({{0, mesh.block(0).cell_extent().size()}}, 1);
        distribution.apply(mesh);
        const auto topology = DistributedTopology::build(mesh, distribution);
        std::vector<StructuredBlock> local_storage {mesh.block(0)};
        LocalBlockSet local(0, std::move(local_storage), distribution);

        GasModelInput gas_input;
        gas_input.specific_gas_constant = 287.0;
        const auto gas = GasModel::from_input(gas_input);
        const auto reference
            = ReferenceScales::derive({340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
        const NumericalFloors floors;
        const auto profile = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
        const TemperaturePrimitiveState freestream {{1.0, 0.2, 0.0, 0.0, 1.0}};
        const auto conservative = thermodynamic_conservative(freestream, gas, reference, floors, 2);

        BlockMetricMap metrics;
        BlockBoundaryDataMap boundary_data;
        for (auto& block : local.blocks()) {
            const auto cells = block.cell_extent();
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    store_state(block.flow.conservative, {i, j, 0}, conservative);
                }
            }
            metrics.emplace(block.id(), initialize_metric_field(block, profile).metric);
            BoundaryDataMap data;
            for (const auto& patch : block.boundaries) {
                BoundaryData value;
                value.target_state = freestream;
                data.emplace(patch.name, value);
            }
            boundary_data.emplace(block.id(), std::move(data));
        }

        ViscousWcnsConfig config;
        config.inviscid.reconstruction.scheme = "linear5";
        ViscousWcnsSolver solver(mpi,
                                 local,
                                 mesh,
                                 topology,
                                 1,
                                 metrics,
                                 boundary_data,
                                 profile,
                                 gas,
                                 reference,
                                 floors,
                                 config);
        const auto first = measure([&] { solver.compute_residuals(0.0, 1); });
        const auto second = measure([&] { solver.compute_residuals(0.0, 2); });
        if (first.count == 0 || second.count == 0) {
            throw std::runtime_error("stage T allocation probe observed no allocations");
        }
        constexpr std::size_t baseline_allocations = 946214;
        constexpr std::size_t maximum_candidate_allocations = 94621;
        if (first.count > maximum_candidate_allocations
            || second.count > maximum_candidate_allocations) {
            throw std::runtime_error(
                "stage T residual allocations did not decrease by at least 90 percent");
        }
        if (second.count > first.count) {
            throw std::runtime_error("stage T residual workspace expanded on its second use");
        }
        std::cout << "stage_t_allocation_probe"
                  << " baseline_allocations=" << baseline_allocations
                  << " first_allocations=" << first.count << " first_bytes=" << first.bytes
                  << " second_allocations=" << second.count << " second_bytes=" << second.bytes
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        allocation_probe_enabled = false;
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
