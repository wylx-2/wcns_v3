#include <wcns/solver/euler_solver.hpp>

#include <wcns/solver/boundary_conditions.hpp>
#include <wcns/solver/time_integrator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace wcns {

EulerSolver::EulerSolver(
    const MpiRuntime& mpi,
    LocalBlockSet& local_blocks,
    const DistributedTopology& topology,
    int distribution_rank_count,
    PrimitiveState prescribed_state,
    SpatialParameters parameters)
    : mpi_(mpi)
    , local_blocks_(local_blocks)
    , topology_(topology)
    , exchanger_(mpi, topology, distribution_rank_count)
    , prescribed_state_(prescribed_state)
    , parameters_(parameters)
{
    parameters_.validate();
    static_cast<void>(to_conservative(prescribed_state_, parameters_.gas));
    if (local_blocks_.rank() != mpi_.rank()) {
        throw std::invalid_argument("local block set rank differs from MPI rank");
    }
}

void EulerSolver::compute_residuals()
{
    BlockFieldRegistry conservative_fields(euler_components);
    for (auto& block : local_blocks_.blocks()) {
        update_primitive_interior(block, parameters_.gas);
        conservative_fields.add(block.id(), block.flow.conservative);
    }
    exchanger_.exchange(conservative_fields);

    for (const auto& exchange : topology_.exchanges()) {
        if (exchange.receiver_rank != mpi_.rank()) {
            continue;
        }
        auto& receiver = local_blocks_.block(exchange.halo.receiver_block);
        for (const auto& pair : exchange.halo.cell_pairs) {
            update_primitive_cell(receiver, pair.receiver_ghost, parameters_.gas);
        }
    }
    for (auto& block : local_blocks_.blocks()) {
        fill_physical_boundaries(block, prescribed_state_, parameters_.gas);
        compute_euler_residual(block, parameters_);
    }
}

Real EulerSolver::global_time_step() const
{
    Real local = std::numeric_limits<Real>::infinity();
    for (const auto& block : local_blocks_.blocks()) {
        local = std::min(local, stable_time_step(block, parameters_.cfl, parameters_.gas));
    }
    const Real global = mpi_.min(local);
    if (!std::isfinite(global) || global <= 0.0) {
        throw PhysicsError("global CFL time step is not positive and finite");
    }
    return global;
}

Real EulerSolver::global_residual_l2() const
{
    Real local_sum = 0.0;
    Real local_count = 0.0;
    for (const auto& block : local_blocks_.blocks()) {
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        const Real value = block.flow.residual(i, j, k, component);
                        local_sum += value * value;
                        local_count += 1.0;
                    }
                }
            }
        }
    }
    const Real count = mpi_.sum(local_count);
    if (count <= 0.0) {
        throw PhysicsError("global residual norm has no cells");
    }
    return std::sqrt(mpi_.sum(local_sum) / count);
}

void EulerSolver::advance(Real time_step)
{
    std::vector<StructuredBlock*> blocks;
    blocks.reserve(local_blocks_.blocks().size());
    for (auto& block : local_blocks_.blocks()) {
        blocks.push_back(&block);
    }
    advance_ssprk3(blocks, time_step, [this] { compute_residuals(); });
}

Real EulerSolver::advance_cfl()
{
    compute_residuals();
    const Real time_step = global_time_step();
    advance(time_step);
    return time_step;
}

} // namespace wcns
