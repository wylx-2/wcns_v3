#pragma once

#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/spatial_operator.hpp>

namespace wcns {

class EulerSolver {
public:
    EulerSolver(
        const MpiRuntime& mpi,
        LocalBlockSet& local_blocks,
        const DistributedTopology& topology,
        int distribution_rank_count,
        PrimitiveState prescribed_state,
        SpatialParameters parameters = {});

    // Refreshes interior primitive variables, exchanges conservative halos,
    // converts received ghosts, applies physical BCs, and computes residuals.
    void compute_residuals();

    [[nodiscard]] Real global_time_step() const;
    [[nodiscard]] Real global_residual_l2() const;

    void advance(Real time_step);
    [[nodiscard]] Real advance_cfl();

private:
    const MpiRuntime& mpi_;
    LocalBlockSet& local_blocks_;
    const DistributedTopology& topology_;
    HaloExchanger exchanger_;
    PrimitiveState prescribed_state_ {};
    SpatialParameters parameters_ {};
};

} // namespace wcns
