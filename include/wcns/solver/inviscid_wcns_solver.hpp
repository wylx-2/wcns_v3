#pragma once

#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/source_operator.hpp>

#include <unordered_map>

namespace wcns {

using BlockMetricMap = std::unordered_map<BlockId, MetricField>;
using BlockBoundaryDataMap = std::unordered_map<BlockId, BoundaryDataMap>;

struct InviscidWcnsConfig {
    ReconstructionConfig reconstruction {};
    RiemannConfig riemann {};
    InviscidBoundaryOptions boundary {};
    SourceTermConfig source_terms {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

class InviscidWcnsSolver {
public:
    InviscidWcnsSolver(
        const MpiRuntime& mpi,
        LocalBlockSet& local_blocks,
        const StructuredMesh& global_mesh,
        const DistributedTopology& topology,
        int distribution_rank_count,
        BlockMetricMap& metrics,
        const BlockBoundaryDataMap& boundary_data,
        AlgorithmProfile profile,
        GasModel gas,
        ReferenceScales reference,
        NumericalFloors floors,
        InviscidWcnsConfig config = {});

    void compute_residuals(Real stage_time);
    void advance(Real time_step, Real initial_time);

    [[nodiscard]] Real global_residual_l2() const;
    [[nodiscard]] const ReconstructionDiagnostics& reconstruction_diagnostics() const noexcept
    {
        return reconstruction_diagnostics_;
    }

private:
    const MpiRuntime& mpi_;
    LocalBlockSet& local_blocks_;
    const StructuredMesh& global_mesh_;
    const DistributedTopology& topology_;
    HaloExchanger state_exchanger_;
    BlockMetricMap& metrics_;
    const BlockBoundaryDataMap& boundary_data_;
    AlgorithmProfile profile_;
    GasModel gas_;
    ReferenceScales reference_;
    NumericalFloors floors_;
    InviscidWcnsConfig config_;
    SourceTermRegistry source_registry_;
    RiemannSolver riemann_ {};
    std::uint64_t version_ = 0;
    ReconstructionDiagnostics reconstruction_diagnostics_ {};
};

} // namespace wcns
