#pragma once

#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/viscous_operator.hpp>

namespace wcns {

struct ViscousStabilityCoefficients {
    Real phenglei_2d_ssprk3 = 4.0;
    Real phenglei_3d_ssprk3 = 4.0;
    Real scmm6_2d_ssprk3 = 4.0;
    Real scmm6_3d_ssprk3 = 4.0;

    void validate() const;
    [[nodiscard]] Real for_ssprk3(
        AlgorithmProfileKind profile, int dimension) const;
    [[nodiscard]] std::string summary() const;
};

struct ViscousWcnsConfig {
    InviscidWcnsConfig inviscid {};
    TransportConfig transport {};
    ViscousStabilityCoefficients stability {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

class ViscousWcnsSolver {
public:
    ViscousWcnsSolver(
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
        ViscousWcnsConfig config = {});

    void compute_residuals(Real stage_time, int rk_stage = 0);
    void advance(Real time_step, Real initial_time);
    [[nodiscard]] Real global_time_step(Real cfl);
    [[nodiscard]] Real global_residual_l2() const;
    [[nodiscard]] const ReconstructionDiagnostics& reconstruction_diagnostics() const noexcept
    {
        return reconstruction_diagnostics_;
    }
    [[nodiscard]] const RiemannDiagnostics& riemann_diagnostics() const noexcept
    {
        return riemann_diagnostics_;
    }
    [[nodiscard]] std::size_t global_reconstruction_fallback_count() const;
    [[nodiscard]] std::size_t global_riemann_face_count() const;
    [[nodiscard]] std::size_t global_riemann_fallback_count() const;

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
    ViscousWcnsConfig config_;
    SourceTermRegistry source_registry_;
    TransportModel transport_;
    RiemannSolver riemann_ {};
    std::uint64_t version_ = 0;
    ReconstructionDiagnostics reconstruction_diagnostics_ {};
    RiemannDiagnostics riemann_diagnostics_ {};
};

} // namespace wcns
