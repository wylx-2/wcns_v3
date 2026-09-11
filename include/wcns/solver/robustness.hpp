#pragma once

#include <wcns/solver/inviscid_flux.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

struct RobustnessConfig {
    bool enabled = false;
    int max_local_recomputations = 3;
    int max_step_retries = 4;
    Real time_step_reduction = 0.5;
    Real minimum_time_step = 1.0e-12;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct RobustFluxStrategy {
    ReconstructionConfig reconstruction;
    bool force_rusanov = false;
};

class RobustnessLadder {
public:
    [[nodiscard]] static RobustnessLadder build(const ReconstructionConfig& reconstruction,
                                                const RiemannConfig& riemann);

    [[nodiscard]] int maximum_level() const noexcept
    {
        return static_cast<int>(strategies_.size()) - 1;
    }
    [[nodiscard]] const RobustFluxStrategy& strategy(int level) const;
    [[nodiscard]] std::size_t size() const noexcept { return strategies_.size(); }
    [[nodiscard]] std::string summary() const;

private:
    std::vector<RobustFluxStrategy> strategies_;
};

class FaceRobustnessField {
public:
    FaceRobustnessField(Extent3 cell_extent, int dimension, AlgorithmProfileKind profile);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] int halo_layers() const noexcept { return halo_layers_; }
    [[nodiscard]] Field<int>& field(Axis axis);
    [[nodiscard]] const Field<int>& field(Axis axis) const;
    [[nodiscard]] int level(Axis axis, Index3 face) const;
    [[nodiscard]] bool request_next(Axis axis, Index3 face, int maximum_level);
    [[nodiscard]] bool merge_max(Axis axis, Index3 face, int requested_level);

private:
    AlgorithmProfileKind profile_;
    int dimension_ = 0;
    int halo_layers_ = 0;
    Field<int> i_;
    Field<int> j_;
    Field<int> k_;
};

class FaceRobustnessRegistry {
public:
    void add(BlockId block, FaceRobustnessField& field);
    [[nodiscard]] FaceRobustnessField& field(BlockId block) const;

private:
    std::unordered_map<BlockId, FaceRobustnessField*> fields_;
};

using BlockFaceRobustnessMap = std::unordered_map<BlockId, FaceRobustnessField>;
using RobustResidualEvaluator = std::function<void(Real, int, const BlockFaceRobustnessMap&)>;

class FaceRobustnessExchanger {
public:
    FaceRobustnessExchanger(const MpiRuntime& mpi, const FaceFluxHaloPlan& plan)
        : mpi_(mpi)
        , plan_(plan)
    { }

    [[nodiscard]] bool exchange_requests(const FaceRobustnessRegistry& fields,
                                         int rk_stage,
                                         int recomputation_round) const;

private:
    const MpiRuntime& mpi_;
    const FaceFluxHaloPlan& plan_;
};

struct TroubledCell {
    BlockId block = invalid_block_id;
    RankId rank = invalid_rank_id;
    Index3 cell {};
    int rk_stage = 0;
    Real stage_time = 0.0;
    Real density = std::numeric_limits<Real>::quiet_NaN();
    Real pressure = std::numeric_limits<Real>::quiet_NaN();
    Real temperature = std::numeric_limits<Real>::quiet_NaN();
    Real internal_energy = std::numeric_limits<Real>::quiet_NaN();
    std::string reason;
};

struct BlockStateBuffer {
    BlockId block = invalid_block_id;
    Extent3 extent {};
    std::vector<Real> values;
};

using StateSnapshot = std::vector<BlockStateBuffer>;

struct CandidateValidation {
    std::vector<TroubledCell> troubled_cells;
    Real minimum_density = std::numeric_limits<Real>::infinity();
    Real minimum_pressure = std::numeric_limits<Real>::infinity();
    Real minimum_temperature = std::numeric_limits<Real>::infinity();
    Real minimum_internal_energy = std::numeric_limits<Real>::infinity();

    [[nodiscard]] bool valid() const noexcept { return troubled_cells.empty(); }
};

struct RobustnessDiagnostics {
    std::array<std::size_t, 4> face_levels {};
    std::size_t troubled_cells = 0;
    std::size_t local_recomputations = 0;
    std::size_t step_retries = 0;
    Real proposed_time_step = 0.0;
    Real accepted_time_step = 0.0;
    Real minimum_density = std::numeric_limits<Real>::infinity();
    Real minimum_pressure = std::numeric_limits<Real>::infinity();
    Real minimum_temperature = std::numeric_limits<Real>::infinity();
    Real minimum_internal_energy = std::numeric_limits<Real>::infinity();
    std::vector<TroubledCell> failures;

    void observe(const CandidateValidation& validation);
};

[[nodiscard]] StateSnapshot capture_conservative_state(const std::vector<StructuredBlock*>& blocks);

void restore_conservative_state(const std::vector<StructuredBlock*>& blocks,
                                const StateSnapshot& snapshot);

[[nodiscard]] StateSnapshot form_ssprk_candidate(const std::vector<StructuredBlock*>& blocks,
                                                 const StateSnapshot& initial,
                                                 Real initial_weight,
                                                 Real stage_weight,
                                                 Real residual_weight);

[[nodiscard]] CandidateValidation
validate_candidate_state(const StateSnapshot& candidate,
                         const std::vector<StructuredBlock*>& blocks,
                         const GasModel& gas,
                         const ReferenceScales& reference,
                         const NumericalFloors& floors,
                         int rk_stage,
                         Real stage_time);

void commit_candidate_state(const std::vector<StructuredBlock*>& blocks,
                            const StateSnapshot& candidate);

[[nodiscard]] bool request_troubled_cell_support(const StructuredBlock& block,
                                                 const AlgorithmProfile& profile,
                                                 FluxDifferenceMode mode,
                                                 const std::vector<TroubledCell>& troubled_cells,
                                                 FaceRobustnessField& levels,
                                                 int maximum_level);

[[nodiscard]] std::array<std::size_t, 4> count_owned_face_levels(const StructuredBlock& block,
                                                                 const FaceRobustnessField& levels);

[[nodiscard]] Real advance_ssprk3_with_robustness(const MpiRuntime& mpi,
                                                  const std::vector<StructuredBlock*>& blocks,
                                                  const StructuredMesh& global_mesh,
                                                  const AlgorithmProfile& profile,
                                                  FluxDifferenceMode flux_difference,
                                                  const GasModel& gas,
                                                  const ReferenceScales& reference,
                                                  const NumericalFloors& floors,
                                                  const RobustnessConfig& config,
                                                  const RobustnessLadder& ladder,
                                                  Real proposed_time_step,
                                                  Real initial_time,
                                                  const RobustResidualEvaluator& evaluate_residuals,
                                                  RobustnessDiagnostics& diagnostics);

} // namespace wcns
