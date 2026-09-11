#include <wcns/solver/robustness.hpp>

#include <wcns/parallel/mpi_runtime.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <set>
#include <stdexcept>
#include <tuple>

namespace wcns {
namespace {

constexpr int robustness_message_tag_base = 24576;

bool same_strategy(
    const RobustFluxStrategy& lhs,
    const RobustFluxStrategy& rhs)
{
    return lhs.reconstruction.scheme == rhs.reconstruction.scheme
        && lhs.reconstruction.variables == rhs.reconstruction.variables
        && lhs.force_rusanov == rhs.force_rusanov;
}

std::size_t expected_values(Extent3 extent)
{
    return extent.size() * static_cast<std::size_t>(euler_components);
}

void require_matching_snapshot(
    const std::vector<StructuredBlock*>& blocks,
    const StateSnapshot& snapshot)
{
    if (blocks.empty() || blocks.size() != snapshot.size()) {
        throw std::invalid_argument("state snapshot does not match block set");
    }
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        if (blocks[index] == nullptr
            || blocks[index]->id() != snapshot[index].block
            || blocks[index]->cell_extent() != snapshot[index].extent
            || snapshot[index].values.size()
                != expected_values(snapshot[index].extent)) {
            throw std::invalid_argument("state snapshot block metadata mismatch");
        }
    }
}

std::size_t candidate_offset(Extent3 extent, Index3 cell, int component)
{
    return (((static_cast<std::size_t>(cell.k)
                  * static_cast<std::size_t>(extent.nj)
              + static_cast<std::size_t>(cell.j))
                 * static_cast<std::size_t>(extent.ni)
             + static_cast<std::size_t>(cell.i))
            * static_cast<std::size_t>(euler_components))
        + static_cast<std::size_t>(component);
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("robustness request message exceeds MPI int count");
    }
    return static_cast<int>(count);
}
#endif

} // namespace

void RobustnessConfig::validate() const
{
    if (max_local_recomputations < 0 || max_step_retries < 0) {
        throw std::invalid_argument("robustness retry limits must be non-negative");
    }
    if (!std::isfinite(time_step_reduction)
        || time_step_reduction <= 0.0 || time_step_reduction >= 1.0) {
        throw std::invalid_argument("robustness time-step reduction must lie in (0,1)");
    }
    if (!std::isfinite(minimum_time_step) || minimum_time_step <= 0.0) {
        throw std::invalid_argument("robustness minimum time step must be positive and finite");
    }
}

std::string RobustnessConfig::summary() const
{
    validate();
    std::ostringstream result;
    result << "robustness(enabled=" << (enabled ? "true" : "false")
           << ",max_local_recomputations=" << max_local_recomputations
           << ",max_step_retries=" << max_step_retries
           << ",time_step_reduction=" << std::setprecision(17)
           << time_step_reduction
           << ",minimum_time_step=" << minimum_time_step << ')';
    return result.str();
}

std::string RobustnessConfig::restart_signature() const
{
    return "robustness_v1;" + summary();
}

RobustnessLadder RobustnessLadder::build(
    const ReconstructionConfig& reconstruction,
    const RiemannConfig& riemann)
{
    reconstruction.validate();
    riemann.validate();
    RobustnessLadder result;
    const auto append = [&](RobustFluxStrategy strategy) {
        strategy.reconstruction.validate();
        if (result.strategies_.empty()
            || !same_strategy(result.strategies_.back(), strategy)) {
            result.strategies_.push_back(std::move(strategy));
        }
    };
    append({reconstruction, false});
    auto primitive = reconstruction;
    primitive.variables = ReconstructionVariables::Primitive;
    append({primitive, false});
    auto linear = primitive;
    linear.scheme = "linear5";
    append({linear, false});
    auto first_order = linear;
    first_order.scheme = "zero_order";
    first_order.variables = ReconstructionVariables::Conservative;
    append({first_order, riemann.scheme != "rusanov"});
    if (result.strategies_.empty()) {
        throw std::logic_error("robustness ladder is empty");
    }
    return result;
}

const RobustFluxStrategy& RobustnessLadder::strategy(int level) const
{
    if (level < 0 || level >= static_cast<int>(strategies_.size())) {
        throw std::out_of_range("robustness level is outside the effective ladder");
    }
    return strategies_[static_cast<std::size_t>(level)];
}

std::string RobustnessLadder::summary() const
{
    std::ostringstream result;
    for (std::size_t index = 0; index < strategies_.size(); ++index) {
        if (index != 0) result << "->";
        result << index << ':' << strategies_[index].reconstruction.scheme << ':';
        switch (strategies_[index].reconstruction.variables) {
        case ReconstructionVariables::Conservative: result << "conservative"; break;
        case ReconstructionVariables::Primitive: result << "primitive"; break;
        case ReconstructionVariables::Characteristic: result << "characteristic"; break;
        }
        if (strategies_[index].force_rusanov) result << ":rusanov";
    }
    return result.str();
}

FaceRobustnessField::FaceRobustnessField(
    Extent3 cells,
    int dimension,
    AlgorithmProfileKind profile)
    : profile_(profile)
    , dimension_(dimension)
    , halo_layers_(profile == AlgorithmProfileKind::PhengleiWcns ? 1 : 2)
    , i_({cells.ni + 1, cells.nj, cells.nk}, 1, halo_layers_, 0)
    , j_({cells.ni, cells.nj + 1, cells.nk}, 1, halo_layers_, 0)
    , k_({cells.ni, cells.nj, cells.nk + 1}, 1, halo_layers_, 0)
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("face robustness field dimension must be 2 or 3");
    }
}

Field<int>& FaceRobustnessField::field(Axis axis)
{
    if (axis == Axis::I) return i_;
    if (axis == Axis::J) return j_;
    if (axis == Axis::K && dimension_ == 3) return k_;
    throw std::out_of_range("face robustness axis is not active");
}

const Field<int>& FaceRobustnessField::field(Axis axis) const
{
    return const_cast<FaceRobustnessField*>(this)->field(axis);
}

int FaceRobustnessField::level(Axis axis, Index3 face) const
{
    return field(axis)(face.i, face.j, face.k, 0);
}

bool FaceRobustnessField::request_next(
    Axis axis,
    Index3 face,
    int maximum_level)
{
    if (maximum_level < 0) {
        throw std::invalid_argument("maximum robustness level must be non-negative");
    }
    auto& value = field(axis)(face.i, face.j, face.k, 0);
    if (value < 0 || value > maximum_level) {
        throw std::logic_error("face robustness level is invalid");
    }
    if (value == maximum_level) return false;
    ++value;
    return true;
}

bool FaceRobustnessField::merge_max(
    Axis axis,
    Index3 face,
    int requested_level)
{
    if (requested_level < 0) {
        throw std::invalid_argument("requested robustness level is negative");
    }
    auto& value = field(axis)(face.i, face.j, face.k, 0);
    if (requested_level <= value) return false;
    value = requested_level;
    return true;
}

void FaceRobustnessRegistry::add(BlockId block, FaceRobustnessField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("face robustness registry has invalid or duplicate block");
    }
}

FaceRobustnessField& FaceRobustnessRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("face robustness field is not registered");
    }
    return *iterator->second;
}

bool FaceRobustnessExchanger::exchange_requests(
    const FaceRobustnessRegistry& fields,
    int rk_stage,
    int recomputation_round) const
{
    if (rk_stage < 1 || rk_stage > 3 || recomputation_round < 1) {
        throw std::invalid_argument("robustness exchange stage/round is invalid");
    }
    struct Pending {
        const FaceFluxExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };
    const RankId rank = mpi_.rank();
    bool changed = false;
    std::vector<Pending> receives;
    std::vector<Pending> sends;
    for (const auto& descriptor : plan_.exchanges()) {
        const std::size_t count = 3 + descriptor.pairs.size();
        if (descriptor.receiver_rank == rank && descriptor.donor_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            auto& donor = fields.field(descriptor.donor_block);
            for (const auto& pair : descriptor.pairs) {
                changed = donor.merge_max(
                    descriptor.donor_axis,
                    pair.donor,
                    receiver.level(descriptor.receiver_axis, pair.receiver)) || changed;
            }
        } else if (descriptor.receiver_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            Pending pending {&descriptor, std::vector<Real>(count)};
            pending.values[0] = static_cast<Real>(descriptor.version);
            pending.values[1] = static_cast<Real>(rk_stage);
            pending.values[2] = static_cast<Real>(recomputation_round);
            std::size_t offset = 3;
            for (const auto& pair : descriptor.pairs) {
                pending.values[offset++] = static_cast<Real>(
                    receiver.level(descriptor.receiver_axis, pair.receiver));
            }
            sends.push_back(std::move(pending));
        } else if (descriptor.donor_rank == rank) {
            static_cast<void>(fields.field(descriptor.donor_block));
            receives.push_back({&descriptor, std::vector<Real>(count)});
        }
    }

#if WCNS_HAS_MPI
    std::vector<MPI_Request> requests(receives.size() + sends.size(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives) {
        check_mpi(MPI_Irecv(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->receiver_rank,
            pending.descriptor->message_tag(robustness_message_tag_base),
            mpi_.communicator(), &requests[request++]),
            "MPI_Irecv robustness request");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->donor_rank,
            pending.descriptor->message_tag(robustness_message_tag_base),
            mpi_.communicator(), &requests[request++]),
            "MPI_Isend robustness request");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall robustness request");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote robustness exchange requires WCNS_ENABLE_MPI");
    }
#endif

    for (const auto& pending : receives) {
        if (pending.values.size() < 3
            || pending.values[0] != static_cast<Real>(pending.descriptor->version)
            || pending.values[1] != static_cast<Real>(rk_stage)
            || pending.values[2] != static_cast<Real>(recomputation_round)) {
            throw MpiError("robustness request message header mismatch");
        }
        auto& donor = fields.field(pending.descriptor->donor_block);
        std::size_t offset = 3;
        for (const auto& pair : pending.descriptor->pairs) {
            const Real encoded = pending.values[offset++];
            const int level = static_cast<int>(encoded);
            if (!std::isfinite(encoded) || encoded != static_cast<Real>(level)
                || level < 0) {
                throw MpiError("robustness request contains an invalid level");
            }
            changed = donor.merge_max(
                pending.descriptor->donor_axis, pair.donor, level) || changed;
        }
    }
    return !mpi_.all_true(!changed);
}

StateSnapshot capture_conservative_state(
    const std::vector<StructuredBlock*>& blocks)
{
    if (blocks.empty()) {
        throw std::invalid_argument("cannot capture an empty block set");
    }
    StateSnapshot result;
    result.reserve(blocks.size());
    for (const auto* block : blocks) {
        if (block == nullptr) throw std::invalid_argument("block pointer is null");
        BlockStateBuffer buffer {block->id(), block->cell_extent(), {}};
        buffer.values.reserve(expected_values(buffer.extent));
        for (int k = 0; k < buffer.extent.nk; ++k) {
            for (int j = 0; j < buffer.extent.nj; ++j) {
                for (int i = 0; i < buffer.extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        buffer.values.push_back(
                            block->flow.conservative(i, j, k, component));
                    }
                }
            }
        }
        result.push_back(std::move(buffer));
    }
    return result;
}

void restore_conservative_state(
    const std::vector<StructuredBlock*>& blocks,
    const StateSnapshot& snapshot)
{
    require_matching_snapshot(blocks, snapshot);
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        std::size_t offset = 0;
        const auto extent = snapshot[index].extent;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        blocks[index]->flow.conservative(i, j, k, component)
                            = snapshot[index].values[offset++];
                    }
                }
            }
        }
    }
}

StateSnapshot form_ssprk_candidate(
    const std::vector<StructuredBlock*>& blocks,
    const StateSnapshot& initial,
    Real initial_weight,
    Real stage_weight,
    Real residual_weight)
{
    require_matching_snapshot(blocks, initial);
    if (!std::isfinite(initial_weight) || !std::isfinite(stage_weight)
        || !std::isfinite(residual_weight)) {
        throw std::invalid_argument("SSPRK candidate weights must be finite");
    }
    StateSnapshot result = initial;
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        std::size_t offset = 0;
        const auto extent = initial[index].extent;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        result[index].values[offset]
                            = initial_weight * initial[index].values[offset]
                            + stage_weight
                                * blocks[index]->flow.conservative(i, j, k, component)
                            + residual_weight
                                * blocks[index]->flow.residual(i, j, k, component);
                        ++offset;
                    }
                }
            }
        }
    }
    return result;
}

CandidateValidation validate_candidate_state(
    const StateSnapshot& candidate,
    const std::vector<StructuredBlock*>& blocks,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int rk_stage,
    Real stage_time)
{
    require_matching_snapshot(blocks, candidate);
    floors.validate();
    if (rk_stage < 1 || rk_stage > 3 || !std::isfinite(stage_time)) {
        throw std::invalid_argument("candidate validation stage/time is invalid");
    }
    CandidateValidation result;
    const Real gamma_minus_one = gas.gamma() - 1.0;
    const Real temperature_scale = gas.gamma() * reference.mach() * reference.mach();
    for (std::size_t block_index = 0; block_index < blocks.size(); ++block_index) {
        const auto& values = candidate[block_index].values;
        const auto extent = candidate[block_index].extent;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    const auto base = candidate_offset(extent, cell, 0);
                    bool components_finite = true;
                    for (int component = 0; component < euler_components; ++component) {
                        components_finite = components_finite
                            && std::isfinite(values[base + static_cast<std::size_t>(component)]);
                    }
                    const Real density = values[base];
                    Real pressure = std::numeric_limits<Real>::quiet_NaN();
                    Real temperature = std::numeric_limits<Real>::quiet_NaN();
                    Real internal_energy = std::numeric_limits<Real>::quiet_NaN();
                    if (components_finite && std::isfinite(density) && density > 0.0) {
                        const Real mx = values[base + 1];
                        const Real my = values[base + 2];
                        const Real mz = values[base + 3];
                        const Real internal_density = values[base + 4]
                            - (mx * mx + my * my + mz * mz) / (2.0 * density);
                        pressure = gamma_minus_one * internal_density;
                        internal_energy = internal_density / density;
                        temperature = temperature_scale * pressure / density;
                    }
                    result.minimum_density = std::min(result.minimum_density, density);
                    result.minimum_pressure = std::min(result.minimum_pressure, pressure);
                    result.minimum_temperature = std::min(result.minimum_temperature, temperature);
                    result.minimum_internal_energy = std::min(
                        result.minimum_internal_energy, internal_energy);
                    std::string reason;
                    if (!components_finite) reason = "non_finite_conservative";
                    else if (!std::isfinite(density) || density < floors.density)
                        reason = "density_floor";
                    else if (!std::isfinite(pressure) || pressure < floors.pressure)
                        reason = "pressure_floor";
                    else if (!std::isfinite(temperature)
                        || temperature < floors.temperature)
                        reason = "temperature_floor";
                    else if (!std::isfinite(internal_energy) || internal_energy <= 0.0)
                        reason = "internal_energy";
                    if (!reason.empty()) {
                        result.troubled_cells.push_back({
                            blocks[block_index]->id(), blocks[block_index]->owner_rank(),
                            cell, rk_stage, stage_time, density, pressure, temperature,
                            internal_energy, std::move(reason)});
                    }
                }
            }
        }
    }
    return result;
}

void commit_candidate_state(
    const std::vector<StructuredBlock*>& blocks,
    const StateSnapshot& candidate)
{
    restore_conservative_state(blocks, candidate);
}

bool request_troubled_cell_support(
    const StructuredBlock& block,
    const AlgorithmProfile& profile,
    FluxDifferenceMode mode,
    const std::vector<TroubledCell>& troubled_cells,
    FaceRobustnessField& levels,
    int maximum_level)
{
    if (levels.profile() != profile.kind()
        || levels.dimension() != block.cell_dimension()) {
        throw ProfileError("troubled-cell support uses a mismatched face field");
    }
    bool changed = false;
    using FaceKey = std::tuple<int, int, int, int>;
    using CellKey = std::tuple<int, int, int>;
    std::set<FaceKey> direct_faces;
    const auto add_cell_support = [&](Index3 cell, std::set<FaceKey>& faces) {
        for (int logical = 0; logical < block.cell_dimension(); ++logical) {
            const auto axis = static_cast<Axis>(logical);
            for (const auto& [face_index, coefficient] : inviscid_residual_stencil(
                     block, profile, mode, axis, cell)) {
                if (coefficient == 0.0) continue;
                auto face = cell;
                face[static_cast<std::size_t>(axis)] = face_index;
                faces.emplace(
                    logical, face.i, face.j, face.k);
            }
        }
    };
    for (const auto& troubled : troubled_cells) {
        if (troubled.block == block.id()) {
            add_cell_support(troubled.cell, direct_faces);
        }
    }
    if (direct_faces.empty()) return false;

    // One transpose-support guard layer makes every cell whose residual is
    // changed by a direct downgrade use a coherent downgraded face patch.
    std::set<CellKey> guard_cells;
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                std::set<FaceKey> support;
                add_cell_support({i, j, k}, support);
                bool intersects = false;
                for (const auto& face : support) {
                    if (direct_faces.find(face) != direct_faces.end()) {
                        intersects = true;
                        break;
                    }
                }
                if (intersects) guard_cells.emplace(i, j, k);
            }
        }
    }
    std::set<FaceKey> requested_faces = direct_faces;
    for (const auto& [i, j, k] : guard_cells) {
        add_cell_support({i, j, k}, requested_faces);
    }
    for (const auto& [logical, i, j, k] : requested_faces) {
        changed = levels.request_next(
            static_cast<Axis>(logical), {i, j, k}, maximum_level) || changed;
    }
    return changed;
}

std::array<std::size_t, 4> count_owned_face_levels(
    const StructuredBlock& block,
    const FaceRobustnessField& levels)
{
    std::array<std::size_t, 4> result {};
    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
        const auto axis = static_cast<Axis>(logical);
        const auto extent = levels.field(axis).interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 face {i, j, k};
                    if (is_non_owned_connection_face(block, axis, face)) continue;
                    const int level = levels.level(axis, face);
                    if (level < 0 || level >= static_cast<int>(result.size())) {
                        throw std::logic_error("owned face robustness level is invalid");
                    }
                    ++result[static_cast<std::size_t>(level)];
                }
            }
        }
    }
    return result;
}

void RobustnessDiagnostics::observe(const CandidateValidation& validation)
{
    troubled_cells += validation.troubled_cells.size();
    minimum_density = std::min(minimum_density, validation.minimum_density);
    minimum_pressure = std::min(minimum_pressure, validation.minimum_pressure);
    minimum_temperature = std::min(minimum_temperature, validation.minimum_temperature);
    minimum_internal_energy = std::min(
        minimum_internal_energy, validation.minimum_internal_energy);
    failures.insert(
        failures.end(), validation.troubled_cells.begin(),
        validation.troubled_cells.end());
}

Real advance_ssprk3_with_robustness(
    const MpiRuntime& mpi,
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
    RobustnessDiagnostics& diagnostics)
{
    config.validate();
    if (!config.enabled || blocks.empty() || !evaluate_residuals
        || !std::isfinite(proposed_time_step) || proposed_time_step <= 0.0
        || !std::isfinite(initial_time)) {
        throw std::invalid_argument("robust SSPRK3 inputs are invalid");
    }
    diagnostics = {};
    diagnostics.proposed_time_step = proposed_time_step;
    const StateSnapshot initial = capture_conservative_state(blocks);
    Real time_step = proposed_time_step;
    std::uint64_t request_version = 1;
    std::array<std::size_t, 4> attempt_face_levels {};

    const auto run_stage = [&](int rk_stage, Real stage_time,
                               Real initial_weight, Real stage_weight,
                               Real residual_weight) {
        BlockFaceRobustnessMap levels;
        FaceRobustnessRegistry registry;
        for (const auto* block : blocks) {
            auto [iterator, inserted] = levels.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block->id()),
                std::forward_as_tuple(
                    block->cell_extent(), block->cell_dimension(), profile.kind()));
            if (!inserted) throw std::logic_error("duplicate robustness block");
            registry.add(block->id(), iterator->second);
        }
        for (int round = 0;; ++round) {
            evaluate_residuals(stage_time, rk_stage, levels);
            auto candidate = form_ssprk_candidate(
                blocks, initial, initial_weight, stage_weight, residual_weight);
            const auto validation = validate_candidate_state(
                candidate, blocks, gas, reference, floors, rk_stage, stage_time);
            diagnostics.observe(validation);
            if (mpi.all_true(validation.valid())) {
                commit_candidate_state(blocks, candidate);
                for (const auto* block : blocks) {
                    const auto counts = count_owned_face_levels(
                        *block, levels.at(block->id()));
                    for (std::size_t level = 0;
                         level < attempt_face_levels.size(); ++level) {
                        attempt_face_levels[level] += counts[level];
                    }
                }
                return true;
            }
            if (round >= config.max_local_recomputations) return false;
            bool local_changed = false;
            for (const auto* block : blocks) {
                local_changed = request_troubled_cell_support(
                    *block, profile, flux_difference,
                    validation.troubled_cells, levels.at(block->id()),
                    ladder.maximum_level()) || local_changed;
            }
            if (request_version == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("robustness request version overflow");
            }
            const auto plan = FaceFluxHaloPlan::build(
                global_mesh, profile, request_version++);
            const bool exchange_changed = FaceRobustnessExchanger(mpi, plan)
                .exchange_requests(registry, rk_stage, round + 1);
            const bool global_changed = !mpi.all_true(
                !(local_changed || exchange_changed));
            if (!global_changed) return false;
            ++diagnostics.local_recomputations;
        }
    };

    for (int retry = 0; retry <= config.max_step_retries; ++retry) {
        restore_conservative_state(blocks, initial);
        attempt_face_levels.fill(0);
        bool accepted = run_stage(
            1, initial_time, 1.0, 0.0, time_step);
        if (accepted) {
            accepted = run_stage(
                2, initial_time + time_step,
                0.75, 0.25, 0.25 * time_step);
        }
        if (accepted) {
            accepted = run_stage(
                3, initial_time + 0.5 * time_step,
                1.0 / 3.0, 2.0 / 3.0, 2.0 * time_step / 3.0);
        }
        if (accepted) {
            diagnostics.face_levels = attempt_face_levels;
            diagnostics.accepted_time_step = time_step;
            return time_step;
        }
        restore_conservative_state(blocks, initial);
        if (retry == config.max_step_retries
            || time_step * config.time_step_reduction
                < config.minimum_time_step) {
            std::ostringstream message;
            message << "robust SSPRK3 exhausted retries at dt="
                    << std::setprecision(17) << time_step;
            if (!diagnostics.failures.empty()) {
                const auto& first = diagnostics.failures.front();
                message << "; first troubled cell rank=" << first.rank
                        << " block=" << first.block
                        << " cell=(" << first.cell.i << ',' << first.cell.j
                        << ',' << first.cell.k << ") stage=" << first.rk_stage
                        << " time=" << first.stage_time
                        << " reason=" << first.reason;
            }
            throw PhysicsError(message.str());
        }
        time_step *= config.time_step_reduction;
        ++diagnostics.step_retries;
    }
    throw std::logic_error("robust SSPRK3 retry loop did not terminate");
}

} // namespace wcns
