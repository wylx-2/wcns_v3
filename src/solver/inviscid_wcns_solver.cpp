#include <wcns/solver/inviscid_wcns_solver.hpp>

#include <wcns/solver/time_integrator.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace wcns {

void InviscidWcnsConfig::validate() const
{
    reconstruction.validate();
    source_terms.validate();
}

std::string InviscidWcnsConfig::summary() const
{
    validate();
    return reconstruction.summary() + ';' + boundary.summary() + ';'
        + source_terms.summary();
}

std::string InviscidWcnsConfig::restart_signature() const
{
    return "inviscid_wcns_v1;" + summary();
}

InviscidWcnsSolver::InviscidWcnsSolver(
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
    InviscidWcnsConfig config)
    : mpi_(mpi)
    , local_blocks_(local_blocks)
    , global_mesh_(global_mesh)
    , topology_(topology)
    , state_exchanger_(mpi, topology, distribution_rank_count)
    , metrics_(metrics)
    , boundary_data_(boundary_data)
    , profile_(std::move(profile))
    , gas_(std::move(gas))
    , reference_(std::move(reference))
    , floors_(floors)
    , config_(std::move(config))
    , source_registry_(SourceTermRegistry::create_stage_j(config_.source_terms))
{
    config_.validate();
    floors_.validate();
    ProfileFactory::validate_bundle(profile_.components());
    if (local_blocks_.rank() != mpi_.rank()) {
        throw std::invalid_argument("WCNS local block rank differs from MPI rank");
    }
    for (const auto& block : local_blocks_.blocks()) {
        const auto metric = metrics_.find(block.id());
        if (metric == metrics_.end() || metric->second.profile() != profile_.kind()) {
            throw std::invalid_argument("WCNS metric is missing or has a different profile");
        }
        if (boundary_data_.find(block.id()) == boundary_data_.end()) {
            throw std::invalid_argument("WCNS boundary data is missing for a local block");
        }
    }
}

void InviscidWcnsSolver::compute_residuals(Real stage_time)
{
    if (!std::isfinite(stage_time)) {
        throw std::invalid_argument("WCNS stage time must be finite");
    }
    if (version_ == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("WCNS state version overflow");
    }
    ++version_;
    BlockFieldRegistry conservative_fields(euler_components);
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(
            block, gas_, reference_, floors_);
        conservative_fields.add(block.id(), block.flow.conservative);
    }
    state_exchanger_.exchange(conservative_fields);
    for (const auto& exchange : topology_.exchanges()) {
        if (exchange.receiver_rank != mpi_.rank()) continue;
        auto& receiver = local_blocks_.block(exchange.halo.receiver_block);
        for (const auto& pair : exchange.halo.cell_pairs) {
            update_temperature_primitive_cell(
                receiver, pair.receiver_ghost, gas_, reference_, floors_);
        }
    }

    std::unordered_map<BlockId, InviscidFaceFluxField> fluxes;
    FaceFluxFieldRegistry flux_registry;
    reconstruction_diagnostics_ = {};
    for (auto& block : local_blocks_.blocks()) {
        const auto data = boundary_data_.find(block.id());
        const auto ghost = PhysicalGhostStateOperator::fill(
            block, data->second, gas_, reference_, floors_, version_);
        if (ghost.version != version_) {
            throw std::logic_error("WCNS physical ghost version mismatch");
        }
        auto [iterator, inserted] = fluxes.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(compute_inviscid_face_fluxes(
                block, metrics_.at(block.id()), profile_, config_.reconstruction,
                riemann_, gas_, reference_, floors_, data->second,
                config_.boundary, version_, reconstruction_diagnostics_)));
        if (!inserted) {
            throw std::logic_error("duplicate local WCNS face-flux block");
        }
        flux_registry.add(block.id(), iterator->second);
    }
    const auto current_plan
        = FaceFluxHaloPlan::build(global_mesh_, profile_, version_);
    FaceFluxHaloExchanger(mpi_, current_plan).exchange(flux_registry);
    for (auto& block : local_blocks_.blocks()) {
        compute_wcns_inviscid_residual(
            block, metrics_.at(block.id()), fluxes.at(block.id()), profile_);
        add_source_terms(
            block, metrics_.at(block.id()), source_registry_, stage_time);
    }
}

void InviscidWcnsSolver::advance(Real time_step, Real initial_time)
{
    std::vector<StructuredBlock*> blocks;
    for (auto& block : local_blocks_.blocks()) blocks.push_back(&block);
    advance_ssprk3(blocks, time_step, initial_time,
        [this](Real stage_time) { compute_residuals(stage_time); });
}

Real InviscidWcnsSolver::global_residual_l2() const
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
    if (count <= 0.0) throw PhysicsError("WCNS residual norm has no cells");
    return std::sqrt(mpi_.sum(local_sum) / count);
}

} // namespace wcns
