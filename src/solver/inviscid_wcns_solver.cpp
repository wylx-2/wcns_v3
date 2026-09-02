#include <wcns/solver/inviscid_wcns_solver.hpp>

#include <wcns/solver/time_integrator.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace wcns {

namespace {

constexpr std::size_t maximum_exact_diagnostic_count
    = static_cast<std::size_t>(9007199254740992ULL);

bool same_floors(const NumericalFloors& lhs, const NumericalFloors& rhs)
{
    return lhs.density == rhs.density && lhs.pressure == rhs.pressure
        && lhs.temperature == rhs.temperature
        && lhs.jacobian_absolute == rhs.jacobian_absolute
        && lhs.jacobian_relative == rhs.jacobian_relative
        && lhs.face_area_absolute == rhs.face_area_absolute
        && lhs.face_area_relative == rhs.face_area_relative
        && lhs.reconstruction_scale == rhs.reconstruction_scale
        && lhs.reconstruction_epsilon == rhs.reconstruction_epsilon;
}

std::size_t global_diagnostic_count(
    const MpiRuntime& mpi, std::size_t local, const char* label)
{
    if (local > maximum_exact_diagnostic_count) {
        throw std::overflow_error(std::string(label) + " exceeds exact MPI reduction range");
    }
    const Real global = mpi.sum(static_cast<Real>(local));
    if (!std::isfinite(global) || global < 0.0
        || global > static_cast<Real>(maximum_exact_diagnostic_count)) {
        throw std::overflow_error(std::string(label) + " global reduction is invalid");
    }
    return static_cast<std::size_t>(global);
}

std::array<Real, 3> metric_area_vector(
    const MetricField& metric,
    Axis axis,
    Index3 face)
{
    const FaceAreaVectors* vectors = nullptr;
    switch (axis) {
    case Axis::I: vectors = &metric.i_faces(); break;
    case Axis::J: vectors = &metric.j_faces(); break;
    case Axis::K: vectors = &metric.k_faces(); break;
    }
    return {{
        vectors->x(face.i, face.j, face.k),
        vectors->y(face.i, face.j, face.k),
        vectors->z(face.i, face.j, face.k),
    }};
}

} // namespace

void InviscidWcnsConfig::validate() const
{
    reconstruction.validate();
    riemann.validate();
    source_terms.validate();
}

std::string InviscidWcnsConfig::summary() const
{
    validate();
    return reconstruction.summary() + ';' + riemann.summary() + ';'
        + boundary.summary() + ';' + source_terms.summary();
}

std::string InviscidWcnsConfig::restart_signature() const
{
    return "inviscid_wcns_v2;" + summary();
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
    const auto riemann_registry
        = RiemannSolverRegistry::with_builtins(config_.riemann.parameters);
    config_.riemann.validate(riemann_registry);
    riemann_ = RiemannSolver(
        config_.riemann.scheme, riemann_registry, config_.riemann.parameters);
    floors_.validate();
    if (!same_floors(config_.reconstruction.floors, floors_)) {
        throw std::invalid_argument(
            "WCNS reconstruction and solver must share one NumericalFloors instance");
    }
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

void InviscidWcnsSolver::compute_residuals(Real stage_time, int rk_stage)
{
    if (!std::isfinite(stage_time)) {
        throw std::invalid_argument("WCNS stage time must be finite");
    }
    if (rk_stage < 0 || rk_stage > 3) {
        throw std::invalid_argument("WCNS RK stage must lie in [0,3]");
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
    riemann_diagnostics_ = {};
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
                config_.boundary, version_, reconstruction_diagnostics_,
                &riemann_diagnostics_, rk_stage)));
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
    int rk_stage = 0;
    advance_ssprk3(blocks, time_step, initial_time,
        [this, &rk_stage](Real stage_time) {
            compute_residuals(stage_time, ++rk_stage);
        });
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(
            block, gas_, reference_, floors_);
    }
}

Real InviscidWcnsSolver::global_time_step(Real cfl)
{
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw std::invalid_argument("WCNS CFL must be positive and finite");
    }
    Real local_minimum = std::numeric_limits<Real>::infinity();
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(
            block, gas_, reference_, floors_);
        const auto& metric = metrics_.at(block.id());
        const auto cells = block.cell_extent();
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    TemperaturePrimitiveState state {};
                    for (int component = 0; component < fluid_components; ++component) {
                        state[static_cast<std::size_t>(component)]
                            = block.flow.temperature_primitive(
                                i, j, k, component);
                    }
                    const std::array<Real, 3> velocity {{
                        state[temperature_velocity_x],
                        state[temperature_velocity_y],
                        state[temperature_velocity_z],
                    }};
                    const Real sound = thermodynamic_sound_speed(
                        state,
                        gas_,
                        reference_,
                        floors_,
                        block.cell_dimension());
                    Real spectral_sum = 0.0;
                    for (int logical = 0;
                         logical < block.cell_dimension();
                         ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        for (int side = 0; side <= 1; ++side) {
                            Index3 face {i, j, k};
                            face[static_cast<std::size_t>(axis)] += side;
                            const auto area = metric_area_vector(
                                metric, axis, face);
                            const Real magnitude = std::sqrt(
                                area[0] * area[0]
                                + area[1] * area[1]
                                + area[2] * area[2]);
                            spectral_sum += std::abs(
                                velocity[0] * area[0]
                                + velocity[1] * area[1]
                                + velocity[2] * area[2])
                                + sound * magnitude;
                        }
                    }
                    const Real jacobian = metric.jacobian()(i, j, k);
                    const Real denominator = spectral_sum / (2.0 * jacobian);
                    if (!std::isfinite(denominator) || denominator <= 0.0) {
                        throw PhysicsError(
                            "WCNS time-step denominator is invalid");
                    }
                    local_minimum = std::min(
                        local_minimum, cfl / denominator);
                }
            }
        }
    }
    const Real result = mpi_.min(local_minimum);
    if (!std::isfinite(result) || result <= 0.0) {
        throw PhysicsError("WCNS global time step is invalid");
    }
    return result;
}

std::size_t InviscidWcnsSolver::global_reconstruction_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, reconstruction_diagnostics_.fallback_events.size(),
        "reconstruction fallback count");
}

std::size_t InviscidWcnsSolver::global_riemann_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.fallback_count(),
        "Riemann fallback count");
}

std::size_t InviscidWcnsSolver::global_riemann_face_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.total_faces, "Riemann face count");
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
