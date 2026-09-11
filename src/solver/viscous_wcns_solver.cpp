#include <wcns/solver/viscous_wcns_solver.hpp>

#include <wcns/solver/time_integrator.hpp>

#include <algorithm>
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
        && lhs.temperature == rhs.temperature && lhs.jacobian_absolute == rhs.jacobian_absolute
        && lhs.jacobian_relative == rhs.jacobian_relative
        && lhs.face_area_absolute == rhs.face_area_absolute
        && lhs.face_area_relative == rhs.face_area_relative
        && lhs.reconstruction_scale == rhs.reconstruction_scale
        && lhs.reconstruction_epsilon == rhs.reconstruction_epsilon;
}

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

const FaceAreaVectors& faces(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(const MetricField& metric, Axis axis, Index3 face)
{
    const auto& field = faces(metric, axis);
    return {{
        field.x(face.i, face.j, face.k),
        field.y(face.i, face.j, face.k),
        field.z(face.i, face.j, face.k),
    }};
}

Real area_squared(const std::array<Real, 3>& area)
{
    return area[0] * area[0] + area[1] * area[1] + area[2] * area[2];
}

std::size_t global_diagnostic_count(const MpiRuntime& mpi, std::size_t local, const char* label)
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

} // namespace

void ViscousStabilityCoefficients::validate() const
{
    if (!positive_finite(phenglei_2d_ssprk3) || !positive_finite(phenglei_3d_ssprk3)
        || !positive_finite(scmm6_2d_ssprk3) || !positive_finite(scmm6_3d_ssprk3)) {
        throw std::invalid_argument("all profile/dimension SSPRK3 viscous stability coefficients "
                                    "must be positive and finite");
    }
}

Real ViscousStabilityCoefficients::for_ssprk3(AlgorithmProfileKind profile, int dimension) const
{
    validate();
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument(
            "viscous SSPRK3 stability coefficient requires dimension 2 or 3");
    }
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return dimension == 2 ? phenglei_2d_ssprk3 : phenglei_3d_ssprk3;
    }
    if (profile == AlgorithmProfileKind::Scmm6Wcns) {
        return dimension == 2 ? scmm6_2d_ssprk3 : scmm6_3d_ssprk3;
    }
    throw std::invalid_argument("unknown viscous SSPRK3 algorithm profile");
}

std::string ViscousStabilityCoefficients::summary() const
{
    validate();
    return "Cv_phenglei_2d_ssprk3=" + std::to_string(phenglei_2d_ssprk3)
        + " Cv_phenglei_3d_ssprk3=" + std::to_string(phenglei_3d_ssprk3)
        + " Cv_scmm6_2d_ssprk3=" + std::to_string(scmm6_2d_ssprk3)
        + " Cv_scmm6_3d_ssprk3=" + std::to_string(scmm6_3d_ssprk3);
}

void ViscousWcnsConfig::validate() const
{
    inviscid.validate();
    transport.validate();
    stability.validate();
}

std::string ViscousWcnsConfig::summary() const
{
    validate();
    return inviscid.summary() + ';' + transport.summary() + "; " + stability.summary();
}

std::string ViscousWcnsConfig::restart_signature() const
{
    return "viscous_wcns_v1;" + summary();
}

ViscousWcnsSolver::ViscousWcnsSolver(const MpiRuntime& mpi,
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
                                     ViscousWcnsConfig config)
    : mpi_(mpi)
    , local_blocks_(local_blocks)
    , global_mesh_(global_mesh)
    , topology_(topology)
    , state_exchanger_(mpi, topology, distribution_rank_count, euler_components)
    , metrics_(metrics)
    , boundary_data_(boundary_data)
    , profile_(std::move(profile))
    , gas_(std::move(gas))
    , reference_(std::move(reference))
    , floors_(floors)
    , config_(std::move(config))
    , source_registry_(SourceTermRegistry::create_stage_j(config_.inviscid.source_terms))
    , transport_(config_.transport)
    , inviscid_flux_plan_(FaceFluxHaloPlan::build(global_mesh_, profile_, 1))
    , inviscid_flux_exchanger_(mpi_, inviscid_flux_plan_)
    , operand_plan_(GradientOperandFaceHaloPlan::build(global_mesh_, profile_, 1))
    , operand_exchanger_(mpi_, operand_plan_)
    , gradient_plan_(GradientHaloPlan::build(global_mesh_, topology_, profile_, 1))
    , gradient_exchanger_(mpi_, gradient_plan_)
    , viscous_flux_plan_(ViscousFaceFluxHaloPlan::build(global_mesh_, profile_, 1))
    , viscous_flux_exchanger_(mpi_, viscous_flux_plan_)
{
    config_.validate();
    const auto riemann_registry
        = RiemannSolverRegistry::with_builtins(config_.inviscid.riemann.parameters);
    config_.inviscid.riemann.validate(riemann_registry);
    riemann_ = RiemannSolver(
        config_.inviscid.riemann.scheme, riemann_registry, config_.inviscid.riemann.parameters);
    robust_riemann_
        = RiemannSolver("rusanov", riemann_registry, config_.inviscid.riemann.parameters);
    robustness_ladder_
        = RobustnessLadder::build(config_.inviscid.reconstruction, config_.inviscid.riemann);
    floors_.validate();
    if (!same_floors(config_.inviscid.reconstruction.floors, floors_)) {
        throw std::invalid_argument("viscous WCNS reconstruction and solver floors differ");
    }
    ProfileFactory::validate_bundle(profile_.components());
    if (local_blocks_.rank() != mpi_.rank()) {
        throw std::invalid_argument("viscous WCNS local block rank differs from MPI rank");
    }
    for (const auto& block : local_blocks_.blocks()) {
        const auto metric = metrics_.find(block.id());
        if (metric == metrics_.end() || metric->second.profile() != profile_.kind()) {
            throw std::invalid_argument("viscous WCNS metric is missing or incompatible");
        }
        if (boundary_data_.find(block.id()) == boundary_data_.end()) {
            throw std::invalid_argument("viscous WCNS boundary data is missing");
        }
        for (int logical = 0; logical < block.cell_dimension(); ++logical) {
            static_cast<void>(cached_line_operators(
                profile_, block.cell_extent()[static_cast<std::size_t>(logical)]));
        }
    }
    const auto local_count = local_blocks_.blocks().size();
    block_workspace_.reserve(local_count);
    inviscid_flux_workspace_.reserve(local_count);
    operand_workspace_.reserve(local_count);
    gradient_workspace_.reserve(local_count);
    viscous_flux_workspace_.reserve(local_count);
    for (auto& block : local_blocks_.blocks()) {
        block_workspace_.push_back(&block);
        const auto inviscid = inviscid_flux_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto operand = operand_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto gradient = gradient_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto viscous = viscous_flux_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        if (!inviscid.second || !operand.second || !gradient.second || !viscous.second) {
            throw std::logic_error("duplicate viscous workspace block");
        }
        inviscid_flux_registry_.add(block.id(), inviscid.first->second);
        operand_registry_.add(block.id(), operand.first->second);
        gradient_registry_.add(block.id(), gradient.first->second);
        viscous_flux_registry_.add(block.id(), viscous.first->second);
    }
}

void ViscousWcnsSolver::compute_residuals(Real stage_time, int rk_stage)
{
    compute_residuals_impl(stage_time, rk_stage, nullptr);
}

void ViscousWcnsSolver::compute_residuals_impl(Real stage_time,
                                               int rk_stage,
                                               const BlockFaceRobustnessMap* robustness_levels)
{
    if (!std::isfinite(stage_time)) {
        throw std::invalid_argument("viscous WCNS stage time must be finite");
    }
    if (rk_stage < 0 || rk_stage > 3) {
        throw std::invalid_argument("viscous WCNS RK stage must lie in [0,3]");
    }
    if (version_ == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("viscous WCNS state version overflow");
    }
    ++version_;
    BlockFieldRegistry conservative_fields(euler_components);
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
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

    reconstruction_diagnostics_ = {};
    riemann_diagnostics_ = {};
    for (auto& block : local_blocks_.blocks()) {
        const auto& data = boundary_data_.at(block.id());
        const auto ghost = PhysicalGhostStateOperator::fill(
            block, data, gas_, reference_, floors_, version_, stage_time);
        if (ghost.version != version_) {
            throw std::logic_error("viscous WCNS physical ghost version mismatch");
        }
        compute_inviscid_face_fluxes_into(
            inviscid_flux_workspace_.at(block.id()),
            block,
            metrics_.at(block.id()),
            profile_,
            config_.inviscid.reconstruction,
            riemann_,
            gas_,
            reference_,
            floors_,
            data,
            config_.inviscid.boundary,
            version_,
            reconstruction_diagnostics_,
            &riemann_diagnostics_,
            rk_stage,
            stage_time,
            robustness_levels == nullptr ? nullptr : &robustness_levels->at(block.id()),
            robustness_levels == nullptr ? nullptr : &robustness_ladder_,
            robustness_levels == nullptr ? nullptr : &robust_riemann_);
        compute_gradient_face_operands_into(
            operand_workspace_.at(block.id()), block, metrics_.at(block.id()), profile_, version_);
    }
    inviscid_flux_plan_.set_version(version_);
    operand_plan_.set_version(version_);
    inviscid_flux_exchanger_.exchange(inviscid_flux_registry_);
    operand_exchanger_.exchange(operand_registry_);

    for (auto& block : local_blocks_.blocks()) {
        compute_primitive_gradients_into(gradient_workspace_.at(block.id()),
                                         block,
                                         metrics_.at(block.id()),
                                         operand_workspace_.at(block.id()),
                                         profile_);
    }
    gradient_plan_.set_version(version_);
    gradient_exchanger_.exchange(gradient_registry_);

    for (auto& block : local_blocks_.blocks()) {
        compute_viscous_face_fluxes_into(viscous_flux_workspace_.at(block.id()),
                                         block,
                                         metrics_.at(block.id()),
                                         gradient_workspace_.at(block.id()),
                                         profile_,
                                         transport_,
                                         boundary_data_.at(block.id()),
                                         gas_,
                                         reference_,
                                         floors_,
                                         version_);
    }
    viscous_flux_plan_.set_version(version_);
    viscous_flux_exchanger_.exchange(viscous_flux_registry_);

    for (auto& block : local_blocks_.blocks()) {
        compute_wcns_inviscid_residual(block,
                                       metrics_.at(block.id()),
                                       inviscid_flux_workspace_.at(block.id()),
                                       profile_,
                                       config_.inviscid.flux_difference);
        add_wcns_viscous_residual(block,
                                  metrics_.at(block.id()),
                                  viscous_flux_workspace_.at(block.id()),
                                  profile_,
                                  reference_.reynolds());
        add_source_terms(block, metrics_.at(block.id()), source_registry_, stage_time);
    }
}

Real ViscousWcnsSolver::advance(Real time_step, Real initial_time)
{
    Real accepted_time_step = time_step;
    if (!config_.inviscid.robustness.enabled) {
        robustness_diagnostics_ = {};
        robustness_diagnostics_.proposed_time_step = time_step;
        robustness_diagnostics_.accepted_time_step = time_step;
        int rk_stage = 0;
        advance_ssprk3(
            block_workspace_,
            time_workspace_,
            time_step,
            initial_time,
            [this, &rk_stage](Real stage_time) { compute_residuals(stage_time, ++rk_stage); });
    } else {
        accepted_time_step = advance_ssprk3_with_robustness(
            mpi_,
            block_workspace_,
            global_mesh_,
            profile_,
            config_.inviscid.flux_difference,
            gas_,
            reference_,
            floors_,
            config_.inviscid.robustness,
            robustness_ladder_,
            time_step,
            initial_time,
            [this](Real stage_time, int rk_stage, const BlockFaceRobustnessMap& levels) {
                compute_residuals_impl(stage_time, rk_stage, &levels);
            },
            robustness_diagnostics_);
    }
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
    }
    return accepted_time_step;
}

std::size_t ViscousWcnsSolver::global_reconstruction_fallback_count() const
{
    return global_diagnostic_count(mpi_,
                                   reconstruction_diagnostics_.fallback_events.size(),
                                   "viscous reconstruction fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.fallback_count(), "viscous Riemann fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_face_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.total_faces, "viscous Riemann face count");
}

RobustnessDiagnostics ViscousWcnsSolver::global_robustness_diagnostics() const
{
    RobustnessDiagnostics result = robustness_diagnostics_;
    for (std::size_t level = 0; level < result.face_levels.size(); ++level) {
        result.face_levels[level]
            = global_diagnostic_count(mpi_,
                                      robustness_diagnostics_.face_levels[level],
                                      "viscous robustness face-level count");
    }
    result.troubled_cells = global_diagnostic_count(
        mpi_, robustness_diagnostics_.troubled_cells, "viscous robustness troubled-cell count");
    result.local_recomputations = static_cast<std::size_t>(
        mpi_.max(static_cast<Real>(robustness_diagnostics_.local_recomputations)));
    result.step_retries = static_cast<std::size_t>(
        mpi_.max(static_cast<Real>(robustness_diagnostics_.step_retries)));
    result.minimum_density = mpi_.min(robustness_diagnostics_.minimum_density);
    result.minimum_pressure = mpi_.min(robustness_diagnostics_.minimum_pressure);
    result.minimum_temperature = mpi_.min(robustness_diagnostics_.minimum_temperature);
    result.minimum_internal_energy = mpi_.min(robustness_diagnostics_.minimum_internal_energy);
    return result;
}

Real ViscousWcnsSolver::global_time_step(Real cfl)
{
    if (!positive_finite(cfl)) {
        throw std::invalid_argument("viscous WCNS CFL must be positive and finite");
    }
    Real local_minimum = std::numeric_limits<Real>::infinity();
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
        const auto& metric = metrics_.at(block.id());
        const auto cells = block.cell_extent();
        const Real viscous_stability
            = config_.stability.for_ssprk3(profile_.kind(), block.cell_dimension());
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const Index3 cell {i, j, k};
                    TemperaturePrimitiveState state {};
                    for (int component = 0; component < fluid_components; ++component) {
                        state[static_cast<std::size_t>(component)]
                            = block.flow.temperature_primitive(i, j, k, component);
                    }
                    const std::array<Real, 3> velocity {{
                        state[temperature_velocity_x],
                        state[temperature_velocity_y],
                        state[temperature_velocity_z],
                    }};
                    const Real sound = thermodynamic_sound_speed(
                        state, gas_, reference_, floors_, block.cell_dimension());
                    Real inviscid_sum = 0.0;
                    Real area_square_sum = 0.0;
                    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        for (int side = 0; side <= 1; ++side) {
                            auto face = cell;
                            face[static_cast<std::size_t>(axis)] += side;
                            const auto area = area_vector(metric, axis, face);
                            const Real area2 = area_squared(area);
                            inviscid_sum += std::abs(velocity[0] * area[0] + velocity[1] * area[1]
                                                     + velocity[2] * area[2])
                                + sound * std::sqrt(area2);
                            area_square_sum += area2;
                        }
                    }
                    const Real jacobian = metric.jacobian()(i, j, k);
                    const Real mu = transport_.viscosity(state[temperature_value]);
                    const Real rho = state[temperature_density];
                    const Real nu_effective
                        = std::max(mu / rho, mu / (rho * config_.transport.prandtl));
                    const Real denominator = inviscid_sum / (2.0 * jacobian)
                        + viscous_stability * nu_effective * area_square_sum
                            / (2.0 * reference_.reynolds() * jacobian * jacobian);
                    if (!positive_finite(denominator)) {
                        throw PhysicsError("viscous WCNS time-step denominator is invalid");
                    }
                    local_minimum = std::min(local_minimum, cfl / denominator);
                }
            }
        }
    }
    const Real result = mpi_.min(local_minimum);
    if (!positive_finite(result)) {
        throw PhysicsError("viscous WCNS global time step is invalid");
    }
    return result;
}

Real ViscousWcnsSolver::global_residual_l2() const
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
    if (count <= 0.0) throw PhysicsError("viscous residual norm has no cells");
    return std::sqrt(mpi_.sum(local_sum) / count);
}

} // namespace wcns
