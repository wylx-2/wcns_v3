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
        && lhs.temperature == rhs.temperature
        && lhs.jacobian_absolute == rhs.jacobian_absolute
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

std::array<Real, 3> area_vector(
    const MetricField& metric, Axis axis, Index3 face)
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

} // namespace

void ViscousStabilityCoefficients::validate() const
{
    if (!positive_finite(phenglei_2d_ssprk3)
        || !positive_finite(phenglei_3d_ssprk3)
        || !positive_finite(scmm6_2d_ssprk3)
        || !positive_finite(scmm6_3d_ssprk3)) {
        throw std::invalid_argument(
            "all profile/dimension SSPRK3 viscous stability coefficients "
            "must be positive and finite");
    }
}

Real ViscousStabilityCoefficients::for_ssprk3(
    AlgorithmProfileKind profile,
    int dimension) const
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
    return inviscid.summary() + ';' + transport.summary()
        + "; " + stability.summary();
}

std::string ViscousWcnsConfig::restart_signature() const
{
    return "viscous_wcns_v1;" + summary();
}

ViscousWcnsSolver::ViscousWcnsSolver(
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
    ViscousWcnsConfig config)
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
    , source_registry_(SourceTermRegistry::create_stage_j(
          config_.inviscid.source_terms))
    , transport_(config_.transport)
{
    config_.validate();
    const auto riemann_registry = RiemannSolverRegistry::with_builtins(
        config_.inviscid.riemann.parameters);
    config_.inviscid.riemann.validate(riemann_registry);
    riemann_ = RiemannSolver(
        config_.inviscid.riemann.scheme, riemann_registry,
        config_.inviscid.riemann.parameters);
    floors_.validate();
    if (!same_floors(config_.inviscid.reconstruction.floors, floors_)) {
        throw std::invalid_argument(
            "viscous WCNS reconstruction and solver floors differ");
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
    }
}

void ViscousWcnsSolver::compute_residuals(Real stage_time, int rk_stage)
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

    std::unordered_map<BlockId, InviscidFaceFluxField> inviscid_fluxes;
    FaceFluxFieldRegistry inviscid_registry;
    std::unordered_map<BlockId, GradientOperandFaceField> operands;
    GradientOperandFieldRegistry operand_registry;
    reconstruction_diagnostics_ = {};
    riemann_diagnostics_ = {};
    for (auto& block : local_blocks_.blocks()) {
        const auto& data = boundary_data_.at(block.id());
        const auto ghost = PhysicalGhostStateOperator::fill(
            block, data, gas_, reference_, floors_, version_);
        if (ghost.version != version_) {
            throw std::logic_error("viscous WCNS physical ghost version mismatch");
        }
        auto [inviscid, inserted_inviscid] = inviscid_fluxes.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(compute_inviscid_face_fluxes(
                block, metrics_.at(block.id()), profile_,
                config_.inviscid.reconstruction, riemann_, gas_, reference_, floors_,
                data, config_.inviscid.boundary, version_,
                reconstruction_diagnostics_, &riemann_diagnostics_, rk_stage)));
        if (!inserted_inviscid) throw std::logic_error("duplicate inviscid flux block");
        inviscid_registry.add(block.id(), inviscid->second);

        auto [operand, inserted_operand] = operands.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(compute_gradient_face_operands(
                block, metrics_.at(block.id()), profile_, version_)));
        if (!inserted_operand) throw std::logic_error("duplicate gradient operand block");
        operand_registry.add(block.id(), operand->second);
    }
    FaceFluxHaloExchanger(
        mpi_, FaceFluxHaloPlan::build(global_mesh_, profile_, version_))
        .exchange(inviscid_registry);
    GradientOperandFaceHaloExchanger(
        mpi_, GradientOperandFaceHaloPlan::build(
                  global_mesh_, profile_, version_))
        .exchange(operand_registry);

    std::unordered_map<BlockId, PrimitiveGradientField> gradients;
    GradientFieldRegistry gradient_registry;
    for (auto& block : local_blocks_.blocks()) {
        auto [gradient, inserted] = gradients.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(compute_primitive_gradients(
                block, metrics_.at(block.id()), operands.at(block.id()), profile_)));
        if (!inserted) throw std::logic_error("duplicate primitive gradient block");
        gradient_registry.add(block.id(), gradient->second);
    }
    GradientHaloExchanger(
        mpi_, GradientHaloPlan::build(
                  global_mesh_, topology_, profile_, version_))
        .exchange(gradient_registry);

    std::unordered_map<BlockId, ViscousFaceFluxField> viscous_fluxes;
    ViscousFaceFluxFieldRegistry viscous_registry;
    for (auto& block : local_blocks_.blocks()) {
        auto [flux, inserted] = viscous_fluxes.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(compute_viscous_face_fluxes(
                block, metrics_.at(block.id()), gradients.at(block.id()),
                profile_, transport_, boundary_data_.at(block.id()),
                gas_, reference_, floors_, version_)));
        if (!inserted) throw std::logic_error("duplicate viscous flux block");
        viscous_registry.add(block.id(), flux->second);
    }
    ViscousFaceFluxHaloExchanger(
        mpi_, ViscousFaceFluxHaloPlan::build(
                  global_mesh_, profile_, version_))
        .exchange(viscous_registry);

    for (auto& block : local_blocks_.blocks()) {
        compute_wcns_inviscid_residual(
            block, metrics_.at(block.id()), inviscid_fluxes.at(block.id()), profile_);
        add_wcns_viscous_residual(
            block, metrics_.at(block.id()), viscous_fluxes.at(block.id()),
            profile_, reference_.reynolds());
        add_source_terms(
            block, metrics_.at(block.id()), source_registry_, stage_time);
    }
}

void ViscousWcnsSolver::advance(Real time_step, Real initial_time)
{
    std::vector<StructuredBlock*> blocks;
    for (auto& block : local_blocks_.blocks()) blocks.push_back(&block);
    int rk_stage = 0;
    advance_ssprk3(blocks, time_step, initial_time,
        [this, &rk_stage](Real stage_time) {
            compute_residuals(stage_time, ++rk_stage);
        });
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
    }
}

std::size_t ViscousWcnsSolver::global_reconstruction_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, reconstruction_diagnostics_.fallback_events.size(),
        "viscous reconstruction fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.fallback_count(),
        "viscous Riemann fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_face_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.total_faces, "viscous Riemann face count");
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
        const Real viscous_stability = config_.stability.for_ssprk3(
            profile_.kind(), block.cell_dimension());
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
                            inviscid_sum += std::abs(
                                velocity[0] * area[0] + velocity[1] * area[1]
                                + velocity[2] * area[2])
                                + sound * std::sqrt(area2);
                            area_square_sum += area2;
                        }
                    }
                    const Real jacobian = metric.jacobian()(i, j, k);
                    const Real mu = transport_.viscosity(state[temperature_value]);
                    const Real rho = state[temperature_density];
                    const Real nu_effective = std::max(
                        mu / rho, mu / (rho * config_.transport.prandtl));
                    const Real denominator = inviscid_sum / (2.0 * jacobian)
                        + viscous_stability * nu_effective
                            * area_square_sum
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
