#include <wcns/solver/viscous_operator.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace wcns {
namespace {

constexpr std::uint64_t maximum_exact_message_version = 9007199254740992ULL;
constexpr int viscous_flux_tag_base = 24576;

int side_sign(Side side)
{
    return side == Side::Lower ? -1 : 1;
}

bool contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) return false;
    }
    return true;
}

const BoundaryPatch* physical_patch(
    const StructuredBlock& block, Axis axis, Index3 face)
{
    for (const auto& patch : block.boundaries) {
        if (patch.face.axis == axis
            && contains(patch.boundary_face_range.untyped(), face)) {
            return &patch;
        }
    }
    return nullptr;
}

bool connection_covers(
    const StructuredBlock& block, Axis axis, Side side, Index3 face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis
            && connection.receiver_face.side == side
            && contains(connection.shared_face_range.untyped(), face)) {
            return true;
        }
    }
    return false;
}

const FaceAreaVectors& face_metrics(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

Normal3 positive_unit_normal(const FaceAreaVectors& faces, Index3 face)
{
    const Real x = faces.x(face.i, face.j, face.k);
    const Real y = faces.y(face.i, face.j, face.k);
    const Real z = faces.z(face.i, face.j, face.k);
    const Real area = std::sqrt(x * x + y * y + z * z);
    if (!std::isfinite(area) || area <= 0.0) {
        throw PhysicsError("viscous face has invalid area");
    }
    return {x / area, y / area, z / area};
}

Normal3 outward(Normal3 positive, Side side)
{
    const Real sign = static_cast<Real>(side_sign(side));
    return {sign * positive.x, sign * positive.y, sign * positive.z};
}

Real centered_derivative(
    const Field<Real>& field,
    Axis axis,
    Index3 cell,
    int component,
    AlgorithmProfileKind profile)
{
    const auto value = [&](int offset) {
        auto face = cell;
        face[static_cast<std::size_t>(axis)] += offset;
        const Real result = field(face.i, face.j, face.k, component);
        if (!std::isfinite(result)) {
            throw PhysicsError("viscous face-flux halo is incomplete");
        }
        return result;
    };
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return (value(-1) - 27.0 * value(0) + 27.0 * value(1) - value(2))
            / 24.0;
    }
    return (-9.0 * value(-2) + 125.0 * value(-1) - 2250.0 * value(0)
        + 2250.0 * value(1) - 125.0 * value(2) + 9.0 * value(3))
        / 1920.0;
}

ConservativeState load_flux(
    const ViscousFaceFluxField& field, Axis axis, Index3 index)
{
    ConservativeState result {};
    const auto& values = field.field(axis);
    for (int component = 0; component < euler_components; ++component) {
        result[static_cast<std::size_t>(component)]
            = values(index.i, index.j, index.k, component);
    }
    return result;
}

void store_flux(
    ViscousFaceFluxField& field, Axis axis, Index3 index,
    const ConservativeState& state)
{
    auto& values = field.field(axis);
    for (int component = 0; component < euler_components; ++component) {
        values(index.i, index.j, index.k, component)
            = state[static_cast<std::size_t>(component)];
    }
}

void validate_field(
    const ViscousFaceFluxField& field,
    const FaceFluxExchangeDescriptor& descriptor)
{
    if (field.profile() != descriptor.profile || field.version() != descriptor.version) {
        throw std::invalid_argument("viscous face-flux field metadata mismatch");
    }
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("viscous face-flux message exceeds MPI int count");
    }
    return static_cast<int>(count);
}
#endif

} // namespace

ViscousFaceFluxField::ViscousFaceFluxField(
    Extent3 cells,
    int dimension,
    AlgorithmProfileKind profile,
    std::uint64_t version)
    : profile_(profile)
    , version_(version)
    , halo_layers_(profile == AlgorithmProfileKind::PhengleiWcns ? 1 : 2)
    , dimension_(dimension)
    , i_({cells.ni + 1, cells.nj, cells.nk}, euler_components, halo_layers_)
    , j_({cells.ni, cells.nj + 1, cells.nk}, euler_components, halo_layers_)
    , k_({cells.ni, cells.nj, cells.nk + 1}, euler_components, halo_layers_)
{
    if ((dimension != 2 && dimension != 3) || version == 0
        || version > maximum_exact_message_version) {
        throw std::invalid_argument("viscous face-flux metadata is invalid");
    }
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    i_.fill(nan);
    j_.fill(nan);
    k_.fill(nan);
}

Field<Real>& ViscousFaceFluxField::field(Axis axis)
{
    if (axis == Axis::I) return i_;
    if (axis == Axis::J) return j_;
    if (dimension_ != 3) throw std::out_of_range("2D viscous flux has no K field");
    return k_;
}

const Field<Real>& ViscousFaceFluxField::field(Axis axis) const
{
    return const_cast<ViscousFaceFluxField*>(this)->field(axis);
}

ConservativeState transform_viscous_face_flux_for_receiver(
    const ConservativeState& donor,
    const FaceFluxExchangeDescriptor& descriptor)
{
    ConservativeState result = donor;
    const auto rotated = descriptor.periodic.inverse().apply_vector({{
        donor[momentum_x], donor[momentum_y], donor[momentum_z],
    }});
    result[density] *= descriptor.orientation;
    result[momentum_x] = descriptor.orientation * rotated[0];
    result[momentum_y] = descriptor.orientation * rotated[1];
    result[momentum_z] = descriptor.orientation * rotated[2];
    result[total_energy] *= descriptor.orientation;
    return result;
}

ViscousFaceFluxHaloPlan ViscousFaceFluxHaloPlan::build(
    const StructuredMesh& mesh,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    ViscousFaceFluxHaloPlan result;
    const auto base = FaceFluxHaloPlan::build(mesh, profile, version);
    result.exchanges_ = base.exchanges();
    return result;
}

void ViscousFaceFluxFieldRegistry::add(
    BlockId block, ViscousFaceFluxField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("viscous face-flux registry has an invalid block");
    }
}

bool ViscousFaceFluxFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

ViscousFaceFluxField& ViscousFaceFluxFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("viscous face-flux field is not registered");
    }
    return *iterator->second;
}

void ViscousFaceFluxHaloExchanger::exchange(
    const ViscousFaceFluxFieldRegistry& fields) const
{
    struct Pending {
        const FaceFluxExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };
    const RankId rank = mpi_.rank();
    std::vector<Pending> receives;
    std::vector<Pending> sends;
    for (const auto& descriptor : plan_.exchanges()) {
        const std::size_t count = 1 + descriptor.pairs.size()
            * static_cast<std::size_t>(euler_components);
        if (descriptor.receiver_rank == rank && descriptor.donor_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            const auto& donor = fields.field(descriptor.donor_block);
            validate_field(receiver, descriptor);
            validate_field(donor, descriptor);
            for (const auto& pair : descriptor.pairs) {
                store_flux(receiver, descriptor.receiver_axis, pair.receiver,
                    transform_viscous_face_flux_for_receiver(
                        load_flux(donor, descriptor.donor_axis, pair.donor),
                        descriptor));
            }
        } else if (descriptor.receiver_rank == rank) {
            validate_field(fields.field(descriptor.receiver_block), descriptor);
            receives.push_back({&descriptor, std::vector<Real>(count)});
        } else if (descriptor.donor_rank == rank) {
            const auto& donor = fields.field(descriptor.donor_block);
            validate_field(donor, descriptor);
            Pending pending {&descriptor, std::vector<Real>(count)};
            pending.values[0] = static_cast<Real>(descriptor.version);
            std::size_t offset = 1;
            for (const auto& pair : descriptor.pairs) {
                const auto state = load_flux(
                    donor, descriptor.donor_axis, pair.donor);
                for (const Real value : state) pending.values[offset++] = value;
            }
            sends.push_back(std::move(pending));
        }
    }
#if WCNS_HAS_MPI
    std::vector<MPI_Request> requests(receives.size() + sends.size(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives) {
        check_mpi(MPI_Irecv(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->donor_rank,
            pending.descriptor->message_tag(viscous_flux_tag_base),
            mpi_.communicator(), &requests[request++]),
            "MPI_Irecv viscous face flux");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->receiver_rank,
            pending.descriptor->message_tag(viscous_flux_tag_base),
            mpi_.communicator(), &requests[request++]),
            "MPI_Isend viscous face flux");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall viscous face flux");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote viscous face-flux exchange requires MPI");
    }
#endif
    for (const auto& pending : receives) {
        if (pending.values[0] != static_cast<Real>(pending.descriptor->version)) {
            throw MpiError("viscous face-flux message version mismatch");
        }
        auto& receiver = fields.field(pending.descriptor->receiver_block);
        std::size_t offset = 1;
        for (const auto& pair : pending.descriptor->pairs) {
            ConservativeState donor {};
            for (auto& value : donor) value = pending.values[offset++];
            store_flux(receiver, pending.descriptor->receiver_axis, pair.receiver,
                transform_viscous_face_flux_for_receiver(
                    donor, *pending.descriptor));
        }
    }
}

ViscousFaceFluxField compute_viscous_face_fluxes(
    const StructuredBlock& block,
    const MetricField& metric,
    const PrimitiveGradientField& gradients,
    const AlgorithmProfile& profile,
    const TransportModel& transport,
    const BoundaryDataMap& boundary_data,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    std::uint64_t version)
{
    if (metric.profile() != profile.kind() || gradients.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()
        || gradients.dimension() != block.cell_dimension()
        || gradients.version() != version) {
        throw ProfileError("viscous face flux inputs have incompatible metadata");
    }
    ViscousFaceFluxField result(
        block.cell_extent(), block.cell_dimension(), profile.kind(), version);
    const auto compute_axis = [&](Axis axis) {
        const auto& faces = face_metrics(metric, axis);
        const auto extent = faces.x.interior_extent();
        auto& output = result.field(axis);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 face {i, j, k};
                    auto trace = interpolate_viscous_face_trace(
                        block, gradients, profile, axis, face);
                    if (const auto* patch = physical_patch(block, axis, face)) {
                        const auto data = boundary_data.find(patch->name);
                        if (data == boundary_data.end()) {
                            throw PhysicsConfigurationError(
                                "viscous boundary data is missing for patch " + patch->name);
                        }
                        const auto normal = outward(
                            positive_unit_normal(faces, face), patch->face.side);
                        trace = apply_viscous_boundary_trace(
                            block, metric, *patch, face, trace,
                            interpolate_internal_pressure_trace(
                                block, profile, axis, face),
                            normal, data->second, profile, gas, reference, floors);
                    }
                    const auto cartesian = compute_viscous_cartesian_flux(
                        trace, transport, gas, reference, floors,
                        block.cell_dimension());
                    const Real sx = faces.x(i, j, k);
                    const Real sy = faces.y(i, j, k);
                    const Real sz = faces.z(i, j, k);
                    for (int component = 0; component < euler_components; ++component) {
                        const Real value = sx * cartesian.x[static_cast<std::size_t>(component)]
                            + sy * cartesian.y[static_cast<std::size_t>(component)]
                            + sz * cartesian.z[static_cast<std::size_t>(component)];
                        if (!std::isfinite(value)) {
                            throw PhysicsError("viscous computational flux is non-finite");
                        }
                        output(i, j, k, component) = value;
                    }
                }
            }
        }
    };
    compute_axis(Axis::I);
    compute_axis(Axis::J);
    if (block.cell_dimension() == 3) compute_axis(Axis::K);
    return result;
}

void add_wcns_viscous_residual(
    StructuredBlock& block,
    const MetricField& metric,
    const ViscousFaceFluxField& flux,
    const AlgorithmProfile& profile,
    Real reynolds)
{
    if (metric.profile() != profile.kind() || flux.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()
        || flux.dimension() != block.cell_dimension()
        || !std::isfinite(reynolds) || reynolds <= 0.0) {
        throw ProfileError("viscous residual inputs are incompatible or Re is invalid");
    }
    const auto cells = block.cell_extent();
    const auto accumulate_axis = [&](Axis axis) {
        const int count = cells[static_cast<std::size_t>(axis)];
        const auto operators = LineOperators::build(profile, count);
        const auto& values = flux.field(axis);
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const Index3 cell {i, j, k};
                    const int normal = cell[static_cast<std::size_t>(axis)];
                    Index3 lower = cell;
                    lower[static_cast<std::size_t>(axis)] = 0;
                    Index3 upper = cell;
                    upper[static_cast<std::size_t>(axis)] = count;
                    const bool connected_lower
                        = connection_covers(block, axis, Side::Lower, lower);
                    const bool connected_upper
                        = connection_covers(block, axis, Side::Upper, upper);
                    const int width = profile.kind()
                            == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
                    const Real jacobian = metric.jacobian()(i, j, k);
                    if (!std::isfinite(jacobian) || jacobian <= 0.0) {
                        throw PhysicsError("viscous residual has an invalid Jacobian");
                    }
                    for (int component = 0; component < euler_components; ++component) {
                        Real derivative = 0.0;
                        if ((connected_lower && normal < width)
                            || (connected_upper && normal >= count - width)) {
                            derivative = centered_derivative(
                                values, axis, cell, component, profile.kind());
                        } else {
                            const auto& row = operators.derivative_rows()[
                                static_cast<std::size_t>(normal)];
                            for (const auto [face_index, coefficient] : row) {
                                auto face = cell;
                                face[static_cast<std::size_t>(axis)] = face_index;
                                derivative += coefficient * values(
                                    face.i, face.j, face.k, component);
                            }
                        }
                        block.flow.residual(i, j, k, component)
                            += derivative / (reynolds * jacobian);
                    }
                }
            }
        }
    };
    accumulate_axis(Axis::I);
    accumulate_axis(Axis::J);
    if (block.cell_dimension() == 3) accumulate_axis(Axis::K);
}

} // namespace wcns
