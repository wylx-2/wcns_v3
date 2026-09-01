#include <wcns/solver/inviscid_flux.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>

namespace wcns {
namespace {

constexpr std::uint64_t maximum_exact_message_version = 9007199254740992ULL;

int side_sign(Side side)
{
    return side == Side::Lower ? -1 : 1;
}

int inward_sign(Side side)
{
    return side == Side::Lower ? 1 : -1;
}

bool contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) {
            return false;
        }
    }
    return true;
}

const ConnectivityPatch& reciprocal(
    const StructuredMesh& mesh, const ConnectivityPatch& connection)
{
    const auto& donor = mesh.block(connection.donor_block);
    const auto iterator = std::find_if(
        donor.connectivities.begin(), donor.connectivities.end(),
        [&](const ConnectivityPatch& candidate) {
            return candidate.receiver_block == connection.donor_block
                && candidate.donor_block == connection.receiver_block
                && candidate.receiver_face == connection.donor_face
                && candidate.donor_face == connection.receiver_face
                && candidate.transform
                    == connection.transform.inverse(donor.cell_dimension());
        });
    if (iterator == donor.connectivities.end()) {
        throw TopologyError("face-flux plan cannot find reciprocal connectivity");
    }
    return *iterator;
}

auto connection_key(const ConnectivityPatch& connection)
{
    return std::tuple {
        connection.receiver_block,
        connection.donor_block,
        static_cast<int>(connection.receiver_face.axis),
        static_cast<int>(connection.receiver_face.side),
        connection.receiver_vertex_range.begin.i,
        connection.receiver_vertex_range.begin.j,
        connection.receiver_vertex_range.begin.k,
    };
}

FaceFluxExchangeDescriptor make_descriptor(
    const StructuredMesh& mesh,
    const ConnectivityPatch& connection,
    ConnectionId id,
    BlockId owner,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    FaceFluxExchangeDescriptor descriptor;
    descriptor.connection = id;
    descriptor.receiver_block = connection.receiver_block;
    descriptor.donor_block = connection.donor_block;
    descriptor.receiver_rank = mesh.block(connection.receiver_block).owner_rank();
    descriptor.donor_rank = mesh.block(connection.donor_block).owner_rank();
    descriptor.shared_face_owner = owner;
    descriptor.receiver_axis = connection.receiver_face.axis;
    descriptor.donor_axis = connection.donor_face.axis;
    descriptor.orientation = static_cast<Real>(
        -side_sign(connection.receiver_face.side)
        / side_sign(connection.donor_face.side));
    descriptor.periodic = connection.periodic;
    descriptor.profile = profile.kind();
    descriptor.version = version;
    const int maximum_layer
        = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
    const int first_layer = connection.receiver_block == owner ? 1 : 0;
    const auto& reverse = reciprocal(mesh, connection);
    const auto counts = connection.shared_face_range.counts();
    for (int layer = first_layer; layer <= maximum_layer; ++layer) {
        for (int k = 0; k < counts.nk; ++k) {
            for (int j = 0; j < counts.nj; ++j) {
                for (int i = 0; i < counts.ni; ++i) {
                    const Index3 receiver_ordinal {i, j, k};
                    Index3 donor_ordinal;
                    for (int receiver_axis = 0;
                         receiver_axis < mesh.block(connection.receiver_block).cell_dimension();
                         ++receiver_axis) {
                        const int donor_axis = std::abs(
                            connection.transform.receiver_to_donor[
                                static_cast<std::size_t>(receiver_axis)]) - 1;
                        donor_ordinal[static_cast<std::size_t>(donor_axis)]
                            = receiver_ordinal[static_cast<std::size_t>(receiver_axis)];
                    }
                    auto receiver_index
                        = connection.shared_face_range.at(receiver_ordinal);
                    auto donor_index = reverse.shared_face_range.at(donor_ordinal);
                    receiver_index[static_cast<std::size_t>(connection.receiver_face.axis)]
                        += side_sign(connection.receiver_face.side) * layer;
                    donor_index[static_cast<std::size_t>(connection.donor_face.axis)]
                        += inward_sign(connection.donor_face.side) * layer;
                    descriptor.pairs.push_back(
                        {receiver_index, donor_index, layer});
                }
            }
        }
    }
    return descriptor;
}

ConservativeState transform_flux_impl(
    const ConservativeState& donor,
    const FaceFluxExchangeDescriptor& descriptor)
{
    ConservativeState result = donor;
    const std::array<Real, 3> momentum {{donor[1], donor[2], donor[3]}};
    const auto rotated = descriptor.periodic.inverse().apply_vector(momentum);
    result[0] *= descriptor.orientation;
    result[1] = descriptor.orientation * rotated[0];
    result[2] = descriptor.orientation * rotated[1];
    result[3] = descriptor.orientation * rotated[2];
    result[4] *= descriptor.orientation;
    return result;
}

ConservativeState load_flux(
    const InviscidFaceFluxField& field, Axis axis, Index3 index)
{
    const auto& values = field.field(axis);
    ConservativeState result {};
    for (int component = 0; component < euler_components; ++component) {
        result[static_cast<std::size_t>(component)]
            = values(index.i, index.j, index.k, component);
    }
    return result;
}

void store_flux(
    InviscidFaceFluxField& field, Axis axis, Index3 index,
    const ConservativeState& state)
{
    auto& values = field.field(axis);
    for (int component = 0; component < euler_components; ++component) {
        values(index.i, index.j, index.k, component)
            = state[static_cast<std::size_t>(component)];
    }
}

void validate_field(
    const InviscidFaceFluxField& field,
    const FaceFluxExchangeDescriptor& descriptor)
{
    if (field.profile() != descriptor.profile || field.version() != descriptor.version) {
        throw std::invalid_argument("face-flux field profile or version mismatch");
    }
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

const FaceAreaVectors& metric_faces(const MetricField& metric, Axis axis)
{
    switch (axis) {
    case Axis::I: return metric.i_faces();
    case Axis::J: return metric.j_faces();
    case Axis::K: return metric.k_faces();
    }
    throw std::invalid_argument("invalid metric face axis");
}

Normal3 unit_normal(const FaceAreaVectors& metric, Index3 face, Real& area)
{
    area = metric.area(face.i, face.j, face.k);
    if (!std::isfinite(area) || area <= 0.0) {
        throw PhysicsError("inviscid face has invalid area");
    }
    return {
        metric.x(face.i, face.j, face.k) / area,
        metric.y(face.i, face.j, face.k) / area,
        metric.z(face.i, face.j, face.k) / area,
    };
}

Normal3 outward(Normal3 positive, Side side)
{
    const Real sign = side == Side::Lower ? -1.0 : 1.0;
    return {sign * positive.x, sign * positive.y, sign * positive.z};
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

Real centered_derivative(
    const Field<Real>& flux,
    Axis axis,
    Index3 cell,
    int component,
    AlgorithmProfileKind profile)
{
    const auto value = [&](int face_offset) {
        auto index = cell;
        index[static_cast<std::size_t>(axis)] += face_offset;
        const Real result = flux(index.i, index.j, index.k, component);
        if (!std::isfinite(result)) {
            throw PhysicsError("face-flux halo contains a non-finite value");
        }
        return result;
    };
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return (value(-1) - 27.0 * value(0) + 27.0 * value(1) - value(2)) / 24.0;
    }
    return (-9.0 * value(-2) + 125.0 * value(-1) - 2250.0 * value(0)
        + 2250.0 * value(1) - 125.0 * value(2) + 9.0 * value(3)) / 1920.0;
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("face-flux MPI message exceeds int count");
    }
    return static_cast<int>(count);
}
#endif

} // namespace

ConservativeState transform_inviscid_face_flux_for_receiver(
    const ConservativeState& donor,
    const FaceFluxExchangeDescriptor& descriptor)
{
    return transform_flux_impl(donor, descriptor);
}

InviscidFaceFluxField::InviscidFaceFluxField(
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
        throw std::invalid_argument("face-flux field has invalid dimension or version");
    }
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    i_.fill(nan);
    j_.fill(nan);
    k_.fill(nan);
}

Field<Real>& InviscidFaceFluxField::field(Axis axis)
{
    switch (axis) {
    case Axis::I: return i_;
    case Axis::J: return j_;
    case Axis::K:
        if (dimension_ != 3) throw std::out_of_range("2D face flux has no K field");
        return k_;
    }
    throw std::invalid_argument("invalid face-flux axis");
}

const Field<Real>& InviscidFaceFluxField::field(Axis axis) const
{
    return const_cast<InviscidFaceFluxField*>(this)->field(axis);
}

int FaceFluxExchangeDescriptor::message_tag(int tag_base) const
{
    if (tag_base < 0 || connection < 0 || receiver_block < 0 || donor_block < 0) {
        throw TopologyError("face-flux message tag inputs are invalid");
    }
    const int direction = receiver_block < donor_block ? 0 : 1;
    const long long tag = static_cast<long long>(tag_base)
        + 8LL * connection + 2LL * static_cast<int>(profile) + direction;
    if (tag > std::numeric_limits<int>::max()) {
        throw TopologyError("face-flux message tag exceeds int range");
    }
    return static_cast<int>(tag);
}

FaceFluxHaloPlan FaceFluxHaloPlan::build(
    const StructuredMesh& mesh,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    if (version == 0 || version > maximum_exact_message_version) {
        throw TopologyError("face-flux plan version must be non-zero");
    }
    mesh.validate_connectivities();
    std::vector<const ConnectivityPatch*> canonical;
    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            if (connection.receiver_block < connection.donor_block) {
                canonical.push_back(&connection);
            }
        }
    }
    std::sort(canonical.begin(), canonical.end(), [](const auto* lhs, const auto* rhs) {
        return connection_key(*lhs) < connection_key(*rhs);
    });
    FaceFluxHaloPlan result;
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto id = static_cast<ConnectionId>(index);
        const auto& forward = *canonical[index];
        const auto& reverse = reciprocal(mesh, forward);
        const BlockId owner = std::min(forward.receiver_block, forward.donor_block);
        result.exchanges_.push_back(
            make_descriptor(mesh, forward, id, owner, profile, version));
        result.exchanges_.push_back(
            make_descriptor(mesh, reverse, id, owner, profile, version));
    }
    std::set<int> tags;
    for (const auto& descriptor : result.exchanges_) {
        if (!tags.insert(descriptor.message_tag()).second) {
            throw TopologyError("face-flux plan generated duplicate message tags");
        }
    }
    return result;
}

void FaceFluxFieldRegistry::add(BlockId block, InviscidFaceFluxField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("face-flux registry contains invalid or duplicate block");
    }
}

bool FaceFluxFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

InviscidFaceFluxField& FaceFluxFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("face-flux field is not registered");
    }
    return *iterator->second;
}

void FaceFluxHaloExchanger::exchange(const FaceFluxFieldRegistry& fields) const
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
                    transform_inviscid_face_flux_for_receiver(
                        load_flux(donor, descriptor.donor_axis, pair.donor), descriptor));
            }
        } else if (descriptor.receiver_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            validate_field(receiver, descriptor);
            receives.push_back({&descriptor, std::vector<Real>(count)});
        } else if (descriptor.donor_rank == rank) {
            const auto& donor = fields.field(descriptor.donor_block);
            validate_field(donor, descriptor);
            Pending pending {&descriptor, std::vector<Real>(count)};
            pending.values[0] = static_cast<Real>(descriptor.version);
            std::size_t offset = 1;
            for (const auto& pair : descriptor.pairs) {
                const auto value = load_flux(donor, descriptor.donor_axis, pair.donor);
                for (const auto component : value) pending.values[offset++] = component;
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
            pending.descriptor->donor_rank, pending.descriptor->message_tag(),
            mpi_.communicator(), &requests[request++]), "MPI_Irecv face flux");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->receiver_rank, pending.descriptor->message_tag(),
            mpi_.communicator(), &requests[request++]), "MPI_Isend face flux");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall face flux");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote face-flux exchange requires WCNS_ENABLE_MPI");
    }
#endif
    for (const auto& pending : receives) {
        auto& receiver = fields.field(pending.descriptor->receiver_block);
        if (pending.values.empty()
            || pending.values[0] != static_cast<Real>(pending.descriptor->version)) {
            throw MpiError("face-flux message version mismatch");
        }
        std::size_t offset = 1;
        for (const auto& pair : pending.descriptor->pairs) {
            ConservativeState donor {};
            for (auto& component : donor) component = pending.values[offset++];
            store_flux(receiver, pending.descriptor->receiver_axis, pair.receiver,
                transform_inviscid_face_flux_for_receiver(
                    donor, *pending.descriptor));
        }
    }
}

InviscidFaceFluxField compute_inviscid_face_fluxes(
    const StructuredBlock& block,
    const MetricField& metric,
    const AlgorithmProfile& profile,
    const ReconstructionConfig& reconstruction,
    const RiemannSolver& riemann,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const BoundaryDataMap& boundary_data,
    const InviscidBoundaryOptions& boundary_options,
    std::uint64_t version,
    ReconstructionDiagnostics& diagnostics)
{
    ProfileFactory::validate_bundle(profile.components());
    if (metric.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()) {
        throw ProfileError("inviscid flux metric belongs to another profile or dimension");
    }
    InviscidFaceFluxField result(
        block.cell_extent(), block.cell_dimension(), profile.kind(), version);
    const auto compute_axis = [&](Axis axis) {
        const auto& faces = metric_faces(metric, axis);
        const auto extent = faces.x.interior_extent();
        auto& output = result.field(axis);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 face {i, j, k};
                    auto states = reconstruct_thermodynamic_face(
                        block.flow.conservative, block.flow.primitive,
                        axis, face, reconstruction, gas, reference,
                        diagnostics, block.cell_dimension());
                    Real area = 0.0;
                    const auto normal = unit_normal(faces, face, area);
                    if (const auto* patch = physical_patch(block, axis, face)) {
                        const auto data_iterator = boundary_data.find(patch->name);
                        if (data_iterator == boundary_data.end()) {
                            throw PhysicsConfigurationError(
                                "inviscid face boundary data is missing for patch " + patch->name);
                        }
                        if (patch->face.side == Side::Lower) {
                            states.left = apply_inviscid_boundary_face_state(
                                *patch, states.right, states.left,
                                outward(normal, Side::Lower), data_iterator->second,
                                boundary_options, gas, reference, floors,
                                block.cell_dimension());
                        } else {
                            states.right = apply_inviscid_boundary_face_state(
                                *patch, states.left, states.right,
                                outward(normal, Side::Upper), data_iterator->second,
                                boundary_options, gas, reference, floors,
                                block.cell_dimension());
                        }
                    }
                    const auto numerical
                        = riemann.flux(states.left, states.right, normal, gas, floors);
                    for (int component = 0; component < euler_components; ++component) {
                        output(i, j, k, component)
                            = area * numerical[static_cast<std::size_t>(component)];
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

void compute_wcns_inviscid_residual(
    StructuredBlock& block,
    const MetricField& metric,
    const InviscidFaceFluxField& flux,
    const AlgorithmProfile& profile)
{
    if (metric.profile() != profile.kind() || flux.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()
        || flux.dimension() != block.cell_dimension()) {
        throw ProfileError("WCNS residual inputs belong to different profiles or dimensions");
    }
    block.flow.residual.fill(0.0);
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
                    Index3 lower_face = cell;
                    lower_face[static_cast<std::size_t>(axis)] = 0;
                    Index3 upper_face = cell;
                    upper_face[static_cast<std::size_t>(axis)] = count;
                    const bool lower_connection
                        = connection_covers(block, axis, Side::Lower, lower_face);
                    const bool upper_connection
                        = connection_covers(block, axis, Side::Upper, upper_face);
                    for (int component = 0; component < euler_components; ++component) {
                        Real derivative = 0.0;
                        const int boundary_width
                            = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
                        if ((lower_connection && normal < boundary_width)
                            || (upper_connection && normal >= count - boundary_width)) {
                            derivative = centered_derivative(
                                values, axis, cell, component, profile.kind());
                        } else {
                            const auto& row = operators.derivative_rows()[
                                static_cast<std::size_t>(normal)];
                            for (const auto [face_index, coefficient] : row) {
                                auto face = cell;
                                face[static_cast<std::size_t>(axis)] = face_index;
                                derivative += coefficient
                                    * values(face.i, face.j, face.k, component);
                            }
                        }
                        const Real jacobian = metric.jacobian()(i, j, k);
                        if (!std::isfinite(derivative) || !std::isfinite(jacobian)
                            || jacobian <= 0.0) {
                            throw PhysicsError("WCNS flux divergence is non-finite or has invalid Jacobian");
                        }
                        block.flow.residual(i, j, k, component)
                            -= derivative / jacobian;
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
