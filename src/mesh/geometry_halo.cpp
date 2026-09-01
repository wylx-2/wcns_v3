#include <wcns/mesh/geometry_halo.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <tuple>

namespace wcns {

struct GeometryHaloPlanBuilderAccess {
    static std::vector<GeometryExchangeDescriptor>& exchanges(
        GeometryHaloPlan& plan)
    {
        return plan.exchanges_;
    }
};

namespace {

Side opposite(Side side)
{
    return side == Side::Lower ? Side::Upper : Side::Lower;
}

int side_sign(Side side)
{
    return side == Side::Lower ? -1 : 1;
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
        connection.receiver_vertex_range.end.i,
        connection.receiver_vertex_range.end.j,
        connection.receiver_vertex_range.end.k,
    };
}

const ConnectivityPatch& reciprocal_connection(
    const StructuredMesh& mesh,
    const ConnectivityPatch& connection)
{
    const auto& donor = mesh.block(connection.donor_block);
    const auto iterator = std::find_if(
        donor.connectivities.begin(),
        donor.connectivities.end(),
        [&](const ConnectivityPatch& candidate) {
            return candidate.receiver_block == connection.donor_block
                && candidate.donor_block == connection.receiver_block
                && candidate.receiver_face == connection.donor_face
                && candidate.donor_face == connection.receiver_face
                && candidate.transform
                    == connection.transform.inverse(donor.cell_dimension());
        });
    if (iterator == donor.connectivities.end()) {
        throw TopologyError("geometry plan cannot find a reciprocal connection");
    }
    return *iterator;
}

std::vector<BlockId> resolve_donor_path(
    const StructuredMesh& mesh,
    const ConnectivityPatch& connection,
    int required_width)
{
    if (required_width <= 0) {
        return {connection.donor_block};
    }
    std::vector<BlockId> path;
    const ConnectivityPatch* current = &connection;
    int remaining = required_width;
    for (int hop = 0; hop < 64; ++hop) {
        const auto& donor = mesh.block(current->donor_block);
        path.push_back(donor.id());
        const int available = donor.cell_extent()[
            static_cast<std::size_t>(current->donor_face.axis)];
        remaining -= available;
        if (remaining <= 0) {
            return path;
        }
        const FaceLocation exit_face {
            current->donor_face.axis,
            opposite(current->donor_face.side),
        };
        std::vector<const ConnectivityPatch*> candidates;
        for (const auto& candidate : donor.connectivities) {
            if (candidate.receiver_face == exit_face) {
                candidates.push_back(&candidate);
            }
        }
        if (candidates.empty()) {
            throw TopologyError(
                "geometry operand stencil reaches a physical boundary before its halo is complete");
        }
        if (candidates.size() != 1) {
            throw TopologyError(
                "geometry operand stencil has more than one donor propagation path");
        }
        current = candidates.front();
    }
    throw TopologyError("geometry donor path did not terminate within 64 connections");
}

void append_descriptor(
    GeometryHaloPlan& plan,
    const StructuredMesh& mesh,
    const ConnectivityPatch& connection,
    ConnectionId id,
    BlockId owner,
    GeometryMessageKind kind,
    GeometryOperandStage stage,
    int width)
{
    GeometryExchangeDescriptor descriptor;
    descriptor.connection = id;
    descriptor.receiver_block = connection.receiver_block;
    descriptor.donor_block = connection.donor_block;
    descriptor.donor_rank = mesh.block(connection.donor_block).owner_rank();
    descriptor.shared_face_owner = owner;
    descriptor.kind = kind;
    descriptor.stage = stage;
    descriptor.halo_width = width;
    descriptor.index_transform = connection.transform;
    descriptor.periodic = connection.periodic;
    descriptor.donor_path = resolve_donor_path(mesh, connection, width);
    GeometryHaloPlanBuilderAccess::exchanges(plan).push_back(std::move(descriptor));
}

} // namespace

int GeometryExchangeDescriptor::message_tag(int tag_base) const
{
    if (tag_base < 0 || connection < 0 || receiver_block < 0 || donor_block < 0) {
        throw TopologyError("geometry message tag inputs must be non-negative");
    }
    const int direction = receiver_block < donor_block ? 0 : 1;
    const long long tag = static_cast<long long>(tag_base)
        + 32LL * static_cast<long long>(connection)
        + 10LL * static_cast<int>(kind) + 2LL * static_cast<int>(stage)
        + direction;
    if (tag > std::numeric_limits<int>::max()) {
        throw TopologyError("geometry message tag exceeds int range");
    }
    return static_cast<int>(tag);
}

GeometryHaloPlan GeometryHaloPlan::build(
    const StructuredMesh& mesh,
    const AlgorithmProfile& profile)
{
    mesh.validate_connectivities();
    std::vector<const ConnectivityPatch*> canonical;
    std::set<std::tuple<BlockId, int, int, int, int, int, int, int, int>> occupied;
    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            const auto key = std::tuple {
                connection.receiver_block,
                static_cast<int>(connection.receiver_face.axis),
                static_cast<int>(connection.receiver_face.side),
                connection.shared_face_range.begin.i,
                connection.shared_face_range.begin.j,
                connection.shared_face_range.begin.k,
                connection.shared_face_range.end.i,
                connection.shared_face_range.end.j,
                connection.shared_face_range.end.k,
            };
            if (!occupied.insert(key).second) {
                throw TopologyError("two geometry connections occupy the same shared face range");
            }
            if (connection.receiver_block < connection.donor_block) {
                canonical.push_back(&connection);
            }
        }
    }
    std::sort(canonical.begin(), canonical.end(), [](const auto* lhs, const auto* rhs) {
        return connection_key(*lhs) < connection_key(*rhs);
    });

    GeometryHaloPlan plan;
    plan.profile_ = profile.kind();
    const int operand_width
        = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 3 : 5;
    const int final_width
        = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 3 : 3;
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        if (index > static_cast<std::size_t>(std::numeric_limits<ConnectionId>::max())) {
            throw TopologyError("geometry connection count exceeds ConnectionId range");
        }
        const auto id = static_cast<ConnectionId>(index);
        const auto& forward = *canonical[index];
        const auto& reverse = reciprocal_connection(mesh, forward);
        const BlockId owner = std::min(forward.receiver_block, forward.donor_block);
        for (const auto* directed : {&forward, &reverse}) {
            append_descriptor(
                plan, mesh, *directed, id, owner,
                GeometryMessageKind::GeometryVertex,
                GeometryOperandStage::None, 2);
            for (const auto stage : {
                     GeometryOperandStage::CenterCoordinates,
                     GeometryOperandStage::FirstDerivative,
                     GeometryOperandStage::MetricProduct,
                     GeometryOperandStage::JacobianProduct}) {
                append_descriptor(
                    plan, mesh, *directed, id, owner,
                    GeometryMessageKind::GeometryOperand,
                    stage, operand_width);
            }
            append_descriptor(
                plan, mesh, *directed, id, owner,
                GeometryMessageKind::SharedMetric,
                GeometryOperandStage::None, final_width);
        }
    }
    std::set<int> tags;
    for (const auto& exchange : plan.exchanges_) {
        if (!tags.insert(exchange.message_tag()).second) {
            throw TopologyError("geometry message descriptors generated duplicate MPI tags");
        }
    }
    return plan;
}

FaceAreaVectors& SharedMetricSynchronizer::face_vectors(
    MetricField& metric,
    Axis axis)
{
    switch (axis) {
    case Axis::I:
        return metric.i_faces_;
    case Axis::J:
        return metric.j_faces_;
    case Axis::K:
        if (metric.dimension_ != 3) {
            throw TopologyError("2D shared metric cannot be K-face located");
        }
        return metric.k_faces_;
    }
    throw TopologyError("invalid shared metric face axis");
}

void SharedMetricSynchronizer::synchronize(
    const StructuredMesh& mesh,
    std::unordered_map<BlockId, MetricField>& metrics)
{
    mesh.validate_connectivities();
    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            if (connection.receiver_block >= connection.donor_block) {
                continue;
            }
            auto owner_iterator = metrics.find(connection.receiver_block);
            auto donor_iterator = metrics.find(connection.donor_block);
            if (owner_iterator == metrics.end() || donor_iterator == metrics.end()) {
                throw TopologyError("shared metric synchronization is missing a block metric");
            }
            auto& owner_faces
                = face_vectors(owner_iterator->second, connection.receiver_face.axis);
            auto& donor_faces
                = face_vectors(donor_iterator->second, connection.donor_face.axis);
            const auto& reciprocal = reciprocal_connection(mesh, connection);
            const auto counts = connection.shared_face_range.counts();
            const Real orientation = static_cast<Real>(
                -side_sign(connection.receiver_face.side)
                / side_sign(connection.donor_face.side));
            for (int k = 0; k < counts.nk; ++k) {
                for (int j = 0; j < counts.nj; ++j) {
                    for (int i = 0; i < counts.ni; ++i) {
                        const Index3 owner_ordinal {i, j, k};
                        Index3 donor_ordinal;
                        for (int receiver_axis = 0;
                             receiver_axis < block.cell_dimension();
                             ++receiver_axis) {
                            const int donor_axis = std::abs(
                                connection.transform.receiver_to_donor[
                                    static_cast<std::size_t>(receiver_axis)])
                                - 1;
                            donor_ordinal[static_cast<std::size_t>(donor_axis)]
                                = owner_ordinal[static_cast<std::size_t>(receiver_axis)];
                        }
                        const auto owner_index
                            = connection.shared_face_range.at(owner_ordinal);
                        const auto donor_index
                            = reciprocal.shared_face_range.at(donor_ordinal);
                        const std::array<Real, 3> owner_vector {{
                            owner_faces.x(owner_index.i, owner_index.j, owner_index.k),
                            owner_faces.y(owner_index.i, owner_index.j, owner_index.k),
                            owner_faces.z(owner_index.i, owner_index.j, owner_index.k),
                        }};
                        const auto transformed
                            = connection.periodic.apply_vector(owner_vector);
                        donor_faces.x(donor_index.i, donor_index.j, donor_index.k)
                            = orientation * transformed[0];
                        donor_faces.y(donor_index.i, donor_index.j, donor_index.k)
                            = orientation * transformed[1];
                        donor_faces.z(donor_index.i, donor_index.j, donor_index.k)
                            = orientation * transformed[2];
                        const Real owner_area = std::sqrt(
                            owner_vector[0] * owner_vector[0]
                            + owner_vector[1] * owner_vector[1]
                            + owner_vector[2] * owner_vector[2]);
                        const Real donor_area
                            = donor_faces.area(donor_index.i, donor_index.j, donor_index.k);
                        if (std::abs(owner_area - donor_area)
                            > 256.0 * std::numeric_limits<Real>::epsilon()
                                * std::max(Real {1}, owner_area)) {
                            throw TopologyError("published shared metric has an inconsistent area");
                        }
                    }
                }
            }
        }
    }
}

} // namespace wcns
