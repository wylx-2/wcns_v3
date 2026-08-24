#include <wcns/mesh/structured_mesh.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace wcns {
namespace {

bool same_undirected_range(const IndexRange3& lhs, const IndexRange3& rhs)
{
    return lhs == rhs || (lhs.begin == rhs.end && lhs.end == rhs.begin);
}

bool is_reciprocal(
    const ConnectivityPatch& connection,
    const ConnectivityPatch& candidate,
    int dimension)
{
    return candidate.receiver_block == connection.donor_block
        && candidate.donor_block == connection.receiver_block
        && candidate.receiver_face == connection.donor_face
        && candidate.donor_face == connection.receiver_face
        && same_undirected_range(
            candidate.receiver_vertex_range, connection.donor_vertex_range)
        && same_undirected_range(
            candidate.donor_vertex_range, connection.receiver_vertex_range)
        && same_undirected_range(
            candidate.receiver_cell_range, connection.donor_cell_range)
        && same_undirected_range(
            candidate.donor_cell_range, connection.receiver_cell_range)
        && candidate.transform == connection.transform.inverse(dimension);
}

bool coordinates_match(Real receiver, Real donor)
{
    const Real scale = std::max({Real {1}, std::abs(receiver), std::abs(donor)});
    return std::abs(receiver - donor)
        <= Real {256} * std::numeric_limits<Real>::epsilon() * scale;
}

void validate_interface_coordinates(
    const StructuredBlock& receiver,
    const StructuredBlock& donor,
    const ConnectivityPatch& connection)
{
    const auto counts = connection.receiver_vertex_range.counts();
    for (int k = 0; k < counts.nk; ++k) {
        for (int j = 0; j < counts.nj; ++j) {
            for (int i = 0; i < counts.ni; ++i) {
                const auto receiver_index
                    = connection.receiver_vertex_range.at({i, j, k});
                const auto donor_index = connection.transform.map(
                    receiver_index,
                    connection.receiver_vertex_range.begin,
                    connection.donor_vertex_range.begin,
                    receiver.cell_dimension());
                if (!coordinates_match(
                        receiver.coordinates.x(
                            receiver_index.i, receiver_index.j, receiver_index.k),
                        donor.coordinates.x(donor_index.i, donor_index.j, donor_index.k))
                    || !coordinates_match(
                        receiver.coordinates.y(
                            receiver_index.i, receiver_index.j, receiver_index.k),
                        donor.coordinates.y(donor_index.i, donor_index.j, donor_index.k))
                    || !coordinates_match(
                        receiver.coordinates.z(
                            receiver_index.i, receiver_index.j, receiver_index.k),
                        donor.coordinates.z(donor_index.i, donor_index.j, donor_index.k))) {
                    throw TopologyError(
                        "connectivity " + connection.name
                        + " maps vertices with different physical coordinates");
                }
            }
        }
    }
}

} // namespace

StructuredMesh::StructuredMesh(std::vector<StructuredBlock> blocks)
    : blocks_(std::move(blocks))
{
    block_index_.reserve(blocks_.size());
    for (std::size_t index = 0; index < blocks_.size(); ++index) {
        const auto [iterator, inserted] = block_index_.emplace(blocks_[index].id(), index);
        static_cast<void>(iterator);
        if (!inserted) {
            throw TopologyError(
                "duplicate structured block id " + std::to_string(blocks_[index].id()));
        }
    }
}

bool StructuredMesh::contains(BlockId id) const noexcept
{
    return block_index_.find(id) != block_index_.end();
}

StructuredBlock& StructuredMesh::block(BlockId id)
{
    const auto iterator = block_index_.find(id);
    if (iterator == block_index_.end()) {
        throw std::out_of_range("structured block id is not present in the mesh");
    }
    return blocks_[iterator->second];
}

const StructuredBlock& StructuredMesh::block(BlockId id) const
{
    const auto iterator = block_index_.find(id);
    if (iterator == block_index_.end()) {
        throw std::out_of_range("structured block id is not present in the mesh");
    }
    return blocks_[iterator->second];
}

void StructuredMesh::validate_connectivities() const
{
    for (const auto& receiver : blocks_) {
        for (const auto& connection : receiver.connectivities) {
            if (connection.receiver_block != receiver.id()) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " has a receiver id inconsistent with its owning block");
            }
            if (!contains(connection.donor_block)) {
                throw TopologyError(
                    "connectivity " + connection.name + " references an unknown donor block");
            }
            const auto& donor = block(connection.donor_block);
            if (receiver.cell_dimension() != donor.cell_dimension()) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " joins blocks with different cell dimensions");
            }
            const int dimension = receiver.cell_dimension();
            if (connection.donor_rank != donor.owner_rank()) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " has a donor rank inconsistent with its donor block");
            }
            if (connection.ghost_width != receiver.ghost_width()) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " has a ghost width inconsistent with its receiver block");
            }
            if (!connection.transform.valid(dimension)) {
                throw TopologyError(
                    "connectivity " + connection.name + " has an invalid index transform");
            }
            if (connection.transform.map(
                    connection.receiver_vertex_range.end,
                    connection.receiver_vertex_range.begin,
                    connection.donor_vertex_range.begin,
                    dimension)
                != connection.donor_vertex_range.end) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " transform does not map the receiver range onto the donor range");
            }
            const auto reciprocal = std::find_if(
                donor.connectivities.begin(),
                donor.connectivities.end(),
                [&](const ConnectivityPatch& candidate) {
                    return is_reciprocal(connection, candidate, dimension);
                });
            if (reciprocal == donor.connectivities.end()) {
                throw TopologyError(
                    "connectivity " + connection.name + " has no reciprocal donor record");
            }
            validate_interface_coordinates(receiver, donor, connection);
        }
    }
}

} // namespace wcns
