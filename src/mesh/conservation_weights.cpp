#include <wcns/mesh/conservation_weights.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace wcns {
namespace {

std::vector<Real> solve_dense(
    std::vector<std::vector<Real>> matrix,
    std::vector<Real> rhs)
{
    const int count = static_cast<int>(rhs.size());
    for (int column = 0; column < count; ++column) {
        int pivot_row = column;
        for (int row = column + 1; row < count; ++row) {
            if (std::abs(matrix[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)])
                > std::abs(matrix[static_cast<std::size_t>(pivot_row)]
                                  [static_cast<std::size_t>(column)])) {
                pivot_row = row;
            }
        }
        const Real pivot = matrix[static_cast<std::size_t>(pivot_row)]
                                 [static_cast<std::size_t>(column)];
        if (!std::isfinite(pivot)
            || std::abs(pivot) < 1024.0 * std::numeric_limits<Real>::epsilon()) {
            throw ConservationWeightError(
                "conservation-weight system does not have a unique solution");
        }
        std::swap(matrix[static_cast<std::size_t>(column)],
            matrix[static_cast<std::size_t>(pivot_row)]);
        std::swap(rhs[static_cast<std::size_t>(column)],
            rhs[static_cast<std::size_t>(pivot_row)]);
        for (int row = column + 1; row < count; ++row) {
            const Real factor = matrix[static_cast<std::size_t>(row)]
                                      [static_cast<std::size_t>(column)]
                / matrix[static_cast<std::size_t>(column)]
                        [static_cast<std::size_t>(column)];
            for (int inner = column; inner < count; ++inner) {
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(inner)]
                    -= factor * matrix[static_cast<std::size_t>(column)]
                                       [static_cast<std::size_t>(inner)];
            }
            rhs[static_cast<std::size_t>(row)]
                -= factor * rhs[static_cast<std::size_t>(column)];
        }
    }
    std::vector<Real> result(static_cast<std::size_t>(count), 0.0);
    for (int row = count - 1; row >= 0; --row) {
        Real value = rhs[static_cast<std::size_t>(row)];
        for (int column = row + 1; column < count; ++column) {
            value -= matrix[static_cast<std::size_t>(row)]
                           [static_cast<std::size_t>(column)]
                * result[static_cast<std::size_t>(column)];
        }
        result[static_cast<std::size_t>(row)]
            = value / matrix[static_cast<std::size_t>(row)]
                             [static_cast<std::size_t>(row)];
    }
    return result;
}

LineConservationWeights build_uncached_line(
    const AlgorithmProfile& profile,
    int cell_count,
    bool periodic)
{
    LineConservationWeights result;
    result.profile = profile.kind();
    result.cell_count = cell_count;
    result.periodic = periodic;
    if (periodic) {
        if (cell_count <= 0) {
            throw ConservationWeightError("periodic line must contain cells");
        }
        result.cell_weights.assign(static_cast<std::size_t>(cell_count), 1.0);
        return result;
    }

    const auto operators = LineOperators::build(profile, cell_count);
    const auto& rows = operators.derivative_rows();
    std::vector<std::vector<Real>> transpose(
        static_cast<std::size_t>(cell_count + 1),
        std::vector<Real>(static_cast<std::size_t>(cell_count), 0.0));
    for (int cell = 0; cell < cell_count; ++cell) {
        for (const auto [face, coefficient] : rows[static_cast<std::size_t>(cell)]) {
            transpose[static_cast<std::size_t>(face)]
                     [static_cast<std::size_t>(cell)] = coefficient;
        }
    }
    std::vector<Real> boundary(static_cast<std::size_t>(cell_count + 1), 0.0);
    boundary.front() = -1.0;
    boundary.back() = 1.0;

    std::vector<std::vector<Real>> normal(
        static_cast<std::size_t>(cell_count),
        std::vector<Real>(static_cast<std::size_t>(cell_count), 0.0));
    std::vector<Real> normal_rhs(static_cast<std::size_t>(cell_count), 0.0);
    for (int row = 0; row < cell_count; ++row) {
        for (int column = 0; column < cell_count; ++column) {
            for (int face = 0; face <= cell_count; ++face) {
                normal[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)]
                    += transpose[static_cast<std::size_t>(face)]
                                [static_cast<std::size_t>(row)]
                    * transpose[static_cast<std::size_t>(face)]
                               [static_cast<std::size_t>(column)];
            }
        }
        for (int face = 0; face <= cell_count; ++face) {
            normal_rhs[static_cast<std::size_t>(row)]
                += transpose[static_cast<std::size_t>(face)]
                            [static_cast<std::size_t>(row)]
                * boundary[static_cast<std::size_t>(face)];
        }
    }
    result.cell_weights = solve_dense(std::move(normal), std::move(normal_rhs));
    for (const Real weight : result.cell_weights) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            throw ConservationWeightError(
                "conservation-weight system produced a non-positive cell weight");
        }
    }
    for (int face = 0; face <= cell_count; ++face) {
        Real actual = 0.0;
        for (int cell = 0; cell < cell_count; ++cell) {
            actual += transpose[static_cast<std::size_t>(face)]
                               [static_cast<std::size_t>(cell)]
                * result.cell_weights[static_cast<std::size_t>(cell)];
        }
        result.maximum_residual = std::max(
            result.maximum_residual,
            std::abs(actual - boundary[static_cast<std::size_t>(face)]));
    }
    if (result.maximum_residual > 1.0e-10) {
        throw ConservationWeightError(
            "conservation-weight residual exceeds the initialization tolerance");
    }
    return result;
}

const ConnectivityPatch& reciprocal_connection(
    const StructuredMesh& mesh,
    const ConnectivityPatch& connection)
{
    const auto& donor = mesh.block(connection.donor_block);
    const auto iterator = std::find_if(
        donor.connectivities.begin(), donor.connectivities.end(),
        [&](const ConnectivityPatch& candidate) {
            return candidate.receiver_block == connection.donor_block
                && candidate.donor_block == connection.receiver_block
                && candidate.receiver_face == connection.donor_face;
        });
    if (iterator == donor.connectivities.end()) {
        throw ConservationWeightError(
            "global conservation assembly cannot find a reciprocal connection");
    }
    return *iterator;
}

Real tangential_weight(
    const std::array<LineConservationWeights, 3>& lines,
    FaceLocation face,
    Index3 index,
    int dimension)
{
    Real result = 1.0;
    for (int axis = 0; axis < dimension; ++axis) {
        if (axis != static_cast<int>(face.axis)) {
            result *= lines[static_cast<std::size_t>(axis)]
                          .cell_weights[static_cast<std::size_t>(
                              index[static_cast<std::size_t>(axis)])];
        }
    }
    return result;
}

Array3D<Real>& boundary_field(BlockConservationWeights& block, Axis axis)
{
    switch (axis) {
    case Axis::I:
        return block.physical_i_face;
    case Axis::J:
        return block.physical_j_face;
    case Axis::K:
        return block.physical_k_face;
    }
    throw ConservationWeightError("invalid boundary-weight face axis");
}

} // namespace

LineConservationWeights build_line_conservation_weights(
    const AlgorithmProfile& profile,
    int cell_count,
    bool periodic)
{
    using CacheKey = std::tuple<int, int, bool>;
    static std::map<CacheKey, LineConservationWeights> cache;
    const CacheKey key {
        static_cast<int>(profile.kind()), cell_count, periodic};
    const auto iterator = cache.find(key);
    if (iterator != cache.end()) {
        return iterator->second;
    }
    auto result = build_uncached_line(profile, cell_count, periodic);
    cache.emplace(key, result);
    return result;
}

GlobalConservationWeights GlobalConservationWeights::build(
    const StructuredMesh& mesh,
    const AlgorithmProfile& profile)
{
    mesh.validate_connectivities();
    GlobalConservationWeights result;
    result.profile_ = profile.kind();
    std::unordered_map<BlockId, std::array<LineConservationWeights, 3>> line_weights;
    for (const auto& block : mesh.blocks()) {
        std::array<LineConservationWeights, 3> lines;
        for (int axis = 0; axis < block.cell_dimension(); ++axis) {
            lines[static_cast<std::size_t>(axis)] = build_line_conservation_weights(
                profile,
                block.cell_extent()[static_cast<std::size_t>(axis)]);
            result.maximum_line_residual_ = std::max(
                result.maximum_line_residual_,
                lines[static_cast<std::size_t>(axis)].maximum_residual);
        }
        if (block.cell_dimension() == 2) {
            lines[2] = build_line_conservation_weights(profile, 1, true);
        }
        line_weights.emplace(block.id(), lines);
        auto [iterator, inserted]
            = result.blocks_.emplace(block.id(), BlockConservationWeights(block.cell_extent()));
        static_cast<void>(inserted);
        auto& block_weights = iterator->second;
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    block_weights.cell(i, j, k)
                        = lines[0].cell_weights[static_cast<std::size_t>(i)]
                        * lines[1].cell_weights[static_cast<std::size_t>(j)]
                        * lines[2].cell_weights[static_cast<std::size_t>(k)];
                }
            }
        }
        for (const auto& patch : block.boundaries) {
            auto& faces = boundary_field(block_weights, patch.face.axis);
            const auto counts = patch.boundary_face_range.counts();
            for (int k = 0; k < counts.nk; ++k) {
                for (int j = 0; j < counts.nj; ++j) {
                    for (int i = 0; i < counts.ni; ++i) {
                        const auto face = patch.boundary_face_range.at({i, j, k});
                        faces(face.i, face.j, face.k) = tangential_weight(
                            lines, patch.face, face, block.cell_dimension());
                    }
                }
            }
        }
    }

    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            if (connection.receiver_block >= connection.donor_block) {
                continue;
            }
            const auto& reciprocal = reciprocal_connection(mesh, connection);
            const auto& receiver_lines = line_weights.at(connection.receiver_block);
            const auto& donor_lines = line_weights.at(connection.donor_block);
            const auto counts = connection.shared_face_range.counts();
            for (int k = 0; k < counts.nk; ++k) {
                for (int j = 0; j < counts.nj; ++j) {
                    for (int i = 0; i < counts.ni; ++i) {
                        const Index3 receiver_ordinal {i, j, k};
                        Index3 donor_ordinal;
                        for (int receiver_axis = 0;
                             receiver_axis < block.cell_dimension();
                             ++receiver_axis) {
                            const int donor_axis = std::abs(
                                connection.transform.receiver_to_donor[
                                    static_cast<std::size_t>(receiver_axis)])
                                - 1;
                            donor_ordinal[static_cast<std::size_t>(donor_axis)]
                                = receiver_ordinal[static_cast<std::size_t>(receiver_axis)];
                        }
                        const auto receiver_face
                            = connection.shared_face_range.at(receiver_ordinal);
                        const auto donor_face
                            = reciprocal.shared_face_range.at(donor_ordinal);
                        const Real receiver_weight = tangential_weight(
                            receiver_lines,
                            connection.receiver_face,
                            receiver_face,
                            block.cell_dimension());
                        const Real donor_weight = tangential_weight(
                            donor_lines,
                            connection.donor_face,
                            donor_face,
                            block.cell_dimension());
                        result.maximum_shared_face_mismatch_ = std::max(
                            result.maximum_shared_face_mismatch_,
                            std::abs(receiver_weight - donor_weight));
                    }
                }
            }
    }
    }
    if (result.maximum_shared_face_mismatch_ > 1.0e-10) {
        throw ConservationWeightError(
            "global conservation weights do not cancel on a shared connection face");
    }
    return result;
}

const BlockConservationWeights& GlobalConservationWeights::block(BlockId id) const
{
    const auto iterator = blocks_.find(id);
    if (iterator == blocks_.end()) {
        throw std::out_of_range("conservation weights do not contain the requested block");
    }
    return iterator->second;
}

} // namespace wcns
