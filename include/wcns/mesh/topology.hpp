#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/tagged_range.hpp>
#include <wcns/core/types.hpp>

#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace wcns {

class TopologyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class Axis : int {
    I = 0,
    J = 1,
    K = 2,
};

enum class Side : int {
    Lower = -1,
    Upper = 1,
};

struct FaceLocation {
    Axis axis = Axis::I;
    Side side = Side::Lower;

    friend constexpr bool operator==(const FaceLocation& lhs, const FaceLocation& rhs)
    {
        return lhs.axis == rhs.axis && lhs.side == rhs.side;
    }
};

enum class BoundaryType {
    Undefined,
    Farfield,
    Inflow,
    Outflow,
    SlipWall,
    NoSlipAdiabaticWall,
    NoSlipIsothermalWall,
    Symmetry,
    Periodic,
};

struct BoundaryPatch {
    std::string name;
    BoundaryType type = BoundaryType::Undefined;
    FaceLocation face;
    VertexRange vertex_range;
    AdjacentCellRange adjacent_cell_range;
    BoundaryFaceRange boundary_face_range;
    std::unordered_map<std::string, Real> parameters;
};

// CGNS-style signed receiver-to-donor axis permutation, e.g. {2, -1, 3}.
struct IndexTransform {
    std::array<int, 3> receiver_to_donor {{1, 2, 3}};

    [[nodiscard]] bool valid(int dimension = 3) const
    {
        if (dimension != 2 && dimension != 3) {
            return false;
        }
        std::array<bool, 3> used {{false, false, false}};
        for (int receiver_axis = 0; receiver_axis < dimension; ++receiver_axis) {
            const int entry = receiver_to_donor[static_cast<std::size_t>(receiver_axis)];
            const int axis = std::abs(entry);
            if (axis < 1 || axis > dimension
                || used[static_cast<std::size_t>(axis - 1)]) {
                return false;
            }
            used[static_cast<std::size_t>(axis - 1)] = true;
        }
        for (int axis = dimension; axis < 3; ++axis) {
            if (receiver_to_donor[static_cast<std::size_t>(axis)] != axis + 1) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] Index3 map(
        Index3 receiver,
        Index3 receiver_origin,
        Index3 donor_origin,
        int dimension) const
    {
        if (!valid(dimension)) {
            throw std::invalid_argument("cannot apply an invalid index transform");
        }
        Index3 donor = donor_origin;
        for (int receiver_axis = 0; receiver_axis < dimension; ++receiver_axis) {
            const int entry = receiver_to_donor[static_cast<std::size_t>(receiver_axis)];
            const auto donor_axis = static_cast<std::size_t>(std::abs(entry) - 1);
            const int sign = entry < 0 ? -1 : 1;
            donor[donor_axis] += sign
                * (receiver[static_cast<std::size_t>(receiver_axis)]
                    - receiver_origin[static_cast<std::size_t>(receiver_axis)]);
        }
        return donor;
    }

    [[nodiscard]] IndexTransform inverse(int dimension) const
    {
        if (!valid(dimension)) {
            throw std::invalid_argument("cannot invert an invalid index transform");
        }
        IndexTransform result;
        for (int receiver_axis = 0; receiver_axis < dimension; ++receiver_axis) {
            const int entry = receiver_to_donor[static_cast<std::size_t>(receiver_axis)];
            const auto donor_axis = static_cast<std::size_t>(std::abs(entry) - 1);
            result.receiver_to_donor[donor_axis]
                = (entry < 0 ? -1 : 1) * (receiver_axis + 1);
        }
        return result;
    }

    friend bool operator==(const IndexTransform& lhs, const IndexTransform& rhs)
    {
        return lhs.receiver_to_donor == rhs.receiver_to_donor;
    }
};

struct PeriodicTransform {
    std::array<std::array<Real, 3>, 3> rotation {{
        {{1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 1.0}},
    }};
    std::array<Real, 3> translation {{0.0, 0.0, 0.0}};

    [[nodiscard]] bool valid(int dimension) const;
    [[nodiscard]] std::array<Real, 3> apply_point(
        const std::array<Real, 3>& point) const;
    [[nodiscard]] std::array<Real, 3> apply_vector(
        const std::array<Real, 3>& vector) const;
    [[nodiscard]] PeriodicTransform inverse() const;

    friend bool operator==(const PeriodicTransform& lhs, const PeriodicTransform& rhs)
    {
        return lhs.rotation == rhs.rotation && lhs.translation == rhs.translation;
    }
};

struct ConnectivityPatch {
    std::string name;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId donor_rank = invalid_rank_id;
    FaceLocation receiver_face;
    FaceLocation donor_face;
    ReceiverVertexRange receiver_vertex_range;
    DonorVertexRange donor_vertex_range;
    ReceiverAdjacentCellRange receiver_adjacent_cell_range;
    DonorAdjacentCellRange donor_adjacent_cell_range;
    SharedFaceRange shared_face_range;
    IndexTransform transform;
    int ghost_width = 0;
    ConnectionId id = invalid_connection_id;
    PeriodicTransform periodic {};
};

} // namespace wcns
