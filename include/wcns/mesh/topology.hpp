#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/types.hpp>

#include <array>
#include <string>
#include <unordered_map>

namespace wcns {

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
    IndexRange3 vertex_range;
    IndexRange3 cell_face_range;
    std::unordered_map<std::string, Real> parameters;
};

// CGNS-style signed receiver-to-donor axis permutation, e.g. {2, -1, 3}.
struct IndexTransform {
    std::array<int, 3> receiver_to_donor {{1, 2, 3}};

    [[nodiscard]] bool valid() const
    {
        std::array<bool, 3> used {{false, false, false}};
        for (const int entry : receiver_to_donor) {
            const int axis = entry < 0 ? -entry : entry;
            if (axis < 1 || axis > 3 || used[static_cast<std::size_t>(axis - 1)]) {
                return false;
            }
            used[static_cast<std::size_t>(axis - 1)] = true;
        }
        return true;
    }
};

struct ConnectivityPatch {
    std::string name;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId donor_rank = invalid_rank_id;
    FaceLocation receiver_face;
    FaceLocation donor_face;
    IndexRange3 receiver_vertex_range;
    IndexRange3 donor_vertex_range;
    IndexRange3 receiver_cell_range;
    IndexRange3 donor_cell_range;
    IndexTransform transform;
    int ghost_width = 0;
};

} // namespace wcns

