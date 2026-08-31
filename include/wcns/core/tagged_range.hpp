#pragma once

#include <wcns/core/index.hpp>

namespace wcns {

// A directed inclusive range whose tag prevents accidental interchange with a
// range from another topological location.
template<class Tag>
class TaggedIndexRange3 {
public:
    constexpr TaggedIndexRange3() = default;

    constexpr TaggedIndexRange3(Index3 begin, Index3 end)
        : begin(begin)
        , end(end)
    {
    }

    explicit constexpr TaggedIndexRange3(IndexRange3 range)
        : begin(range.begin)
        , end(range.end)
    {
    }

    [[nodiscard]] constexpr Index3 step() const
    {
        return untyped().step();
    }

    [[nodiscard]] constexpr Extent3 counts() const
    {
        return untyped().counts();
    }

    [[nodiscard]] constexpr std::size_t size() const
    {
        return untyped().size();
    }

    [[nodiscard]] constexpr Index3 at(Index3 ordinal) const
    {
        return untyped().at(ordinal);
    }

    [[nodiscard]] constexpr IndexRange3 untyped() const
    {
        return {begin, end};
    }

    friend constexpr bool operator==(
        const TaggedIndexRange3& lhs,
        const TaggedIndexRange3& rhs)
    {
        return lhs.begin == rhs.begin && lhs.end == rhs.end;
    }

    Index3 begin {};
    Index3 end {};
};

struct VertexRangeTag;
struct AdjacentCellRangeTag;
struct BoundaryFaceRangeTag;
struct ReceiverVertexRangeTag;
struct DonorVertexRangeTag;
struct ReceiverAdjacentCellRangeTag;
struct DonorAdjacentCellRangeTag;
struct SharedFaceRangeTag;

using VertexRange = TaggedIndexRange3<VertexRangeTag>;
using AdjacentCellRange = TaggedIndexRange3<AdjacentCellRangeTag>;
using BoundaryFaceRange = TaggedIndexRange3<BoundaryFaceRangeTag>;
using ReceiverVertexRange = TaggedIndexRange3<ReceiverVertexRangeTag>;
using DonorVertexRange = TaggedIndexRange3<DonorVertexRangeTag>;
using ReceiverAdjacentCellRange = TaggedIndexRange3<ReceiverAdjacentCellRangeTag>;
using DonorAdjacentCellRange = TaggedIndexRange3<DonorAdjacentCellRangeTag>;
using SharedFaceRange = TaggedIndexRange3<SharedFaceRangeTag>;

} // namespace wcns
