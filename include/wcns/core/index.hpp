#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace wcns {

struct Index3 {
    int i = 0;
    int j = 0;
    int k = 0;

    constexpr int& operator[](std::size_t axis)
    {
        switch (axis) {
        case 0:
            return i;
        case 1:
            return j;
        case 2:
            return k;
        default:
            throw std::out_of_range("Index3 axis must be in [0, 2]");
        }
    }

    constexpr const int& operator[](std::size_t axis) const
    {
        switch (axis) {
        case 0:
            return i;
        case 1:
            return j;
        case 2:
            return k;
        default:
            throw std::out_of_range("Index3 axis must be in [0, 2]");
        }
    }

    friend constexpr bool operator==(const Index3& lhs, const Index3& rhs)
    {
        return lhs.i == rhs.i && lhs.j == rhs.j && lhs.k == rhs.k;
    }

    friend constexpr bool operator!=(const Index3& lhs, const Index3& rhs)
    {
        return !(lhs == rhs);
    }
};

struct Extent3 {
    int ni = 0;
    int nj = 0;
    int nk = 0;

    constexpr int& operator[](std::size_t axis)
    {
        switch (axis) {
        case 0:
            return ni;
        case 1:
            return nj;
        case 2:
            return nk;
        default:
            throw std::out_of_range("Extent3 axis must be in [0, 2]");
        }
    }

    constexpr const int& operator[](std::size_t axis) const
    {
        switch (axis) {
        case 0:
            return ni;
        case 1:
            return nj;
        case 2:
            return nk;
        default:
            throw std::out_of_range("Extent3 axis must be in [0, 2]");
        }
    }

    [[nodiscard]] constexpr bool valid() const
    {
        return ni >= 0 && nj >= 0 && nk >= 0;
    }

    [[nodiscard]] constexpr bool empty() const
    {
        return ni == 0 || nj == 0 || nk == 0;
    }

    [[nodiscard]] constexpr std::size_t size() const
    {
        if (!valid()) {
            throw std::invalid_argument("Extent3 dimensions must be non-negative");
        }
        const auto ij = checked_product(
            static_cast<std::size_t>(ni), static_cast<std::size_t>(nj));
        return checked_product(ij, static_cast<std::size_t>(nk));
    }

    friend constexpr bool operator==(const Extent3& lhs, const Extent3& rhs)
    {
        return lhs.ni == rhs.ni && lhs.nj == rhs.nj && lhs.nk == rhs.nk;
    }

    friend constexpr bool operator!=(const Extent3& lhs, const Extent3& rhs)
    {
        return !(lhs == rhs);
    }

private:
    static constexpr std::size_t checked_product(std::size_t lhs, std::size_t rhs)
    {
        if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
            throw std::overflow_error("Extent3 size exceeds addressable range");
        }
        return lhs * rhs;
    }
};

// Inclusive, directed index range. Each axis may increase, decrease, or stay fixed.
struct IndexRange3 {
    Index3 begin;
    Index3 end;

    [[nodiscard]] constexpr Index3 step() const
    {
        return {
            direction(begin.i, end.i),
            direction(begin.j, end.j),
            direction(begin.k, end.k),
        };
    }

    [[nodiscard]] constexpr Extent3 counts() const
    {
        return {
            distance(begin.i, end.i) + 1,
            distance(begin.j, end.j) + 1,
            distance(begin.k, end.k) + 1,
        };
    }

    [[nodiscard]] constexpr std::size_t size() const
    {
        return counts().size();
    }

    [[nodiscard]] constexpr Index3 at(Index3 ordinal) const
    {
        const auto count = counts();
        if (ordinal.i < 0 || ordinal.i >= count.ni || ordinal.j < 0
            || ordinal.j >= count.nj || ordinal.k < 0 || ordinal.k >= count.nk) {
            throw std::out_of_range("IndexRange3 ordinal is outside the range");
        }
        const auto stride = step();
        return {
            begin.i + ordinal.i * stride.i,
            begin.j + ordinal.j * stride.j,
            begin.k + ordinal.k * stride.k,
        };
    }

    friend constexpr bool operator==(const IndexRange3& lhs, const IndexRange3& rhs)
    {
        return lhs.begin == rhs.begin && lhs.end == rhs.end;
    }

    friend constexpr bool operator!=(const IndexRange3& lhs, const IndexRange3& rhs)
    {
        return !(lhs == rhs);
    }

private:
    static constexpr int direction(int first, int last)
    {
        return (last > first) - (last < first);
    }

    static constexpr int distance(int first, int last)
    {
        return last >= first ? last - first : first - last;
    }
};

} // namespace wcns
