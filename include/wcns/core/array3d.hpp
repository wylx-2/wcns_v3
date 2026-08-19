#pragma once

#include <wcns/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wcns {

template<class T>
class Array3D {
public:
    Array3D() = default;

    explicit Array3D(Extent3 interior, int ghost_width = 0)
        : interior_(interior)
        , storage_(make_storage_extent(interior, ghost_width))
        , ghost_width_(ghost_width)
        , data_(storage_.size())
    {
    }

    Array3D(Extent3 interior, int ghost_width, const T& initial_value)
        : Array3D(interior, ghost_width)
    {
        fill(initial_value);
    }

    [[nodiscard]] const Extent3& interior_extent() const noexcept
    {
        return interior_;
    }

    [[nodiscard]] const Extent3& storage_extent() const noexcept
    {
        return storage_;
    }

    [[nodiscard]] int ghost_width() const noexcept
    {
        return ghost_width_;
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return data_.size();
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return data_.empty();
    }

    T* data() noexcept
    {
        return data_.data();
    }

    const T* data() const noexcept
    {
        return data_.data();
    }

    T& operator()(int i, int j, int k)
    {
        return data_[linear_index(i, j, k)];
    }

    const T& operator()(int i, int j, int k) const
    {
        return data_[linear_index(i, j, k)];
    }

    [[nodiscard]] std::size_t linear_index(int i, int j, int k) const
    {
        check_index(i, j, k);
        const auto storage_i = static_cast<std::size_t>(i + ghost_width_);
        const auto storage_j = static_cast<std::size_t>(j + ghost_width_);
        const auto storage_k = static_cast<std::size_t>(k + ghost_width_);
        return (storage_k * static_cast<std::size_t>(storage_.nj) + storage_j)
                * static_cast<std::size_t>(storage_.ni)
            + storage_i;
    }

    void fill(const T& value)
    {
        std::fill(data_.begin(), data_.end(), value);
    }

private:
    static Extent3 make_storage_extent(Extent3 interior, int ghost_width)
    {
        if (!interior.valid()) {
            throw std::invalid_argument("Array3D interior extent must be non-negative");
        }
        if (ghost_width < 0) {
            throw std::invalid_argument("Array3D ghost width must be non-negative");
        }

        return {
            add_ghost_layers(interior.ni, ghost_width),
            add_ghost_layers(interior.nj, ghost_width),
            add_ghost_layers(interior.nk, ghost_width),
        };
    }

    static int add_ghost_layers(int interior_size, int ghost_width)
    {
        if (ghost_width > (std::numeric_limits<int>::max() - interior_size) / 2) {
            throw std::overflow_error("Array3D storage extent exceeds integer range");
        }
        return interior_size + 2 * ghost_width;
    }

    void check_index(int i, int j, int k) const
    {
        if (i < -ghost_width_ || i >= interior_.ni + ghost_width_
            || j < -ghost_width_ || j >= interior_.nj + ghost_width_
            || k < -ghost_width_ || k >= interior_.nk + ghost_width_) {
            throw std::out_of_range("Array3D index is outside interior and ghost storage");
        }
    }

    Extent3 interior_ {};
    Extent3 storage_ {};
    int ghost_width_ = 0;
    std::vector<T> data_;
};

} // namespace wcns

