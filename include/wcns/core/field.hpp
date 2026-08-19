#pragma once

#include <wcns/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace wcns {

// Cell-major field storage: component is the fastest-varying index.
template<class T>
class Field {
public:
    Field() = default;

    Field(Extent3 interior, int components, int ghost_width = 0)
        : interior_(interior)
        , storage_(make_storage_extent(interior, ghost_width))
        , components_(checked_components(components))
        , ghost_width_(ghost_width)
        , data_(checked_value_count(storage_, components_))
    {
    }

    Field(Extent3 interior, int components, int ghost_width, const T& initial_value)
        : Field(interior, components, ghost_width)
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

    [[nodiscard]] int components() const noexcept
    {
        return components_;
    }

    [[nodiscard]] int ghost_width() const noexcept
    {
        return ghost_width_;
    }

    [[nodiscard]] std::size_t cell_count() const noexcept
    {
        return components_ == 0 ? 0 : data_.size() / static_cast<std::size_t>(components_);
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

    T& operator()(int i, int j, int k, int component)
    {
        return data_[linear_index(i, j, k, component)];
    }

    const T& operator()(int i, int j, int k, int component) const
    {
        return data_[linear_index(i, j, k, component)];
    }

    [[nodiscard]] std::size_t cell_linear_index(int i, int j, int k) const
    {
        check_cell_index(i, j, k);
        const auto storage_i = static_cast<std::size_t>(i + ghost_width_);
        const auto storage_j = static_cast<std::size_t>(j + ghost_width_);
        const auto storage_k = static_cast<std::size_t>(k + ghost_width_);
        return (storage_k * static_cast<std::size_t>(storage_.nj) + storage_j)
                * static_cast<std::size_t>(storage_.ni)
            + storage_i;
    }

    [[nodiscard]] std::size_t linear_index(int i, int j, int k, int component) const
    {
        if (component < 0 || component >= components_) {
            throw std::out_of_range("Field component is outside the valid range");
        }
        return cell_linear_index(i, j, k) * static_cast<std::size_t>(components_)
            + static_cast<std::size_t>(component);
    }

    void fill(const T& value)
    {
        std::fill(data_.begin(), data_.end(), value);
    }

private:
    static int checked_components(int components)
    {
        if (components <= 0) {
            throw std::invalid_argument("Field component count must be positive");
        }
        return components;
    }

    static Extent3 make_storage_extent(Extent3 interior, int ghost_width)
    {
        if (!interior.valid()) {
            throw std::invalid_argument("Field interior extent must be non-negative");
        }
        if (ghost_width < 0) {
            throw std::invalid_argument("Field ghost width must be non-negative");
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
            throw std::overflow_error("Field storage extent exceeds integer range");
        }
        return interior_size + 2 * ghost_width;
    }

    static std::size_t checked_value_count(Extent3 storage, int components)
    {
        const auto cells = storage.size();
        const auto component_count = static_cast<std::size_t>(components);
        if (cells > std::numeric_limits<std::size_t>::max() / component_count) {
            throw std::overflow_error("Field value count exceeds addressable range");
        }
        return cells * component_count;
    }

    void check_cell_index(int i, int j, int k) const
    {
        if (i < -ghost_width_ || i >= interior_.ni + ghost_width_
            || j < -ghost_width_ || j >= interior_.nj + ghost_width_
            || k < -ghost_width_ || k >= interior_.nk + ghost_width_) {
            throw std::out_of_range("Field cell index is outside interior and ghost storage");
        }
    }

    Extent3 interior_ {};
    Extent3 storage_ {};
    int components_ = 0;
    int ghost_width_ = 0;
    std::vector<T> data_;
};

} // namespace wcns

