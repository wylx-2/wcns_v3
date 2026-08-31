#pragma once

#include <wcns/core/field.hpp>

#include <stdexcept>

namespace wcns {

enum class TopologyLocation {
    Vertex,
    Cell,
    FaceI,
    FaceJ,
    FaceK,
};

enum class AccessRegion {
    Interior,
    ConnectionHalo,
    PhysicalBoundarySlab,
};

struct FieldAccessPermissions {
    bool connection_halo = false;
    bool physical_boundary_slab = false;
};

class FieldDomain {
public:
    FieldDomain(
        Extent3 interior,
        int ghost_width,
        int dimension,
        TopologyLocation location,
        FieldAccessPermissions permissions = {})
        : interior_(interior)
        , ghost_width_(ghost_width)
        , dimension_(dimension)
        , location_(location)
        , permissions_(permissions)
    {
        if (dimension != 2 && dimension != 3) {
            throw std::invalid_argument("field domain dimension must be 2 or 3");
        }
        if (!interior.valid() || interior.empty()) {
            throw std::invalid_argument("field domain interior extent must be positive");
        }
        if (ghost_width < 0) {
            throw std::invalid_argument("field domain ghost width must be non-negative");
        }
        if (dimension == 2 && interior.nk != 1) {
            throw std::invalid_argument("a 2D field domain must have one K plane");
        }
        if (location == TopologyLocation::FaceK && dimension == 2) {
            throw std::invalid_argument("a 2D field domain cannot be K-face located");
        }
    }

    [[nodiscard]] const Extent3& interior_extent() const noexcept
    {
        return interior_;
    }

    [[nodiscard]] int ghost_width() const noexcept { return ghost_width_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] TopologyLocation location() const noexcept { return location_; }

    void validate(Index3 index, AccessRegion region) const
    {
        if (dimension_ == 2 && index.k != 0) {
            throw std::out_of_range(
                "2D field access cannot use a nonzero inactive K index");
        }

        int outside_active_axes = 0;
        for (int axis = 0; axis < dimension_; ++axis) {
            const auto axis_index = static_cast<std::size_t>(axis);
            const int coordinate = index[axis_index];
            const int extent = interior_[axis_index];
            if (coordinate < -ghost_width_ || coordinate >= extent + ghost_width_) {
                throw std::out_of_range("field access is outside allocated ghost storage");
            }
            if (coordinate < 0 || coordinate >= extent) {
                ++outside_active_axes;
            }
        }

        switch (region) {
        case AccessRegion::Interior:
            if (outside_active_axes != 0) {
                throw std::out_of_range("interior field access addresses a ghost index");
            }
            return;
        case AccessRegion::ConnectionHalo:
            if (!permissions_.connection_halo) {
                throw std::logic_error("connection-halo access is disabled for this field");
            }
            break;
        case AccessRegion::PhysicalBoundarySlab:
            if (!permissions_.physical_boundary_slab) {
                throw std::logic_error(
                    "physical-boundary ghost access is disabled for this field");
            }
            break;
        }

        // WCNS face stencils use face slabs only. Edge and corner ghosts are
        // deliberately outside the legal domain even if Field owns storage there.
        if (outside_active_axes != 1) {
            throw std::out_of_range(
                "halo access must address a face slab, not an interior, edge, or corner");
        }
    }

private:
    Extent3 interior_ {};
    int ghost_width_ = 0;
    int dimension_ = 0;
    TopologyLocation location_ = TopologyLocation::Cell;
    FieldAccessPermissions permissions_ {};
};

template<class T>
class TopologyField {
public:
    TopologyField(FieldDomain domain, int components)
        : domain_(domain)
        , storage_(domain.interior_extent(), components, domain.ghost_width())
    {
    }

    TopologyField(FieldDomain domain, int components, const T& initial_value)
        : domain_(domain)
        , storage_(
              domain.interior_extent(),
              components,
              domain.ghost_width(),
              initial_value)
    {
    }

    [[nodiscard]] const FieldDomain& domain() const noexcept { return domain_; }
    [[nodiscard]] int components() const noexcept { return storage_.components(); }

    T& at(Index3 index, int component, AccessRegion region)
    {
        domain_.validate(index, region);
        return storage_(index.i, index.j, index.k, component);
    }

    const T& at(Index3 index, int component, AccessRegion region) const
    {
        domain_.validate(index, region);
        return storage_(index.i, index.j, index.k, component);
    }

    void fill(const T& value) { storage_.fill(value); }

private:
    FieldDomain domain_;
    Field<T> storage_;
};

} // namespace wcns
