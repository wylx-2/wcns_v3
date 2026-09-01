#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>

#include <cstdint>

namespace wcns {

inline constexpr int euler_components = 5;

struct FlowFields {
    FlowFields(Extent3 cell_extent, int ghost_width)
        : conservative(cell_extent, euler_components, ghost_width)
        , primitive(cell_extent, euler_components, ghost_width)
        , temperature_primitive(cell_extent, euler_components, ghost_width)
        , residual(cell_extent, euler_components, 0)
    {
    }

    Field<Real> conservative;
    Field<Real> primitive;
    Field<Real> temperature_primitive;
    Field<Real> residual;
    std::uint64_t physical_ghost_version = 0;
};

} // namespace wcns
