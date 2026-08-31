#include "test_support.hpp"

#include <wcns/core/tagged_range.hpp>
#include <wcns/core/topology_field.hpp>

#include <stdexcept>
#include <type_traits>

// 验收强类型拓扑范围和不同 ghost 合法域的访问约束。
void test_topology_field()
{
    using namespace wcns;

    static_assert(!std::is_convertible_v<AdjacentCellRange, BoundaryFaceRange>);
    static_assert(!std::is_constructible_v<BoundaryFaceRange, AdjacentCellRange>);
    static_assert(!std::is_convertible_v<ReceiverVertexRange, DonorVertexRange>);

    const FieldDomain state_domain(
        {4, 3, 1},
        3,
        2,
        TopologyLocation::Cell,
        {.connection_halo = true, .physical_boundary_slab = true});
    TopologyField<double> state(state_domain, 5, 0.0);
    state.at({1, 1, 0}, 0, AccessRegion::Interior) = 2.0;
    state.at({-1, 1, 0}, 0, AccessRegion::PhysicalBoundarySlab) = 3.0;
    state.at({4, 1, 0}, 0, AccessRegion::ConnectionHalo) = 4.0;
    WCNS_REQUIRE(state.at({1, 1, 0}, 0, AccessRegion::Interior) == 2.0);
    WCNS_REQUIRE(
        state.at({-1, 1, 0}, 0, AccessRegion::PhysicalBoundarySlab) == 3.0);
    WCNS_REQUIRE_THROWS(
        std::out_of_range,
        state.at({-1, -1, 0}, 0, AccessRegion::PhysicalBoundarySlab));
    WCNS_REQUIRE_THROWS(
        std::out_of_range,
        state.at({-1, 3, 0}, 0, AccessRegion::ConnectionHalo));
    WCNS_REQUIRE_THROWS(
        std::out_of_range,
        state.at({1, 1, -1}, 0, AccessRegion::PhysicalBoundarySlab));
    WCNS_REQUIRE_THROWS(
        std::out_of_range,
        state.at({-1, 1, 0}, 0, AccessRegion::Interior));

    const FieldDomain geometry_domain(
        {5, 4, 1},
        3,
        2,
        TopologyLocation::Vertex,
        {.connection_halo = true, .physical_boundary_slab = false});
    TopologyField<double> coordinates(geometry_domain, 3, 0.0);
    coordinates.at({5, 1, 0}, 0, AccessRegion::ConnectionHalo) = 1.0;
    WCNS_REQUIRE_THROWS(
        std::logic_error,
        coordinates.at({-1, 1, 0}, 0, AccessRegion::PhysicalBoundarySlab));

    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        FieldDomain({4, 3, 2}, 3, 2, TopologyLocation::Cell));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        FieldDomain({4, 3, 1}, 3, 2, TopologyLocation::FaceK));
}
