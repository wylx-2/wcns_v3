#include "test_support.hpp"

#include <wcns/mesh/metrics.hpp>
#include <wcns/physics/double_mach_reflection.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/solver/physical_boundary.hpp>

#include <cmath>

namespace {

wcns::GasModel double_mach_gas(wcns::Real gamma = 1.4)
{
    wcns::GasModelInput input;
    input.gamma = gamma;
    input.specific_gas_constant = 1.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales double_mach_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive({1.0, 1.0, 1.0, 1.0, 1.0e-5, {}, {}}, gas);
}

wcns::StructuredBlock double_mach_top_block()
{
    using namespace wcns;
    StructuredBlock block(0, "double-mach", 0, 2, 2, {7, 6, 1}, 3);
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 7; ++i) {
            block.coordinates.x(i, j, 0) = 4.0 * static_cast<Real>(i) / 6.0;
            block.coordinates.y(i, j, 0) = static_cast<Real>(j) / 5.0;
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    compute_metrics(block);
    block.boundaries.push_back({"top",
                                BoundaryType::DoubleMachReflection,
                                {Axis::J, Side::Upper},
                                {{0, 5, 0}, {6, 5, 0}},
                                {{0, 4, 0}, {5, 4, 0}},
                                {{0, 5, 0}, {5, 5, 0}},
                                {}});
    return block;
}

} // namespace

// 验收经典双马赫反射的固定状态、移动激波、初场以及分段/时变边界。
void test_double_mach_reflection()
{
    using namespace wcns;
    const DoubleMachReflection model;
    const auto gas = double_mach_gas();
    const auto reference = double_mach_reference(gas);
    const NumericalFloors floors;
    model.validate(gas.gamma(), 2);
    WCNS_REQUIRE_THROWS(PhysicsConfigurationError, model.validate(double_mach_gas(1.3).gamma(), 2));
    WCNS_REQUIRE_THROWS(PhysicsConfigurationError, model.validate(gas.gamma(), 3));

    constexpr Real sqrt_three = 1.732050807568877293527446341505872367;
    WCNS_REQUIRE_NEAR(model.shock_x(1.0, 0.2), 1.0 / 6.0 + 5.0 / sqrt_three, 1.0e-14);
    const auto upstream = model.upstream_state();
    const auto downstream = model.post_shock_state();
    WCNS_REQUIRE_NEAR(upstream[0], 1.4, 0.0);
    WCNS_REQUIRE_NEAR(upstream[4], 1.0, 0.0);
    WCNS_REQUIRE_NEAR(downstream[0], 8.0, 0.0);
    WCNS_REQUIRE_NEAR(downstream[1], 8.25 * sqrt_three / 2.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(downstream[2], -4.125, 0.0);
    WCNS_REQUIRE_NEAR(downstream[4], 116.5, 0.0);

    InitialConditionConfig initial;
    initial.type = "double_mach_reflection";
    const auto initialized_post = pressure_primitive(
        FlowInitializer::evaluate(initial, {0.1, 0.0, 0.0}, gas, reference, floors, 2),
        gas,
        reference,
        floors,
        2);
    const auto initialized_upstream = pressure_primitive(
        FlowInitializer::evaluate(initial, {3.0, 1.0, 0.0}, gas, reference, floors, 2),
        gas,
        reference,
        floors,
        2);
    WCNS_REQUIRE(initialized_post == downstream);
    WCNS_REQUIRE(initialized_upstream == upstream);

    BoundaryData data;
    data.double_mach_reflection = model;
    const PressurePrimitiveState interior {1.4, 1.0, 2.0, 0.0, 1.0};
    const PressurePrimitiveState reconstructed {1.2, -1.0, -3.0, 0.0, 0.8};
    BoundaryPatch bottom {
        "bottom", BoundaryType::DoubleMachReflection, {Axis::J, Side::Lower}, {}, {}, {}, {}};
    const auto bottom_inflow = apply_inviscid_boundary_face_state(bottom,
                                                                  interior,
                                                                  reconstructed,
                                                                  {0.0, -1.0, 0.0},
                                                                  data,
                                                                  {},
                                                                  gas,
                                                                  reference,
                                                                  floors,
                                                                  2,
                                                                  {0.1, 0.0, 0.0},
                                                                  0.0);
    WCNS_REQUIRE(bottom_inflow == downstream);
    const auto bottom_wall = apply_inviscid_boundary_face_state(bottom,
                                                                interior,
                                                                reconstructed,
                                                                {0.0, -1.0, 0.0},
                                                                data,
                                                                {},
                                                                gas,
                                                                reference,
                                                                floors,
                                                                2,
                                                                {0.5, 0.0, 0.0},
                                                                0.0);
    WCNS_REQUIRE_NEAR(bottom_wall[1], interior[1], 0.0);
    WCNS_REQUIRE_NEAR(bottom_wall[2], -interior[2], 0.0);

    BoundaryPatch top = bottom;
    top.face = {Axis::J, Side::Upper};
    const auto top_initial = apply_inviscid_boundary_face_state(top,
                                                                interior,
                                                                reconstructed,
                                                                {0.0, 1.0, 0.0},
                                                                data,
                                                                {},
                                                                gas,
                                                                reference,
                                                                floors,
                                                                2,
                                                                {1.0, 1.0, 0.0},
                                                                0.0);
    const auto top_later = apply_inviscid_boundary_face_state(top,
                                                              interior,
                                                              reconstructed,
                                                              {0.0, 1.0, 0.0},
                                                              data,
                                                              {},
                                                              gas,
                                                              reference,
                                                              floors,
                                                              2,
                                                              {1.0, 1.0, 0.0},
                                                              0.1);
    WCNS_REQUIRE(top_initial == upstream);
    WCNS_REQUIRE(top_later == downstream);

    auto block = double_mach_top_block();
    const auto upstream_temperature = temperature_primitive(upstream, gas, reference, floors, 2);
    const auto conservative
        = thermodynamic_conservative(upstream_temperature, gas, reference, floors, 2);
    for (int j = 0; j < block.cell_extent().nj; ++j) {
        for (int i = 0; i < block.cell_extent().ni; ++i) {
            store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
    update_temperature_primitive_interior(block, gas, reference, floors);
    BoundaryDataMap boundary_data {{"top", data}};
    const auto first_fill
        = PhysicalGhostStateOperator::fill(block, boundary_data, gas, reference, floors, 1, 0.0);
    WCNS_REQUIRE(first_fill.version == 1);
    WCNS_REQUIRE(load_primitive(block.flow.primitive, {0, 5, 0}) == downstream);
    WCNS_REQUIRE(load_primitive(block.flow.primitive, {5, 5, 0}) == upstream);
    const auto second_fill
        = PhysicalGhostStateOperator::fill(block, boundary_data, gas, reference, floors, 2, 0.3);
    WCNS_REQUIRE(second_fill.version == 2);
    WCNS_REQUIRE(load_primitive(block.flow.primitive, {5, 5, 0}) == downstream);
}
