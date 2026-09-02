#include "test_support.hpp"

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/solver/viscous_boundary.hpp>

#include <cmath>

namespace {

wcns::GasModel boundary_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales boundary_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive(
        {300.0, 1.0, 300.0, 1.0, 1.8e-5, {}, {}}, gas);
}

wcns::StructuredBlock make_boundary_block()
{
    using namespace wcns;
    StructuredBlock block(0, "viscous-wall", 0, 2, 2, {9, 9, 1}, 3);
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    const auto cells = block.cell_extent();
    block.boundaries.push_back({
        "wall", BoundaryType::NoSlipIsothermalWall,
        {Axis::I, Side::Lower},
        {{0, 0, 0}, {0, vertices.nj - 1, 0}},
        {{0, 0, 0}, {0, cells.nj - 1, 0}},
        {{0, 0, 0}, {0, cells.nj - 1, 0}}, {}});
    return block;
}

void fill_wall_interior(
    wcns::StructuredBlock& block,
    const wcns::GasModel& gas,
    const wcns::ReferenceScales& reference)
{
    using namespace wcns;
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const Real x = static_cast<Real>(i) + 0.5;
            const Real temperature = 1.0 + 0.5 * x;
            const Real pressure_value = 1.0;
            const Real rho = gas.gamma() * reference.mach() * reference.mach()
                * pressure_value / temperature;
            const TemperaturePrimitiveState state {{
                rho, x, 2.0 * x, 0.0, temperature}};
            const Index3 cell {i, j, 0};
            for (int component = 0; component < fluid_components; ++component) {
                block.flow.temperature_primitive(
                    i, j, 0, component) = state[static_cast<std::size_t>(component)];
            }
            store_state(block.flow.primitive, cell,
                pressure_primitive(state, gas, reference, {}, 2));
            store_state(block.flow.conservative, cell,
                thermodynamic_conservative(state, gas, reference, {}, 2));
        }
    }
}

} // namespace

// 验收 PH 四次及 SCMM6 六次壁面 Dirichlet 法向导数矩条件。
void test_wall_dirichlet_derivative()
{
    using namespace wcns;
    for (const auto kind : {AlgorithmProfileKind::PhengleiWcns,
             AlgorithmProfileKind::Scmm6Wcns}) {
        const auto profile = ProfileFactory::create(kind);
        const int maximum = kind == AlgorithmProfileKind::PhengleiWcns ? 4 : 6;
        for (int degree = 0; degree <= maximum; ++degree) {
            std::vector<Real> centers;
            for (int index = 0; index < maximum; ++index) {
                centers.push_back(std::pow(static_cast<Real>(index) + 0.5, degree));
            }
            const Real derivative = wall_dirichlet_computational_derivative(
                degree == 0 ? 1.0 : 0.0, centers, profile);
            WCNS_REQUIRE_NEAR(derivative, degree == 1 ? 1.0 : 0.0, 2.0e-12);
        }
    }
}

// 验收无滑移等温/绝热壁真实面速度、温度及法向导数始终执行强约束。
void test_viscous_boundary_trace()
{
    using namespace wcns;
    const auto gas = boundary_gas();
    const auto reference = boundary_reference(gas);
    const NumericalFloors floors;
    for (const auto kind : {AlgorithmProfileKind::PhengleiWcns,
             AlgorithmProfileKind::Scmm6Wcns}) {
        auto block = make_boundary_block();
        fill_wall_interior(block, gas, reference);
        const auto profile = ProfileFactory::create(kind);
        const auto metric = initialize_metric_field(block, profile).metric;
        const auto& patch = block.boundaries.front();
        const Index3 face {0, 3, 0};
        ViscousFaceTrace raw;
        raw.state = {{1.0, 0.25, 0.5, 0.0, 1.25}};
        raw.gradients[static_cast<int>(ViscousPrimitive::VelocityX)]
            = {{1.0, 0.0, 0.0}};
        raw.gradients[static_cast<int>(ViscousPrimitive::VelocityY)]
            = {{2.0, 0.0, 0.0}};
        raw.gradients[static_cast<int>(ViscousPrimitive::VelocityZ)]
            = {{0.0, 0.0, 0.0}};
        raw.gradients[static_cast<int>(ViscousPrimitive::Temperature)]
            = {{0.5, 0.2, 0.0}};
        BoundaryData data;
        data.wall_temperature = 1.0;
        const auto isothermal = apply_viscous_boundary_trace(
            block, metric, patch, face, raw, 1.0, {-1.0, 0.0, 0.0},
            data, profile, gas, reference, floors);
        WCNS_REQUIRE(isothermal.state[temperature_velocity_x] == 0.0);
        WCNS_REQUIRE(isothermal.state[temperature_velocity_y] == 0.0);
        WCNS_REQUIRE(isothermal.state[temperature_value] == 1.0);
        WCNS_REQUIRE_NEAR(
            isothermal.gradients[0][0], 1.0, 3.0e-12);
        WCNS_REQUIRE_NEAR(
            isothermal.gradients[1][0], 2.0, 3.0e-12);
        WCNS_REQUIRE_NEAR(
            isothermal.gradients[3][0], 0.5, 3.0e-12);
        WCNS_REQUIRE(isothermal.gradients[3][1] == 0.0);
        WCNS_REQUIRE_NEAR(
            isothermal.state[temperature_density],
            gas.gamma() * reference.mach() * reference.mach(), 1.0e-14);

        auto adiabatic_patch = patch;
        adiabatic_patch.type = BoundaryType::NoSlipAdiabaticWall;
        BoundaryData adiabatic_data;
        const auto adiabatic = apply_viscous_boundary_trace(
            block, metric, adiabatic_patch, face, raw, 1.0,
            {-1.0, 0.0, 0.0}, adiabatic_data, profile,
            gas, reference, floors);
        WCNS_REQUIRE(adiabatic.gradients[3][0] == 0.0);
        WCNS_REQUIRE_NEAR(adiabatic.gradients[3][1], 0.2, 1.0e-15);

        data.wall_velocity = {{1.0, 0.0, 0.0}};
        WCNS_REQUIRE_THROWS(
            PhysicsConfigurationError,
            apply_viscous_boundary_trace(
                block, metric, patch, face, raw, 1.0,
                {-1.0, 0.0, 0.0}, data, profile,
                gas, reference, floors));
    }
}
