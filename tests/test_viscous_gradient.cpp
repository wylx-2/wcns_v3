#include "test_support.hpp"

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/solver/viscous_gradient.hpp>
#include <wcns/solver/viscous_halo.hpp>
#include <wcns/solver/viscous_operator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

void fill_identity_coordinates(wcns::StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<wcns::Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<wcns::Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
}

void fill_linear_temperature_state(wcns::StructuredBlock& block)
{
    const auto cells = block.cell_extent();
    auto& state = block.flow.temperature_primitive;
    const int ghost = state.ghost_width();
    const auto store = [&](int i, int j) {
        const wcns::Real x = static_cast<wcns::Real>(i) + 0.5;
        const wcns::Real y = static_cast<wcns::Real>(j) + 0.5;
        state(i, j, 0, wcns::temperature_density) = 1.0;
        state(i, j, 0, wcns::temperature_velocity_x) = 1.0 + x + 2.0 * y;
        state(i, j, 0, wcns::temperature_velocity_y) = -2.0 - 0.5 * x + 0.25 * y;
        state(i, j, 0, wcns::temperature_velocity_z) = 0.0;
        state(i, j, 0, wcns::temperature_value) = 3.0 + 0.1 * x + 0.2 * y;
    };
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i)
            store(i, j);
    }
    for (int layer = 1; layer <= ghost; ++layer) {
        for (int j = 0; j < cells.nj; ++j) {
            store(-layer, j);
            store(cells.ni - 1 + layer, j);
        }
        for (int i = 0; i < cells.ni; ++i) {
            store(i, -layer);
            store(i, cells.nj - 1 + layer);
        }
    }
}

void run_linear_gradient(wcns::AlgorithmProfileKind kind)
{
    using namespace wcns;
    StructuredBlock block(0, "linear-gradient", 0, 2, 2, {9, 9, 1}, 3);
    fill_identity_coordinates(block);
    fill_linear_temperature_state(block);
    const auto profile = ProfileFactory::create(kind);
    const auto metric = initialize_metric_field(block, profile).metric;
    const auto operands = compute_gradient_face_operands(block, metric, profile, 7);
    const auto gradients = compute_primitive_gradients(block, metric, operands, profile);
    WCNS_REQUIRE(gradients.version() == 7);
    WCNS_REQUIRE(gradients.halo_layers() == (kind == AlgorithmProfileKind::PhengleiWcns ? 2 : 3));
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const Index3 cell {i, j, 0};
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::VelocityX, 0), 1.0, 2.0e-13);
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::VelocityX, 1), 2.0, 2.0e-13);
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::VelocityY, 0), -0.5, 2.0e-13);
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::VelocityY, 1), 0.25, 2.0e-13);
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::Temperature, 0), 0.1, 2.0e-13);
            WCNS_REQUIRE_NEAR(gradients(cell, ViscousPrimitive::Temperature, 1), 0.2, 2.0e-13);
            for (int variable = 0; variable < viscous_primitive_components; ++variable) {
                WCNS_REQUIRE(gradients(cell, static_cast<ViscousPrimitive>(variable), 2) == 0.0);
            }
            for (int direction = 0; direction < 3; ++direction) {
                WCNS_REQUIRE(gradients(cell, ViscousPrimitive::VelocityZ, direction) == 0.0);
            }
        }
    }
    WCNS_REQUIRE(std::isnan(gradients.values()(-1, -1, 0, 0)));
}

void run_linear_gradient_3d(wcns::AlgorithmProfileKind kind)
{
    using namespace wcns;
    StructuredBlock block(0, "linear-gradient-3d", 0, 3, 3, {9, 9, 9}, 3);
    const auto vertices = block.vertex_extent();
    for (int k = 0; k < vertices.nk; ++k) {
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                block.coordinates.x(i, j, k) = static_cast<Real>(i);
                block.coordinates.y(i, j, k) = static_cast<Real>(j);
                block.coordinates.z(i, j, k) = static_cast<Real>(k);
            }
        }
    }
    const auto cells = block.cell_extent();
    auto& state = block.flow.temperature_primitive;
    const auto store = [&](int i, int j, int k) {
        const Real x = static_cast<Real>(i) + 0.5;
        const Real y = static_cast<Real>(j) + 0.5;
        const Real z = static_cast<Real>(k) + 0.5;
        state(i, j, k, temperature_density) = 1.0;
        state(i, j, k, temperature_velocity_x) = 1.0 + x + 2.0 * y + 3.0 * z;
        state(i, j, k, temperature_velocity_y) = -2.0 - 0.5 * x + 0.25 * y - 0.75 * z;
        state(i, j, k, temperature_velocity_z) = 0.5 + 0.4 * x - 0.3 * y + 0.2 * z;
        state(i, j, k, temperature_value) = 3.0 + 0.1 * x + 0.2 * y + 0.3 * z;
    };
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i)
                store(i, j, k);
        }
    }
    for (int layer = 1; layer <= state.ghost_width(); ++layer) {
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                store(-layer, j, k);
                store(cells.ni - 1 + layer, j, k);
            }
        }
        for (int k = 0; k < cells.nk; ++k) {
            for (int i = 0; i < cells.ni; ++i) {
                store(i, -layer, k);
                store(i, cells.nj - 1 + layer, k);
            }
        }
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                store(i, j, -layer);
                store(i, j, cells.nk - 1 + layer);
            }
        }
    }
    const auto profile = ProfileFactory::create(kind);
    const auto metric = initialize_metric_field(block, profile).metric;
    const auto operands = compute_gradient_face_operands(block, metric, profile, 9);
    const auto gradients = compute_primitive_gradients(block, metric, operands, profile);
    const std::array<std::array<Real, 3>, 4> exact {{
        {{1.0, 2.0, 3.0}},
        {{-0.5, 0.25, -0.75}},
        {{0.4, -0.3, 0.2}},
        {{0.1, 0.2, 0.3}},
    }};
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                for (int variable = 0; variable < viscous_primitive_components; ++variable) {
                    for (int direction = 0; direction < 3; ++direction) {
                        WCNS_REQUIRE_NEAR(gradients({i, j, k},
                                                    static_cast<ViscousPrimitive>(variable),
                                                    direction),
                                          exact[static_cast<std::size_t>(variable)]
                                               [static_cast<std::size_t>(direction)],
                                          2.0e-12);
                    }
                }
            }
        }
    }
    WCNS_REQUIRE(std::isnan(gradients.values()(-1, -1, -1, 0)));
}

struct ManufacturedErrorNorms {
    wcns::Real interior_l1 = 0.0;
    wcns::Real interior_l2 = 0.0;
    wcns::Real interior_linf = 0.0;
    wcns::Real boundary_l1 = 0.0;
    wcns::Real boundary_l2 = 0.0;
    wcns::Real boundary_linf = 0.0;
    wcns::Real global_l1 = 0.0;
    wcns::Real global_l2 = 0.0;
    wcns::Real global_linf = 0.0;
};

ManufacturedErrorNorms manufactured_temperature_error(wcns::AlgorithmProfileKind kind,
                                                      int cell_count)
{
    using namespace wcns;
    StructuredBlock block(
        0, "viscous-manufactured", 0, 2, 2, {cell_count + 1, cell_count + 1, 1}, 3);
    const Real spacing = 1.0 / static_cast<Real>(cell_count);
    const auto vertices = block.vertex_extent();
    const auto cells = block.cell_extent();
    block.boundaries.push_back({"lower",
                                BoundaryType::NoSlipIsothermalWall,
                                {Axis::J, Side::Lower},
                                {{0, 0, 0}, {vertices.ni - 1, 0, 0}},
                                {{0, 0, 0}, {cells.ni - 1, 0, 0}},
                                {{0, 0, 0}, {cells.ni - 1, 0, 0}},
                                {}});
    block.boundaries.push_back({"upper",
                                BoundaryType::NoSlipIsothermalWall,
                                {Axis::J, Side::Upper},
                                {{0, vertices.nj - 1, 0}, {vertices.ni - 1, vertices.nj - 1, 0}},
                                {{0, cells.nj - 1, 0}, {cells.ni - 1, cells.nj - 1, 0}},
                                {{0, cells.nj, 0}, {cells.ni - 1, cells.nj, 0}},
                                {}});
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i) * spacing;
            block.coordinates.y(i, j, 0) = static_cast<Real>(j) * spacing;
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    constexpr Real amplitude = 0.05;
    const Real wave_number = 2.0 * std::acos(-1.0);
    const auto temperature = [&](int j) {
        const Real y = (static_cast<Real>(j) + 0.5) * spacing;
        return 1.0 + amplitude * std::sin(wave_number * y);
    };
    auto& state = block.flow.temperature_primitive;
    const auto store = [&](int i, int j) {
        const Real value = temperature(j);
        state(i, j, 0, temperature_density) = 1.0;
        state(i, j, 0, temperature_velocity_x) = 0.0;
        state(i, j, 0, temperature_velocity_y) = 0.0;
        state(i, j, 0, temperature_velocity_z) = 0.0;
        state(i, j, 0, temperature_value) = value;
        block.flow.primitive(i, j, 0, density) = 1.0;
        block.flow.primitive(i, j, 0, velocity_x) = 0.0;
        block.flow.primitive(i, j, 0, velocity_y) = 0.0;
        block.flow.primitive(i, j, 0, velocity_z) = 0.0;
        block.flow.primitive(i, j, 0, pressure) = value / 1.4;
    };
    for (int j = 0; j < cell_count; ++j) {
        for (int i = 0; i < cell_count; ++i)
            store(i, j);
    }
    for (int layer = 1; layer <= state.ghost_width(); ++layer) {
        for (int j = 0; j < cell_count; ++j) {
            store(-layer, j);
            store(cell_count - 1 + layer, j);
        }
        for (int i = 0; i < cell_count; ++i) {
            store(i, -layer);
            store(i, cell_count - 1 + layer);
        }
    }

    const auto profile = ProfileFactory::create(kind);
    const auto metric = initialize_metric_field(block, profile).metric;
    const auto operands = compute_gradient_face_operands(block, metric, profile, 17);
    const auto gradients = compute_primitive_gradients(block, metric, operands, profile);
    GasModelInput gas_input;
    gas_input.specific_gas_constant = 1.0 / gas_input.gamma;
    const auto gas = GasModel::from_input(gas_input);
    const auto reference = ReferenceScales::derive({1.0, 1.0, 1.0, 1.0, 1.0, {}, {}}, gas);
    const TransportConfig transport_config;
    const TransportModel transport(transport_config);
    const NumericalFloors floors;
    BoundaryData wall;
    wall.wall_temperature = 1.0;
    const BoundaryDataMap boundary_data {{"lower", wall}, {"upper", wall}};
    const auto flux = compute_viscous_face_fluxes(
        block, metric, gradients, profile, transport, boundary_data, gas, reference, floors, 17);
    block.flow.residual.fill(0.0);
    add_wcns_viscous_residual(block, metric, flux, profile, reference.reynolds());

    const Real chi = transport.thermal_coefficient(1.0, gas, reference);
    const int margin = kind == AlgorithmProfileKind::PhengleiWcns ? 5 : 7;
    struct Accumulator {
        Real absolute = 0.0;
        Real squared = 0.0;
        Real maximum = 0.0;
        int count = 0;
        void add(Real error)
        {
            WCNS_REQUIRE(std::isfinite(error));
            absolute += std::abs(error);
            squared += error * error;
            maximum = std::max(maximum, std::abs(error));
            ++count;
        }
    } interior, boundary, global;
    for (int j = 0; j < cell_count; ++j) {
        const Real y = (static_cast<Real>(j) + 0.5) * spacing;
        const Real exact = -chi * amplitude * wave_number * wave_number * std::sin(wave_number * y)
            / reference.reynolds();
        for (int i = 0; i < cell_count; ++i) {
            const Real error = block.flow.residual(i, j, 0, total_energy) - exact;
            global.add(error);
            if (j >= margin && j < cell_count - margin)
                interior.add(error);
            else
                boundary.add(error);
        }
    }
    WCNS_REQUIRE(interior.count > 0 && boundary.count > 0 && global.count > 0);
    const auto l1 = [](const Accumulator& values) {
        return values.absolute / static_cast<Real>(values.count);
    };
    const auto l2 = [](const Accumulator& values) {
        return std::sqrt(values.squared / static_cast<Real>(values.count));
    };
    return {
        l1(interior),
        l2(interior),
        interior.maximum,
        l1(boundary),
        l2(boundary),
        boundary.maximum,
        l1(global),
        l2(global),
        global.maximum,
    };
}

ManufacturedErrorNorms manufactured_temperature_error_3d(wcns::AlgorithmProfileKind kind,
                                                         int cell_count)
{
    using namespace wcns;
    constexpr int transverse_cells = 8;
    StructuredBlock block(0,
                          "viscous-manufactured-3d",
                          0,
                          3,
                          3,
                          {transverse_cells + 1, cell_count + 1, transverse_cells + 1},
                          3);
    const Real spacing = 1.0 / static_cast<Real>(cell_count);
    const auto vertices = block.vertex_extent();
    const auto cells = block.cell_extent();
    for (int k = 0; k < vertices.nk; ++k) {
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                block.coordinates.x(i, j, k) = static_cast<Real>(i) / transverse_cells;
                block.coordinates.y(i, j, k) = static_cast<Real>(j) * spacing;
                block.coordinates.z(i, j, k) = static_cast<Real>(k) / transverse_cells;
            }
        }
    }
    block.boundaries.push_back({"lower",
                                BoundaryType::NoSlipIsothermalWall,
                                {Axis::J, Side::Lower},
                                {{0, 0, 0}, {vertices.ni - 1, 0, vertices.nk - 1}},
                                {{0, 0, 0}, {cells.ni - 1, 0, cells.nk - 1}},
                                {{0, 0, 0}, {cells.ni - 1, 0, cells.nk - 1}},
                                {}});
    block.boundaries.push_back(
        {"upper",
         BoundaryType::NoSlipIsothermalWall,
         {Axis::J, Side::Upper},
         {{0, vertices.nj - 1, 0}, {vertices.ni - 1, vertices.nj - 1, vertices.nk - 1}},
         {{0, cells.nj - 1, 0}, {cells.ni - 1, cells.nj - 1, cells.nk - 1}},
         {{0, cells.nj, 0}, {cells.ni - 1, cells.nj, cells.nk - 1}},
         {}});
    constexpr Real amplitude = 0.05;
    const Real wave_number = 2.0 * std::acos(-1.0);
    const auto temperature = [&](int j) {
        const Real y = (static_cast<Real>(j) + 0.5) * spacing;
        return 1.0 + amplitude * std::sin(wave_number * y);
    };
    auto& state = block.flow.temperature_primitive;
    const auto store = [&](int i, int j, int k) {
        const Real value = temperature(j);
        state(i, j, k, temperature_density) = 1.0;
        state(i, j, k, temperature_velocity_x) = 0.0;
        state(i, j, k, temperature_velocity_y) = 0.0;
        state(i, j, k, temperature_velocity_z) = 0.0;
        state(i, j, k, temperature_value) = value;
        block.flow.primitive(i, j, k, density) = 1.0;
        block.flow.primitive(i, j, k, velocity_x) = 0.0;
        block.flow.primitive(i, j, k, velocity_y) = 0.0;
        block.flow.primitive(i, j, k, velocity_z) = 0.0;
        block.flow.primitive(i, j, k, pressure) = value / 1.4;
    };
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i)
                store(i, j, k);
        }
    }
    for (int layer = 1; layer <= state.ghost_width(); ++layer) {
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                store(-layer, j, k);
                store(cells.ni - 1 + layer, j, k);
            }
        }
        for (int k = 0; k < cells.nk; ++k) {
            for (int i = 0; i < cells.ni; ++i) {
                store(i, -layer, k);
                store(i, cells.nj - 1 + layer, k);
            }
        }
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                store(i, j, -layer);
                store(i, j, cells.nk - 1 + layer);
            }
        }
    }
    const auto profile = ProfileFactory::create(kind);
    const auto metric = initialize_metric_field(block, profile).metric;
    const auto operands = compute_gradient_face_operands(block, metric, profile, 18);
    const auto gradients = compute_primitive_gradients(block, metric, operands, profile);
    GasModelInput gas_input;
    gas_input.specific_gas_constant = 1.0 / gas_input.gamma;
    const auto gas = GasModel::from_input(gas_input);
    const auto reference = ReferenceScales::derive({1.0, 1.0, 1.0, 1.0, 1.0, {}, {}}, gas);
    const TransportModel transport(TransportConfig {});
    BoundaryData wall;
    wall.wall_temperature = 1.0;
    const BoundaryDataMap boundary_data {{"lower", wall}, {"upper", wall}};
    const auto flux = compute_viscous_face_fluxes(
        block, metric, gradients, profile, transport, boundary_data, gas, reference, {}, 18);
    block.flow.residual.fill(0.0);
    add_wcns_viscous_residual(block, metric, flux, profile, reference.reynolds());

    const Real chi = transport.thermal_coefficient(1.0, gas, reference);
    const int margin = kind == AlgorithmProfileKind::PhengleiWcns ? 5 : 7;
    struct Accumulator {
        Real absolute = 0.0;
        Real squared = 0.0;
        Real maximum = 0.0;
        int count = 0;
        void add(Real error)
        {
            WCNS_REQUIRE(std::isfinite(error));
            absolute += std::abs(error);
            squared += error * error;
            maximum = std::max(maximum, std::abs(error));
            ++count;
        }
    } interior, boundary, global;
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            const Real y = (static_cast<Real>(j) + 0.5) * spacing;
            const Real exact = -chi * amplitude * wave_number * wave_number
                * std::sin(wave_number * y) / reference.reynolds();
            for (int i = 0; i < cells.ni; ++i) {
                const Real error = block.flow.residual(i, j, k, total_energy) - exact;
                global.add(error);
                if (j >= margin && j < cells.nj - margin)
                    interior.add(error);
                else
                    boundary.add(error);
            }
        }
    }
    WCNS_REQUIRE(interior.count > 0 && boundary.count > 0 && global.count > 0);
    const auto l1 = [](const Accumulator& values) {
        return values.absolute / static_cast<Real>(values.count);
    };
    const auto l2 = [](const Accumulator& values) {
        return std::sqrt(values.squared / static_cast<Real>(values.count));
    };
    return {
        l1(interior),
        l2(interior),
        interior.maximum,
        l1(boundary),
        l2(boundary),
        boundary.maximum,
        l1(global),
        l2(global),
        global.maximum,
    };
}

} // namespace

// 验收两套 profile 在二维仿射网格上以守恒形式精确恢复线性 primitive 梯度。
void test_viscous_linear_gradients()
{
    run_linear_gradient(wcns::AlgorithmProfileKind::PhengleiWcns);
    run_linear_gradient(wcns::AlgorithmProfileKind::Scmm6Wcns);
    run_linear_gradient_3d(wcns::AlgorithmProfileKind::PhengleiWcns);
    run_linear_gradient_3d(wcns::AlgorithmProfileKind::Scmm6Wcns);
}

// 验收旋转周期的速度梯度二阶张量和温度梯度矢量接收侧变换。
void test_viscous_gradient_periodic_transform()
{
    using namespace wcns;
    GradientExchangeDescriptor descriptor;
    descriptor.connection = 0;
    descriptor.receiver_block = 0;
    descriptor.donor_block = 1;
    descriptor.dimension = 3;
    descriptor.periodic.rotation = {{{{0.0, -1.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
    PrimitiveGradients donor {{
        {{1.0, 2.0, 3.0}},
        {{4.0, 5.0, 6.0}},
        {{7.0, 8.0, 9.0}},
        {{10.0, 11.0, 12.0}},
    }};
    const auto received = transform_primitive_gradients_for_receiver(donor, descriptor);
    WCNS_REQUIRE_NEAR(received[0][0], 5.0, 0.0);
    WCNS_REQUIRE_NEAR(received[0][1], -4.0, 0.0);
    WCNS_REQUIRE_NEAR(received[1][0], -2.0, 0.0);
    WCNS_REQUIRE_NEAR(received[1][1], 1.0, 0.0);
    WCNS_REQUIRE_NEAR(received[3][0], 11.0, 0.0);
    WCNS_REQUIRE_NEAR(received[3][1], -10.0, 0.0);

    FaceFluxExchangeDescriptor face;
    face.connection = 0;
    face.receiver_block = 0;
    face.donor_block = 1;
    face.orientation = -1.0;
    face.periodic = descriptor.periodic;
    GradientOperandState operand {};
    for (int i = 0; i < gradient_operand_components; ++i) {
        operand[static_cast<std::size_t>(i)] = static_cast<Real>(i + 1);
    }
    const auto transformed = transform_gradient_operand_for_receiver(operand, face, 3);
    WCNS_REQUIRE_NEAR(transformed[0], -5.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[1], 4.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[3], 2.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[4], -1.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[9], -11.0, 0.0);
    WCNS_REQUIRE_NEAR(transformed[10], 10.0, 0.0);
}

// 验收正弦温度黏性制造解在网格加密后达到 profile 对应的误差下降门槛。
void test_viscous_manufactured_convergence()
{
    using namespace wcns;
    for (const int dimension : {2, 3}) {
        for (const auto kind :
             {AlgorithmProfileKind::PhengleiWcns, AlgorithmProfileKind::Scmm6Wcns}) {
            std::vector<std::pair<int, ManufacturedErrorNorms>> results;
            for (const int cells : {24, 48, 96}) {
                results.push_back({cells,
                                   dimension == 2
                                       ? manufactured_temperature_error(kind, cells)
                                       : manufactured_temperature_error_3d(kind, cells)});
            }
            const char* name
                = kind == AlgorithmProfileKind::PhengleiWcns ? "phenglei_wcns" : "scmm6_wcns";
            std::cout << std::setprecision(17);
            for (const auto& [cells, error] : results) {
                std::cout << "stage_s_spatial dimension=" << dimension << " profile=" << name
                          << " cells=" << cells << " interior_l1=" << error.interior_l1
                          << " interior_l2=" << error.interior_l2
                          << " interior_linf=" << error.interior_linf
                          << " boundary_l1=" << error.boundary_l1
                          << " boundary_l2=" << error.boundary_l2
                          << " boundary_linf=" << error.boundary_linf
                          << " global_l1=" << error.global_l1 << " global_l2=" << error.global_l2
                          << " global_linf=" << error.global_linf << '\n';
            }
            for (std::size_t level = 1; level < results.size(); ++level) {
                const auto& coarse = results[level - 1].second;
                const auto& fine = results[level].second;
                const Real interior_order
                    = std::log(coarse.interior_l2 / fine.interior_l2) / std::log(2.0);
                const Real boundary_order
                    = std::log(coarse.boundary_l2 / fine.boundary_l2) / std::log(2.0);
                const Real global_order
                    = std::log(coarse.global_l2 / fine.global_l2) / std::log(2.0);
                std::cout << "stage_s_spatial_order dimension=" << dimension << " profile=" << name
                          << " coarse=" << results[level - 1].first
                          << " fine=" << results[level].first
                          << " interior_l2_order=" << interior_order
                          << " boundary_l2_order=" << boundary_order
                          << " global_l2_order=" << global_order << '\n';
                if (kind == AlgorithmProfileKind::PhengleiWcns) {
                    WCNS_REQUIRE(interior_order >= 3.5);
                    WCNS_REQUIRE(boundary_order >= 0.55);
                    WCNS_REQUIRE(global_order >= 1.0);
                } else {
                    WCNS_REQUIRE(interior_order >= 5.0);
                    WCNS_REQUIRE(boundary_order >= 3.0);
                    WCNS_REQUIRE(global_order >= 3.0);
                }
            }
        }
    }
}
