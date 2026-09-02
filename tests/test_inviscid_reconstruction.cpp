#include "test_support.hpp"

#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <string_view>
#include <vector>

namespace {

wcns::GasModel reconstruction_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales reconstruction_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive(
        {340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
}

class CustomAverageReconstruction final : public wcns::IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override
    {
        return "custom_average";
    }

    [[nodiscard]] wcns::StencilRequirement stencil_requirement() const noexcept override
    {
        return {};
    }

    [[nodiscard]] wcns::Real reconstruct_scalar(
        wcns::ScalarStencilView stencil,
        wcns::TraceSide side,
        const wcns::ReconstructionContext&) const override
    {
        if (stencil.size() != 6) throw std::invalid_argument("custom stencil size");
        return side == wcns::TraceSide::Left
            ? 0.5 * (stencil[2] + stencil[3])
            : 0.5 * (stencil[3] + stencil[2]);
    }
};

class CustomCentralRiemann final : public wcns::IRiemannSolver {
public:
    [[nodiscard]] std::string_view name() const noexcept override
    {
        return "custom_central";
    }

    [[nodiscard]] wcns::RiemannResult solve(
        const wcns::PressurePrimitiveState& left,
        const wcns::PressurePrimitiveState& right,
        wcns::Normal3 unit_normal,
        const wcns::GasModel& gas,
        const wcns::NumericalFloors& floors) const override
    {
        const wcns::IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto left_flux = wcns::euler_flux(left, unit_normal, ideal);
        const auto right_flux = wcns::euler_flux(right, unit_normal, ideal);
        wcns::RiemannResult result;
        for (int component = 0; component < wcns::euler_components; ++component) {
            result.flux_per_unit_area[static_cast<std::size_t>(component)]
                = 0.5 * (left_flux[static_cast<std::size_t>(component)]
                    + right_flux[static_cast<std::size_t>(component)]);
        }
        result.requested_solver = "custom_central";
        result.used_solver = "custom_central";
        return result;
    }
};

} // namespace

// 验收阶段 L 的重构/Riemann 注册表可增加自定义策略且拒绝未知或重复名称。
void test_stage_l_algorithm_registries()
{
    using namespace wcns;
    auto reconstruction_registry = ReconstructionRegistry::with_builtins();
    WCNS_REQUIRE(reconstruction_registry.names()
        == std::vector<std::string>({"linear5", "weno_js"}));
    reconstruction_registry.register_scheme(
        "custom_average",
        [] { return std::make_unique<CustomAverageReconstruction>(); });
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        reconstruction_registry.register_scheme(
            "custom_average",
            [] { return std::make_unique<CustomAverageReconstruction>(); }));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        reconstruction_registry.create("missing_reconstruction"));
    auto custom = reconstruction_registry.create("custom_average");
    const std::array<Real, 6> stencil {{0.0, 1.0, 2.0, 4.0, 8.0, 16.0}};
    WCNS_REQUIRE_NEAR(
        custom->reconstruct_scalar(stencil, TraceSide::Left, {}), 3.0, 0.0);
    const std::array<Real, 5> wrong_stencil {{0.0, 1.0, 2.0, 3.0, 4.0}};
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        reconstruction_registry.create("linear5")->reconstruct_scalar(
            wrong_stencil, TraceSide::Left, {}));
    ReconstructionConfig custom_config;
    custom_config.scheme = "custom_average";
    custom_config.validate(reconstruction_registry);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        custom_config.validate(ReconstructionRegistry::with_builtins()));

    auto riemann_registry = RiemannSolverRegistry::with_builtins();
    riemann_registry.register_solver(
        "custom_central",
        [] { return std::make_unique<CustomCentralRiemann>(); });
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        riemann_registry.register_solver(
            "custom_central",
            [] { return std::make_unique<CustomCentralRiemann>(); }));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        riemann_registry.create("missing_riemann"));
    const auto gas = reconstruction_gas();
    const PressurePrimitiveState state {1.0, 0.2, 0.0, 0.0, 1.0};
    const RiemannSolver custom_riemann("custom_central", riemann_registry);
    const auto result = custom_riemann.solve(
        state, state, {1.0, 0.0, 0.0}, gas, {});
    WCNS_REQUIRE(result.requested_solver == "custom_central");
    WCNS_REQUIRE(result.used_solver == "custom_central");
    WCNS_REQUIRE(custom_riemann.summary() == "riemann_solver=custom_central");
}

// 验收五阶线性重构的四次多项式精确性及 WCNS-JS 的尺度不变性。
void test_stage_j_scalar_reconstruction()
{
    using namespace wcns;
    std::array<Real, 6> polynomial {};
    for (int point = 0; point < 6; ++point) {
        const Real x = static_cast<Real>(point - 2);
        polynomial[static_cast<std::size_t>(point)] = x * x * x * x;
    }
    const auto linear = linear5_reconstruct(polynomial);
    WCNS_REQUIRE_NEAR(linear.left, 0.0625, 1.0e-14);
    WCNS_REQUIRE_NEAR(linear.right, 0.0625, 1.0e-14);

    const std::array<Real, 6> smooth {{0.7, 0.9, 1.2, 1.6, 2.1, 2.7}};
    std::array<Real, 6> scaled {};
    for (std::size_t index = 0; index < smooth.size(); ++index) {
        scaled[index] = 1.0e6 * smooth[index];
    }
    const auto base = wcns5_reconstruct_scaled(smooth, 1.0);
    const auto large = wcns5_reconstruct_scaled(scaled, 1.0e6);
    WCNS_REQUIRE_NEAR(large.left / 1.0e6, base.left, 2.0e-15);
    WCNS_REQUIRE_NEAR(large.right / 1.0e6, base.right, 2.0e-15);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        wcns5_reconstruct_scaled(smooth, 0.0));
    ReconstructionConfig inconsistent;
    inconsistent.nonlinear.epsilon = 1.0e-5;
    WCNS_REQUIRE_THROWS(std::invalid_argument, inconsistent.validate());
}

// 验收重构出现非物理状态时确定性回退到相邻单元一阶状态并记录计数。
void test_reconstruction_positivity_fallback()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const auto reference = reconstruction_reference(gas);
    Field<Real> conservative({4, 1, 1}, euler_components, 3);
    Field<Real> primitive({4, 1, 1}, euler_components, 3);
    const PressurePrimitiveState valid {1.0, 0.2, 0.0, 0.0, 1.0};
    for (int i = -3; i < 7; ++i) {
        store_state(primitive, {i, 0, 0}, valid);
        store_state(conservative, {i, 0, 0}, to_conservative(valid));
    }
    for (const int i : {0, 3}) {
        auto remote = valid;
        remote[0] = 100.0;
        remote[4] = 100.0;
        store_state(primitive, {i, 0, 0}, remote);
    }
    ReconstructionConfig config;
    config.scheme = std::string(reconstruction_name(ReconstructionKind::Linear5));
    ReconstructionDiagnostics diagnostics;
    const auto result = reconstruct_thermodynamic_face(
        conservative, primitive, Axis::I, {2, 0, 0}, config,
        gas, reference, diagnostics, 2);
    WCNS_REQUIRE(result.left == valid);
    WCNS_REQUIRE(result.right == valid);
    WCNS_REQUIRE(diagnostics.linear_faces == 1);
    WCNS_REQUIRE(diagnostics.first_order_fallbacks == 1);
}

// 验收 Rusanov 只接收单位法向、返回单位面积通量且自由流退化为解析 Euler 通量。
void test_rusanov_riemann_solver()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const NumericalFloors floors;
    const PressurePrimitiveState state {1.1, 0.7, -0.2, 0.0, 0.9};
    const Normal3 normal {0.6, 0.8, 0.0};
    const RiemannSolver solver;
    const auto flux = solver.flux(state, state, normal, gas, floors);
    const auto exact = euler_flux(
        state, normal, {gas.gamma(), floors.density, floors.pressure});
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            flux[static_cast<std::size_t>(component)],
            exact[static_cast<std::size_t>(component)], 1.0e-14);
    }
    WCNS_REQUIRE(solver.summary() == "riemann_solver=rusanov");
    WCNS_REQUIRE_THROWS(
        PhysicsError,
        solver.flux(state, state, {2.0, 0.0, 0.0}, gas, floors));
}
