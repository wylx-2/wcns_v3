#include "test_support.hpp"

#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
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
        == std::vector<std::string>({
            "linear5", "mdcd_hybrid", "mdcd_linear", "weno_js", "weno_z"}));
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
    WCNS_REQUIRE(riemann_registry.names()
        == std::vector<std::string>({"hllc", "roe", "rusanov"}));
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
    WCNS_REQUIRE(custom_riemann.summary().find(
        "riemann_solver=custom_central;") == 0);
}

// 验收 WENO-JS/WENO-Z/MDCD-LINEAR/MDCD-HYBRID 的常数保持、尺度及镜像约定。
void test_stage_l_scalar_reconstruction_schemes()
{
    using namespace wcns;
    const auto registry = ReconstructionRegistry::with_builtins();
    const std::array<Real, 6> constant {{2.5, 2.5, 2.5, 2.5, 2.5, 2.5}};
    for (const auto* name : {"weno_js", "weno_z", "mdcd_linear", "mdcd_hybrid"}) {
        const auto scheme = registry.create(name);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(constant, TraceSide::Left, {}), 2.5, 1.0e-14);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(constant, TraceSide::Right, {}), 2.5, 1.0e-14);
    }
    WCNS_REQUIRE_NEAR(mdcd_six_point_smoothness(constant, 2.5), 0.0, 0.0);

    const std::array<Real, 6> smooth {{0.7, 0.9, 1.2, 1.6, 2.1, 2.7}};
    auto reversed = smooth;
    std::reverse(reversed.begin(), reversed.end());
    for (const auto* name : {"weno_js", "weno_z", "mdcd_hybrid"}) {
        const auto scheme = registry.create(name);
        ReconstructionContext base;
        base.scale = 1.0;
        std::array<Real, 6> scaled {};
        for (std::size_t index = 0; index < smooth.size(); ++index) {
            scaled[index] = 1.0e6 * smooth[index];
        }
        ReconstructionContext large = base;
        large.scale = 1.0e6;
        const Real left = scheme->reconstruct_scalar(smooth, TraceSide::Left, base);
        const Real right = scheme->reconstruct_scalar(smooth, TraceSide::Right, base);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(scaled, TraceSide::Left, large) / 1.0e6,
            left, 3.0e-15);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(reversed, TraceSide::Left, base),
            right, 3.0e-15);
    }

    WcnsParameters parameters;
    const auto mdcd = mdcd_linear_reconstruct(smooth, parameters);
    const Real dispersion = parameters.mdcd_dispersion;
    const Real dissipation = parameters.mdcd_dissipation;
    const std::array<Real, 6> coefficients {{
        3.0 * (dispersion + dissipation) / 8.0,
        (-18.0 * dispersion - 30.0 * dissipation - 1.0) / 16.0,
        (12.0 * dispersion + 60.0 * dissipation + 9.0) / 16.0,
        (12.0 * dispersion - 60.0 * dissipation + 9.0) / 16.0,
        (-18.0 * dispersion + 30.0 * dissipation - 1.0) / 16.0,
        3.0 * (dispersion - dissipation) / 8.0,
    }};
    Real expected = 0.0;
    for (std::size_t index = 0; index < smooth.size(); ++index) {
        expected += coefficients[index] * smooth[index];
    }
    WCNS_REQUIRE_NEAR(mdcd.left, expected, 1.0e-15);

    const std::array<Real, 6> linear_data {{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};
    const auto hybrid_smooth = mdcd_hybrid_reconstruct_scaled(linear_data, 1.0, parameters);
    const auto linear_smooth = mdcd_linear_reconstruct(linear_data, parameters);
    WCNS_REQUIRE_NEAR(hybrid_smooth.left, linear_smooth.left, 0.0);
    WCNS_REQUIRE_NEAR(hybrid_smooth.right, linear_smooth.right, 0.0);

    const std::array<Real, 6> discontinuity {{1.0, 1.0, 1.0, 2.0, 2.0, 2.0}};
    WCNS_REQUIRE_NEAR(
        mdcd_six_point_smoothness(discontinuity, 1.0),
        279739.0 / 5040.0, 1.0e-13);
    const auto hybrid_jump = mdcd_hybrid_reconstruct_scaled(
        discontinuity, 1.0, parameters);
    WCNS_REQUIRE_NEAR(hybrid_jump.left, 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(hybrid_jump.right, 2.0, 1.0e-14);

    auto invalid = parameters;
    invalid.mdcd_dissipation = invalid.mdcd_dispersion;
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid.validate());
}

// 验收同一面 Roe 特征基的 L*R=I、守恒量往返变换和特征重构回退次序。
void test_stage_l_characteristic_reconstruction()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const auto reference = reconstruction_reference(gas);
    const NumericalFloors floors;
    const PressurePrimitiveState left_state {1.0, 0.6, -0.2, 0.0, 1.0};
    const PressurePrimitiveState right_state {0.8, -0.1, 0.3, 0.0, 0.7};
    const Normal3 normal {0.6, 0.8, 0.0};
    const auto basis = make_roe_characteristic_basis(
        left_state, right_state, normal, gas, floors, 2);
    for (int row = 0; row < euler_components; ++row) {
        for (int column = 0; column < euler_components; ++column) {
            Real product = 0.0;
            for (int inner = 0; inner < euler_components; ++inner) {
                product += basis.left[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(inner)]
                    * basis.right[static_cast<std::size_t>(inner)]
                        [static_cast<std::size_t>(column)];
            }
            WCNS_REQUIRE_NEAR(product, row == column ? 1.0 : 0.0, 2.0e-14);
        }
    }
    const auto conservative = to_conservative(left_state);
    const auto restored = restore_characteristic(
        project_characteristic(conservative, basis), basis);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            restored[static_cast<std::size_t>(component)],
            conservative[static_cast<std::size_t>(component)], 3.0e-14);
    }

    Field<Real> conservative_field({4, 1, 1}, euler_components, 3);
    Field<Real> primitive_field({4, 1, 1}, euler_components, 3);
    for (int i = -3; i < 7; ++i) {
        store_state(primitive_field, {i, 0, 0}, left_state);
        store_state(conservative_field, {i, 0, 0}, conservative);
    }
    ReconstructionConfig config;
    config.scheme = "weno_z";
    config.variables = ReconstructionVariables::Characteristic;
    ReconstructionDiagnostics diagnostics;
    const auto uniform = reconstruct_thermodynamic_face(
        conservative_field, primitive_field, Axis::I, {2, 0, 0}, config,
        gas, reference, diagnostics, 2, normal);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            uniform.left[static_cast<std::size_t>(component)],
            left_state[static_cast<std::size_t>(component)], 3.0e-14);
        WCNS_REQUIRE_NEAR(
            uniform.right[static_cast<std::size_t>(component)],
            left_state[static_cast<std::size_t>(component)], 3.0e-14);
    }
    WCNS_REQUIRE(diagnostics.characteristic_faces == 1);
    WCNS_REQUIRE(diagnostics.characteristic_fallbacks == 0);

    const PressurePrimitiveState valid {1.0, 0.2, 0.0, 0.0, 1.0};
    for (int i = -3; i < 7; ++i) {
        store_state(primitive_field, {i, 0, 0}, valid);
        store_state(conservative_field, {i, 0, 0}, to_conservative(valid));
    }
    for (const int i : {0, 3}) {
        auto remote = valid;
        remote[0] = 100.0;
        remote[4] = 100.0;
        store_state(primitive_field, {i, 0, 0}, remote);
        store_state(conservative_field, {i, 0, 0}, to_conservative(remote));
    }
    config.scheme = "linear5";
    diagnostics = {};
    const auto fallback = reconstruct_thermodynamic_face(
        conservative_field, primitive_field, Axis::I, {2, 0, 0}, config,
        gas, reference, diagnostics, 2, {1.0, 0.0, 0.0});
    WCNS_REQUIRE(fallback.left == valid);
    WCNS_REQUIRE(fallback.right == valid);
    WCNS_REQUIRE(diagnostics.characteristic_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.primitive_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.first_order_fallbacks == 1);
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

// 验收三种 Riemann 求解器的相容性、接触保持、迎风极限和法向反转对称性。
void test_stage_l_riemann_solvers()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const NumericalFloors floors;
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const Normal3 normal {0.6, 0.8, 0.0};
    const PressurePrimitiveState uniform {1.1, 0.7, -0.2, 0.0, 0.9};
    const auto exact = euler_flux(uniform, normal, ideal);
    for (const auto kind : {RiemannSolverKind::Rusanov,
             RiemannSolverKind::Hllc, RiemannSolverKind::Roe}) {
        const RiemannSolver solver(kind);
        const auto result = solver.solve(uniform, uniform, normal, gas, floors);
        WCNS_REQUIRE(result.requested_solver == solver.name());
        WCNS_REQUIRE(result.used_solver == solver.name());
        WCNS_REQUIRE(result.fallback_reason == RiemannFallbackReason::None);
        WCNS_REQUIRE(result.spectral_radius > 0.0);
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(
                result.flux_per_unit_area[static_cast<std::size_t>(component)],
                exact[static_cast<std::size_t>(component)], 2.0e-13);
        }
    }

    const PressurePrimitiveState contact_left {1.0, 0.0, 0.0, 0.0, 1.0};
    const PressurePrimitiveState contact_right {2.0, 0.0, 0.0, 0.0, 1.0};
    const ConservativeState contact_flux {{0.0, 1.0, 0.0, 0.0, 0.0}};
    const auto hllc_contact = RiemannSolver(RiemannSolverKind::Hllc).flux(
        contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            hllc_contact[static_cast<std::size_t>(component)],
            contact_flux[static_cast<std::size_t>(component)], 3.0e-13);
    }
    const auto rusanov_contact = RiemannSolver(RiemannSolverKind::Rusanov).flux(
        contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(std::abs(hllc_contact[0]) < std::abs(rusanov_contact[0]));

    RiemannSolverParameters strong_entropy_fix;
    strong_entropy_fix.entropy_fix_coefficient = 0.2;
    const auto roe_contact = RiemannSolver(
        RiemannSolverKind::Roe, strong_entropy_fix).flux(
        contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(roe_contact[0] < 0.0);

    const PressurePrimitiveState supersonic_left {1.0, 4.0, 0.1, 0.0, 1.0};
    const PressurePrimitiveState supersonic_right {0.8, 3.5, -0.2, 0.0, 0.7};
    const auto supersonic_exact = euler_flux(
        supersonic_left, {1.0, 0.0, 0.0}, ideal);
    const auto supersonic_hllc = RiemannSolver(RiemannSolverKind::Hllc).flux(
        supersonic_left, supersonic_right, {1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            supersonic_hllc[static_cast<std::size_t>(component)],
            supersonic_exact[static_cast<std::size_t>(component)], 0.0);
    }

    const PressurePrimitiveState left {1.0, 0.9, -0.3, 0.1, 1.2};
    const PressurePrimitiveState right {0.7, -0.2, 0.4, -0.1, 0.6};
    const Normal3 reverse_normal {-normal.x, -normal.y, -normal.z};
    for (const auto kind : {RiemannSolverKind::Rusanov,
             RiemannSolverKind::Hllc, RiemannSolverKind::Roe}) {
        const RiemannSolver solver(kind);
        const auto forward = solver.flux(left, right, normal, gas, floors);
        const auto reverse = solver.flux(
            right, left, reverse_normal, gas, floors);
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(
                forward[static_cast<std::size_t>(component)],
                -reverse[static_cast<std::size_t>(component)], 5.0e-13);
        }
    }

    RiemannConfig config;
    config.scheme = "hllc";
    config.validate(RiemannSolverRegistry::with_builtins(config.parameters));
    WCNS_REQUIRE(config.restart_signature().find("riemann_config_v1;") == 0);
    config.scheme = "missing";
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        config.validate(RiemannSolverRegistry::with_builtins(config.parameters)));
}

// 验收默认 Rusanov 只接收单位法向并保持阶段 J 的解析 Euler 通量基线。
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
    WCNS_REQUIRE(solver.summary().find("riemann_solver=rusanov;") == 0);
    WCNS_REQUIRE_THROWS(
        PhysicsError,
        solver.flux(state, state, {2.0, 0.0, 0.0}, gas, floors));
}
