#include "test_support.hpp"

#include <wcns/physics/source_terms.hpp>
#include <wcns/solver/spatial_operator.hpp>

#include <stdexcept>
#include <string>

// 验收阶段 H 源项配置契约、关闭路径和未实现模型拒绝路径。
void test_source_terms()
{
    using namespace wcns;

    const SourceTermConfig disabled;
    disabled.validate();
    const auto registry = SourceTermRegistry::create_stage_h(disabled);
    WCNS_REQUIRE(registry.empty());
    WCNS_REQUIRE(registry.size() == 0);
    WCNS_REQUIRE(disabled.summary() == "enable_source_terms=false;models=");
    WCNS_REQUIRE(disabled.restart_signature()
                 == "source_terms_v3;enable_source_terms=false;models=");

    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        (SourceTermConfig {false, {SourceModelKind::BodyForce}}.validate()));
    WCNS_REQUIRE_THROWS(std::invalid_argument, (SourceTermConfig {true, {}}.validate()));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (SourceTermConfig {true, {SourceModelKind::BodyForce, SourceModelKind::BodyForce}}
             .validate()));
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        (SourceTermConfig {true, {static_cast<SourceModelKind>(99)}}.validate()));

    const SourceTermConfig valid_but_unavailable {
        true,
        {SourceModelKind::BodyForce},
    };
    valid_but_unavailable.validate();
    WCNS_REQUIRE(
        valid_but_unavailable.summary().find("enable_source_terms=true;models=body_force;uniform=")
        == 0);
    WCNS_REQUIRE_THROWS(std::logic_error,
                        SourceTermRegistry::create_stage_h(valid_but_unavailable));

    SourceTermConfig pressure_gradient;
    pressure_gradient.enable_source_terms = true;
    pressure_gradient.models = {SourceModelKind::PressureGradient};
    pressure_gradient.pressure_gradient = {{0.8, -0.1, 0.2}};
    const auto pressure_registry = SourceTermRegistry::create_stage_j(pressure_gradient);
    const auto pressure_source
        = pressure_registry.evaluate({{2.0, 1.0, 4.0, -2.0, 10.0}}, {{0.0, 0.0, 0.0}}, 0.0, 3);
    WCNS_REQUIRE_NEAR(pressure_source[0], 0.0, 0.0);
    WCNS_REQUIRE_NEAR(pressure_source[1], 0.8, 0.0);
    WCNS_REQUIRE_NEAR(pressure_source[2], -0.1, 0.0);
    WCNS_REQUIRE_NEAR(pressure_source[3], 0.2, 0.0);
    WCNS_REQUIRE_NEAR(pressure_source[4], 0.0, 1.0e-15);

    SpatialParameters parameters;
    parameters.validate();
    parameters.source_terms = valid_but_unavailable;
    WCNS_REQUIRE_THROWS(std::logic_error, parameters.validate());
}
