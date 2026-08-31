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
    WCNS_REQUIRE(
        disabled.restart_signature()
        == "source_terms_v1;enable_source_terms=false;models=");

    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (SourceTermConfig {false, {SourceModelKind::BodyForce}}.validate()));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (SourceTermConfig {true, {}}.validate()));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (SourceTermConfig {
            true,
            {SourceModelKind::BodyForce, SourceModelKind::BodyForce}}
             .validate()));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        (SourceTermConfig {true, {static_cast<SourceModelKind>(99)}}.validate()));

    const SourceTermConfig valid_but_unavailable {
        true,
        {SourceModelKind::BodyForce},
    };
    valid_but_unavailable.validate();
    WCNS_REQUIRE(
        valid_but_unavailable.summary()
        == "enable_source_terms=true;models=body_force");
    WCNS_REQUIRE_THROWS(
        std::logic_error,
        SourceTermRegistry::create_stage_h(valid_but_unavailable));

    SpatialParameters parameters;
    parameters.validate();
    parameters.source_terms = valid_but_unavailable;
    WCNS_REQUIRE_THROWS(std::logic_error, parameters.validate());
}
