#include "test_support.hpp"

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/linear_operators.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

wcns::Real polynomial(wcns::Real x, int degree)
{
    return std::pow(x, degree);
}

} // namespace

// 验收两套算法 profile 只能由工厂成套创建且不可交叉拼装。
void test_algorithm_profile()
{
    using namespace wcns;
    const auto ph = ProfileFactory::from_string("phenglei_wcns");
    const auto scmm = ProfileFactory::from_string("scmm6_wcns");
    WCNS_REQUIRE(ph.name() == "phenglei_wcns");
    WCNS_REQUIRE(scmm.name() == "scmm6_wcns");
    WCNS_REQUIRE(
        ph.restart_signature()
        == "algorithm_profile_v1;name=phenglei_wcns");
    ph.require_compatible(ph.components());
    WCNS_REQUIRE_THROWS(ProfileError, ph.require_compatible(scmm.components()));

    auto mixed = ph.components();
    mixed.derivative = FluxDerivativeKind::ScmmD6D4;
    WCNS_REQUIRE_THROWS(ProfileError, ProfileFactory::validate_bundle(mixed));
    WCNS_REQUIRE_THROWS(
        ProfileError,
        ProfileFactory::from_string("phenglei_metric_scmm_flux"));
}

// 验收 I4/D4-D2、I6/D6-D4、顶点 I6 和 GridDelta 的解析精度及尺寸保护。
void test_geometry_line_operators()
{
    using namespace wcns;
    const auto ph = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto scmm = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);

    const auto ph_line = LineOperators::build(ph, 8);
    std::vector<Real> cubic(8);
    for (int cell = 0; cell < 8; ++cell) {
        cubic[static_cast<std::size_t>(cell)]
            = polynomial(static_cast<Real>(cell) + 0.5, 3);
    }
    const auto ph_faces = ph_line.interpolate(cubic);
    for (int face = 0; face <= 8; ++face) {
        WCNS_REQUIRE_NEAR(
            ph_faces[static_cast<std::size_t>(face)],
            polynomial(static_cast<Real>(face), 3),
            2.0e-12);
    }

    std::vector<Real> linear_faces(9);
    for (int face = 0; face <= 8; ++face) {
        linear_faces[static_cast<std::size_t>(face)] = 2.0 * face - 1.0;
    }
    for (const auto value : ph_line.differentiate(linear_faces)) {
        WCNS_REQUIRE_NEAR(value, 2.0, 1.0e-14);
    }

    const auto scmm_line = LineOperators::build(scmm, 9);
    std::vector<Real> quartic(9);
    for (int cell = 0; cell < 9; ++cell) {
        quartic[static_cast<std::size_t>(cell)]
            = polynomial(static_cast<Real>(cell) + 0.5, 4);
    }
    const auto scmm_faces = scmm_line.interpolate(quartic);
    for (int face = 0; face <= 9; ++face) {
        WCNS_REQUIRE_NEAR(
            scmm_faces[static_cast<std::size_t>(face)],
            polynomial(static_cast<Real>(face), 4),
            2.0e-10);
    }
    WCNS_REQUIRE_THROWS(ProfileError, ph_line.require_profile(scmm));
    WCNS_REQUIRE_THROWS(ProfileError, LineOperators::build(ph, 3));
    WCNS_REQUIRE_THROWS(ProfileError, LineOperators::build(scmm, 4));

    std::vector<Real> vertices(10);
    for (int vertex = 0; vertex < 10; ++vertex) {
        vertices[static_cast<std::size_t>(vertex)]
            = polynomial(static_cast<Real>(vertex), 5);
    }
    const auto centers = interpolate_vertices_to_centers_i6(vertices);
    for (int cell = 0; cell < 9; ++cell) {
        WCNS_REQUIRE_NEAR(
            centers[static_cast<std::size_t>(cell)],
            polynomial(static_cast<Real>(cell) + 0.5, 5),
            2.0e-9);
    }

    std::vector<Real> refined(13);
    for (int index = 0; index < 13; ++index) {
        const Real coordinate = 0.5 * index;
        refined[static_cast<std::size_t>(index)] = coordinate * coordinate;
    }
    const auto derivative = grid_delta(refined);
    for (int index = 0; index < 13; ++index) {
        WCNS_REQUIRE_NEAR(derivative[static_cast<std::size_t>(index)], index, 1.0e-13);
    }
    WCNS_REQUIRE_THROWS(ProfileError, grid_delta({0.0, 1.0, 2.0, 3.0}));
}
