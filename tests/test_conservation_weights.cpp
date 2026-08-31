#include "test_support.hpp"

#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <vector>

// 验收两套 profile 在规定网格尺寸上生成唯一正权重并满足 D^T w=eR-eL。
void test_line_conservation_weights()
{
    using namespace wcns;
    for (const auto kind : {
             AlgorithmProfileKind::PhengleiWcns,
             AlgorithmProfileKind::Scmm6Wcns}) {
        const auto profile = ProfileFactory::create(kind);
        const int minimum
            = kind == AlgorithmProfileKind::PhengleiWcns ? 4 : 5;
        for (int count = minimum; count <= 12; ++count) {
            const auto weights
                = build_line_conservation_weights(profile, count);
            WCNS_REQUIRE(weights.maximum_residual < 1.0e-11);
            WCNS_REQUIRE(std::all_of(
                weights.cell_weights.begin(),
                weights.cell_weights.end(),
                [](Real value) { return value > 0.0; }));

            std::vector<Real> faces(static_cast<std::size_t>(count + 1));
            for (int face = 0; face <= count; ++face) {
                faces[static_cast<std::size_t>(face)]
                    = 0.25 * face * face - 0.75 * face + 2.0;
            }
            const auto derivative
                = LineOperators::build(profile, count).differentiate(faces);
            Real integral = 0.0;
            for (int cell = 0; cell < count; ++cell) {
                integral += weights.cell_weights[static_cast<std::size_t>(cell)]
                    * derivative[static_cast<std::size_t>(cell)];
            }
            WCNS_REQUIRE_NEAR(integral, faces.back() - faces.front(), 1.0e-10);
        }
    }

    const auto ph = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto periodic = build_line_conservation_weights(ph, 3, true);
    WCNS_REQUIRE(
        periodic.cell_weights == (std::vector<Real> {1.0, 1.0, 1.0}));
    WCNS_REQUIRE_THROWS(
        ConservationWeightError,
        build_line_conservation_weights(ph, 0, true));
}
