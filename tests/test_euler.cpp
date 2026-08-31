#include "test_support.hpp"

#include <wcns/solver/euler.hpp>

// 验收 Euler 状态换算、物理通量、声速和 Rusanov 通量。
void test_euler()
{
    using namespace wcns;
    const IdealGas gas;
    const PrimitiveState primitive {1.2, 2.0, -0.5, 0.25, 1.1};
    const auto conservative = to_conservative(primitive, gas);
    const auto recovered = to_primitive(conservative, gas);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            recovered[static_cast<std::size_t>(component)],
            primitive[static_cast<std::size_t>(component)],
            1.0e-14);
    }
    WCNS_REQUIRE_NEAR(sound_speed(primitive, gas), std::sqrt(1.4 * 1.1 / 1.2), 1.0e-14);

    const auto flux = euler_flux(primitive, {1.0, 0.0, 0.0}, gas);
    WCNS_REQUIRE_NEAR(flux[density], 2.4, 1.0e-14);
    WCNS_REQUIRE_NEAR(flux[momentum_x], 5.9, 1.0e-14);
    WCNS_REQUIRE_NEAR(flux[momentum_y], -1.2, 1.0e-14);
    const auto identical = rusanov_flux(primitive, primitive, {2.0, 0.0, 0.0}, gas);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            identical[static_cast<std::size_t>(component)],
            flux[static_cast<std::size_t>(component)],
            1.0e-14);
    }

    WCNS_REQUIRE_THROWS(
        PhysicsError, to_conservative(PrimitiveState {0.0, 0.0, 0.0, 0.0, 1.0}, gas));
    WCNS_REQUIRE_THROWS(
        PhysicsError, to_primitive(ConservativeState {1.0, 0.0, 0.0, 0.0, -1.0}, gas));
    WCNS_REQUIRE_THROWS(PhysicsError, euler_flux(primitive, {0.0, 0.0, 0.0}, gas));
    WCNS_REQUIRE_THROWS(
        PhysicsError,
        to_conservative(primitive, IdealGas {1.0, 1.0e-12, 1.0e-12}));
}
