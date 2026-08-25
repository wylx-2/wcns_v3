#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>
#include <wcns/solver/flow_fields.hpp>

#include <array>
#include <stdexcept>

namespace wcns {

class PhysicsError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum EulerIndex : int {
    density = 0,
    momentum_x = 1,
    momentum_y = 2,
    momentum_z = 3,
    total_energy = 4,
};

enum PrimitiveIndex : int {
    primitive_density = 0,
    velocity_x = 1,
    velocity_y = 2,
    velocity_z = 3,
    pressure = 4,
};

using ConservativeState = std::array<Real, euler_components>;
using PrimitiveState = std::array<Real, euler_components>;

struct Normal3 {
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
};

struct IdealGas {
    Real gamma = 1.4;
    Real density_floor = 1.0e-12;
    Real pressure_floor = 1.0e-12;

    void validate() const;
};

[[nodiscard]] ConservativeState to_conservative(
    const PrimitiveState& primitive,
    const IdealGas& gas = {});

[[nodiscard]] PrimitiveState to_primitive(
    const ConservativeState& conservative,
    const IdealGas& gas = {});

[[nodiscard]] Real sound_speed(
    const PrimitiveState& primitive,
    const IdealGas& gas = {});

[[nodiscard]] ConservativeState euler_flux(
    const PrimitiveState& primitive,
    Normal3 unit_normal,
    const IdealGas& gas = {});

[[nodiscard]] ConservativeState rusanov_flux(
    const PrimitiveState& left,
    const PrimitiveState& right,
    Normal3 unit_normal,
    const IdealGas& gas = {});

[[nodiscard]] ConservativeState load_conservative(
    const Field<Real>& field,
    Index3 index);
[[nodiscard]] PrimitiveState load_primitive(const Field<Real>& field, Index3 index);
void store_state(Field<Real>& field, Index3 index, const ConservativeState& state);

} // namespace wcns

