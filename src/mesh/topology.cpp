#include <wcns/mesh/topology.hpp>

#include <algorithm>
#include <cmath>

namespace wcns {

bool PeriodicTransform::valid(int dimension) const
{
    if (dimension != 2 && dimension != 3) {
        return false;
    }
    constexpr Real tolerance = 1.0e-12;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            if (!std::isfinite(rotation[static_cast<std::size_t>(row)]
                                      [static_cast<std::size_t>(column)])) {
                return false;
            }
            Real product = 0.0;
            for (int inner = 0; inner < 3; ++inner) {
                product += rotation[static_cast<std::size_t>(row)]
                                   [static_cast<std::size_t>(inner)]
                    * rotation[static_cast<std::size_t>(column)]
                              [static_cast<std::size_t>(inner)];
            }
            const Real expected = row == column ? 1.0 : 0.0;
            if (std::abs(product - expected) > tolerance) {
                return false;
            }
        }
        if (!std::isfinite(translation[static_cast<std::size_t>(row)])) {
            return false;
        }
    }
    const Real determinant
        = rotation[0][0]
            * (rotation[1][1] * rotation[2][2]
                - rotation[1][2] * rotation[2][1])
        - rotation[0][1]
            * (rotation[1][0] * rotation[2][2]
                - rotation[1][2] * rotation[2][0])
        + rotation[0][2]
            * (rotation[1][0] * rotation[2][1]
                - rotation[1][1] * rotation[2][0]);
    if (std::abs(determinant - 1.0) > tolerance) {
        return false;
    }
    if (dimension == 2
        && (std::abs(rotation[0][2]) > tolerance
            || std::abs(rotation[1][2]) > tolerance
            || std::abs(rotation[2][0]) > tolerance
            || std::abs(rotation[2][1]) > tolerance
            || std::abs(rotation[2][2] - 1.0) > tolerance
            || std::abs(translation[2]) > tolerance)) {
        return false;
    }
    return true;
}

std::array<Real, 3> PeriodicTransform::apply_point(
    const std::array<Real, 3>& point) const
{
    auto result = apply_vector(point);
    for (int component = 0; component < 3; ++component) {
        result[static_cast<std::size_t>(component)]
            += translation[static_cast<std::size_t>(component)];
    }
    return result;
}

std::array<Real, 3> PeriodicTransform::apply_vector(
    const std::array<Real, 3>& vector) const
{
    std::array<Real, 3> result {};
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result[static_cast<std::size_t>(row)]
                += rotation[static_cast<std::size_t>(row)]
                           [static_cast<std::size_t>(column)]
                * vector[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

PeriodicTransform PeriodicTransform::inverse() const
{
    PeriodicTransform result;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result.rotation[static_cast<std::size_t>(row)]
                           [static_cast<std::size_t>(column)]
                = rotation[static_cast<std::size_t>(column)]
                          [static_cast<std::size_t>(row)];
        }
    }
    const auto rotated_translation = result.apply_vector(translation);
    for (int component = 0; component < 3; ++component) {
        result.translation[static_cast<std::size_t>(component)]
            = -rotated_translation[static_cast<std::size_t>(component)];
    }
    return result;
}

} // namespace wcns
