#pragma once

#include <wcns/core/types.hpp>
#include <wcns/mesh/algorithm_profile.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace wcns {

using StencilEntry = std::pair<int, Real>;
using StencilRow = std::vector<StencilEntry>;

class LineOperators {
public:
    [[nodiscard]] static LineOperators build(
        const AlgorithmProfile& profile,
        int cell_count);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] int cell_count() const noexcept { return cell_count_; }
    [[nodiscard]] const std::vector<StencilRow>& interpolation_rows() const noexcept
    {
        return interpolation_;
    }
    [[nodiscard]] const std::vector<StencilRow>& derivative_rows() const noexcept
    {
        return derivative_;
    }

    [[nodiscard]] std::vector<Real> interpolate(
        const std::vector<Real>& centers) const;
    [[nodiscard]] std::vector<Real> differentiate(
        const std::vector<Real>& faces) const;
    [[nodiscard]] std::vector<Real> delta(
        const std::vector<Real>& centers) const;

    void require_profile(const AlgorithmProfile& profile) const;

private:
    LineOperators(
        AlgorithmProfileKind profile,
        int cell_count,
        std::vector<StencilRow> interpolation,
        std::vector<StencilRow> derivative)
        : profile_(profile)
        , cell_count_(cell_count)
        , interpolation_(std::move(interpolation))
        , derivative_(std::move(derivative))
    {
    }

    AlgorithmProfileKind profile_;
    int cell_count_;
    std::vector<StencilRow> interpolation_;
    std::vector<StencilRow> derivative_;
};

[[nodiscard]] std::vector<Real> interpolate_vertices_to_centers_i6(
    const std::vector<Real>& vertices);

[[nodiscard]] std::vector<Real> grid_delta(
    const std::vector<Real>& refined_values);

} // namespace wcns
