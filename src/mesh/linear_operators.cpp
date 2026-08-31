#include <wcns/mesh/linear_operators.hpp>

#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

StencilRow row(std::initializer_list<StencilEntry> entries)
{
    return StencilRow(entries);
}

std::vector<Real> apply(
    const std::vector<StencilRow>& rows,
    const std::vector<Real>& values,
    std::size_t expected_size,
    const char* label)
{
    if (values.size() != expected_size) {
        throw std::invalid_argument(std::string(label) + " input has the wrong length");
    }
    std::vector<Real> result(rows.size(), 0.0);
    for (std::size_t output = 0; output < rows.size(); ++output) {
        for (const auto [input, coefficient] : rows[output]) {
            result[output] += coefficient * values[static_cast<std::size_t>(input)];
        }
    }
    return result;
}

std::vector<StencilRow> ph_interpolation(int count)
{
    std::vector<StencilRow> rows(static_cast<std::size_t>(count + 1));
    rows[0] = row({{0, 35.0 / 16.0}, {1, -35.0 / 16.0},
        {2, 21.0 / 16.0}, {3, -5.0 / 16.0}});
    rows[1] = row({{0, 5.0 / 16.0}, {1, 15.0 / 16.0},
        {2, -5.0 / 16.0}, {3, 1.0 / 16.0}});
    for (int face = 2; face <= count - 2; ++face) {
        rows[static_cast<std::size_t>(face)] = row({
            {face - 2, -1.0 / 16.0},
            {face - 1, 9.0 / 16.0},
            {face, 9.0 / 16.0},
            {face + 1, -1.0 / 16.0},
        });
    }
    rows[static_cast<std::size_t>(count - 1)] = row({
        {count - 4, 1.0 / 16.0},
        {count - 3, -5.0 / 16.0},
        {count - 2, 15.0 / 16.0},
        {count - 1, 5.0 / 16.0},
    });
    rows[static_cast<std::size_t>(count)] = row({
        {count - 4, -5.0 / 16.0},
        {count - 3, 21.0 / 16.0},
        {count - 2, -35.0 / 16.0},
        {count - 1, 35.0 / 16.0},
    });
    return rows;
}

std::vector<StencilRow> ph_derivative(int count)
{
    std::vector<StencilRow> rows(static_cast<std::size_t>(count));
    rows[0] = row({{0, -1.0}, {1, 1.0}});
    for (int cell = 1; cell < count - 1; ++cell) {
        rows[static_cast<std::size_t>(cell)] = row({
            {cell - 1, 1.0 / 24.0},
            {cell, -9.0 / 8.0},
            {cell + 1, 9.0 / 8.0},
            {cell + 2, -1.0 / 24.0},
        });
    }
    rows[static_cast<std::size_t>(count - 1)]
        = row({{count - 1, -1.0}, {count, 1.0}});
    return rows;
}

std::vector<StencilRow> scmm_interpolation(int count)
{
    std::vector<StencilRow> rows(static_cast<std::size_t>(count + 1));
    rows[0] = row({{0, 315.0 / 128.0}, {1, -420.0 / 128.0},
        {2, 378.0 / 128.0}, {3, -180.0 / 128.0}, {4, 35.0 / 128.0}});
    rows[1] = row({{0, 35.0 / 128.0}, {1, 140.0 / 128.0},
        {2, -70.0 / 128.0}, {3, 28.0 / 128.0}, {4, -5.0 / 128.0}});
    rows[2] = row({{0, -5.0 / 128.0}, {1, 60.0 / 128.0},
        {2, 90.0 / 128.0}, {3, -20.0 / 128.0}, {4, 3.0 / 128.0}});
    for (int face = 3; face <= count - 3; ++face) {
        rows[static_cast<std::size_t>(face)] = row({
            {face - 3, 3.0 / 256.0},
            {face - 2, -25.0 / 256.0},
            {face - 1, 150.0 / 256.0},
            {face, 150.0 / 256.0},
            {face + 1, -25.0 / 256.0},
            {face + 2, 3.0 / 256.0},
        });
    }
    rows[static_cast<std::size_t>(count - 2)] = row({
        {count - 5, 3.0 / 128.0}, {count - 4, -20.0 / 128.0},
        {count - 3, 90.0 / 128.0}, {count - 2, 60.0 / 128.0},
        {count - 1, -5.0 / 128.0}});
    rows[static_cast<std::size_t>(count - 1)] = row({
        {count - 5, -5.0 / 128.0}, {count - 4, 28.0 / 128.0},
        {count - 3, -70.0 / 128.0}, {count - 2, 140.0 / 128.0},
        {count - 1, 35.0 / 128.0}});
    rows[static_cast<std::size_t>(count)] = row({
        {count - 5, 35.0 / 128.0}, {count - 4, -180.0 / 128.0},
        {count - 3, 378.0 / 128.0}, {count - 2, -420.0 / 128.0},
        {count - 1, 315.0 / 128.0}});
    return rows;
}

std::vector<StencilRow> scmm_derivative(int count)
{
    std::vector<StencilRow> rows(static_cast<std::size_t>(count));
    rows[0] = row({{0, -22.0 / 24.0}, {1, 17.0 / 24.0},
        {2, 9.0 / 24.0}, {3, -5.0 / 24.0}, {4, 1.0 / 24.0}});
    rows[1] = row({{0, 1.0 / 24.0}, {1, -27.0 / 24.0},
        {2, 27.0 / 24.0}, {3, -1.0 / 24.0}});
    for (int cell = 2; cell <= count - 3; ++cell) {
        rows[static_cast<std::size_t>(cell)] = row({
            {cell - 2, -9.0 / 1920.0},
            {cell - 1, 125.0 / 1920.0},
            {cell, -2250.0 / 1920.0},
            {cell + 1, 2250.0 / 1920.0},
            {cell + 2, -125.0 / 1920.0},
            {cell + 3, 9.0 / 1920.0},
        });
    }
    rows[static_cast<std::size_t>(count - 2)] = row({
        {count - 3, 1.0 / 24.0}, {count - 2, -27.0 / 24.0},
        {count - 1, 27.0 / 24.0}, {count, -1.0 / 24.0}});
    rows[static_cast<std::size_t>(count - 1)] = row({
        {count - 4, -1.0 / 24.0}, {count - 3, 5.0 / 24.0},
        {count - 2, -9.0 / 24.0}, {count - 1, -17.0 / 24.0},
        {count, 22.0 / 24.0}});
    return rows;
}

StencilRow vertex_center_row(int cell, int cell_count)
{
    if (cell == 0) {
        return row({{0, 63.0 / 256.0}, {1, 315.0 / 256.0},
            {2, -210.0 / 256.0}, {3, 126.0 / 256.0},
            {4, -45.0 / 256.0}, {5, 7.0 / 256.0}});
    }
    if (cell == 1) {
        return row({{0, -7.0 / 256.0}, {1, 105.0 / 256.0},
            {2, 210.0 / 256.0}, {3, -70.0 / 256.0},
            {4, 21.0 / 256.0}, {5, -3.0 / 256.0}});
    }
    if (cell == cell_count - 2) {
        return row({{cell_count - 5, -3.0 / 256.0},
            {cell_count - 4, 21.0 / 256.0},
            {cell_count - 3, -70.0 / 256.0},
            {cell_count - 2, 210.0 / 256.0},
            {cell_count - 1, 105.0 / 256.0},
            {cell_count, -7.0 / 256.0}});
    }
    if (cell == cell_count - 1) {
        return row({{cell_count - 5, 7.0 / 256.0},
            {cell_count - 4, -45.0 / 256.0},
            {cell_count - 3, 126.0 / 256.0},
            {cell_count - 2, -210.0 / 256.0},
            {cell_count - 1, 315.0 / 256.0},
            {cell_count, 63.0 / 256.0}});
    }
    return row({{cell - 2, 3.0 / 256.0}, {cell - 1, -25.0 / 256.0},
        {cell, 150.0 / 256.0}, {cell + 1, 150.0 / 256.0},
        {cell + 2, -25.0 / 256.0}, {cell + 3, 3.0 / 256.0}});
}

} // namespace

LineOperators LineOperators::build(
    const AlgorithmProfile& profile,
    int cell_count)
{
    if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
        if (cell_count < 4) {
            throw ProfileError("phenglei_wcns line operators require at least four cells");
        }
        return LineOperators(
            profile.kind(),
            cell_count,
            ph_interpolation(cell_count),
            ph_derivative(cell_count));
    }
    if (cell_count < 5) {
        throw ProfileError("scmm6_wcns line operators require at least five cells");
    }
    return LineOperators(
        profile.kind(),
        cell_count,
        scmm_interpolation(cell_count),
        scmm_derivative(cell_count));
}

std::vector<Real> LineOperators::interpolate(
    const std::vector<Real>& centers) const
{
    return apply(
        interpolation_, centers, static_cast<std::size_t>(cell_count_), "interpolation");
}

std::vector<Real> LineOperators::differentiate(
    const std::vector<Real>& faces) const
{
    return apply(
        derivative_, faces, static_cast<std::size_t>(cell_count_ + 1), "derivative");
}

std::vector<Real> LineOperators::delta(const std::vector<Real>& centers) const
{
    return differentiate(interpolate(centers));
}

void LineOperators::require_profile(const AlgorithmProfile& profile) const
{
    if (profile.kind() != profile_) {
        throw ProfileError("line operator belongs to a different algorithm profile");
    }
}

std::vector<Real> interpolate_vertices_to_centers_i6(
    const std::vector<Real>& vertices)
{
    if (vertices.size() < 6) {
        throw ProfileError("I6 vertex-to-center interpolation requires at least six vertices");
    }
    const int cell_count = static_cast<int>(vertices.size()) - 1;
    std::vector<StencilRow> rows;
    rows.reserve(static_cast<std::size_t>(cell_count));
    for (int cell = 0; cell < cell_count; ++cell) {
        rows.push_back(vertex_center_row(cell, cell_count));
    }
    return apply(rows, vertices, vertices.size(), "vertex-to-center interpolation");
}

std::vector<Real> grid_delta(const std::vector<Real>& refined_values)
{
    const int count = static_cast<int>(refined_values.size());
    if (count < 5 || count % 2 == 0) {
        throw ProfileError("GridDelta requires an odd refined line with at least five points");
    }
    std::vector<Real> result(refined_values.size(), 0.0);
    result[0] = -3.0 * refined_values[0] + 4.0 * refined_values[1]
        - refined_values[2];
    result[static_cast<std::size_t>(count - 1)]
        = 3.0 * refined_values[static_cast<std::size_t>(count - 1)]
        - 4.0 * refined_values[static_cast<std::size_t>(count - 2)]
        + refined_values[static_cast<std::size_t>(count - 3)];
    for (int index = 1; index < count - 1; ++index) {
        if (index <= 2 || index >= count - 3) {
            result[static_cast<std::size_t>(index)]
                = refined_values[static_cast<std::size_t>(index + 1)]
                - refined_values[static_cast<std::size_t>(index - 1)];
        } else {
            result[static_cast<std::size_t>(index)]
                = 9.0 / 8.0
                    * (refined_values[static_cast<std::size_t>(index + 1)]
                        - refined_values[static_cast<std::size_t>(index - 1)])
                - 1.0 / 24.0
                    * (refined_values[static_cast<std::size_t>(index + 3)]
                        - refined_values[static_cast<std::size_t>(index - 3)]);
        }
    }
    return result;
}

} // namespace wcns
