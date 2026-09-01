#pragma once

#include <wcns/core/array3d.hpp>
#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/structured_mesh.hpp>

#include <unordered_map>
#include <vector>

namespace wcns {

class ConservationWeightError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct LineConservationWeights {
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    int cell_count = 0;
    bool periodic = false;
    std::vector<Real> cell_weights;
    Real maximum_residual = 0.0;
};

[[nodiscard]] LineConservationWeights build_line_conservation_weights(
    const AlgorithmProfile& profile,
    int cell_count,
    bool periodic = false);

struct BlockConservationWeights {
    explicit BlockConservationWeights(Extent3 cell_extent)
        : cell(cell_extent)
        , physical_i_face({cell_extent.ni + 1, cell_extent.nj, cell_extent.nk})
        , physical_j_face({cell_extent.ni, cell_extent.nj + 1, cell_extent.nk})
        , physical_k_face({cell_extent.ni, cell_extent.nj, cell_extent.nk + 1})
    {
        physical_i_face.fill(0.0);
        physical_j_face.fill(0.0);
        physical_k_face.fill(0.0);
    }

    Array3D<Real> cell;
    Array3D<Real> physical_i_face;
    Array3D<Real> physical_j_face;
    Array3D<Real> physical_k_face;
};

class GlobalConservationWeights {
public:
    [[nodiscard]] static GlobalConservationWeights build(
        const StructuredMesh& mesh,
        const AlgorithmProfile& profile);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] const BlockConservationWeights& block(BlockId id) const;
    [[nodiscard]] Real maximum_line_residual() const noexcept
    {
        return maximum_line_residual_;
    }
    [[nodiscard]] Real maximum_shared_face_mismatch() const noexcept
    {
        return maximum_shared_face_mismatch_;
    }

private:
    AlgorithmProfileKind profile_ = AlgorithmProfileKind::PhengleiWcns;
    std::unordered_map<BlockId, BlockConservationWeights> blocks_;
    Real maximum_line_residual_ = 0.0;
    Real maximum_shared_face_mismatch_ = 0.0;
};

} // namespace wcns
