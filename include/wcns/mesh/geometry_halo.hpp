#pragma once

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/structured_mesh.hpp>

#include <vector>
#include <unordered_map>

namespace wcns {

enum class GeometryMessageKind : int {
    GeometryVertex = 0,
    GeometryOperand = 1,
    SharedMetric = 2,
};

enum class GeometryOperandStage : int {
    None = 0,
    CenterCoordinates = 1,
    FirstDerivative = 2,
    MetricProduct = 3,
    JacobianProduct = 4,
};

struct GeometryExchangeDescriptor {
    ConnectionId connection = invalid_connection_id;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId donor_rank = invalid_rank_id;
    BlockId shared_face_owner = invalid_block_id;
    GeometryMessageKind kind = GeometryMessageKind::GeometryVertex;
    GeometryOperandStage stage = GeometryOperandStage::None;
    int halo_width = 0;
    IndexTransform index_transform {};
    PeriodicTransform periodic {};
    std::vector<BlockId> donor_path;

    [[nodiscard]] int message_tag(int tag_base = 4096) const;
};

class GeometryHaloPlan {
public:
    [[nodiscard]] static GeometryHaloPlan build(
        const StructuredMesh& mesh,
        const AlgorithmProfile& profile);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] const std::vector<GeometryExchangeDescriptor>& exchanges() const noexcept
    {
        return exchanges_;
    }

private:
    friend struct GeometryHaloPlanBuilderAccess;
    AlgorithmProfileKind profile_ = AlgorithmProfileKind::PhengleiWcns;
    std::vector<GeometryExchangeDescriptor> exchanges_;
};

class SharedMetricSynchronizer {
public:
    static void synchronize(
        const StructuredMesh& mesh,
        std::unordered_map<BlockId, MetricField>& metrics);

private:
    static FaceAreaVectors& face_vectors(MetricField& metric, Axis axis);
};

} // namespace wcns
