#pragma once

#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/viscous_gradient.hpp>

#include <unordered_map>
#include <vector>

namespace wcns {

using GradientOperandState = std::array<Real, gradient_operand_components>;

[[nodiscard]] GradientOperandState transform_gradient_operand_for_receiver(
    const GradientOperandState& donor,
    const FaceFluxExchangeDescriptor& descriptor,
    int dimension);

class GradientOperandFaceHaloPlan {
public:
    [[nodiscard]] static GradientOperandFaceHaloPlan build(
        const StructuredMesh& mesh,
        const AlgorithmProfile& profile,
        std::uint64_t version);

    [[nodiscard]] const std::vector<FaceFluxExchangeDescriptor>& exchanges() const noexcept
    {
        return exchanges_;
    }

private:
    std::vector<FaceFluxExchangeDescriptor> exchanges_;
};

class GradientOperandFieldRegistry {
public:
    void add(BlockId block, GradientOperandFaceField& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] GradientOperandFaceField& field(BlockId block) const;

private:
    std::unordered_map<BlockId, GradientOperandFaceField*> fields_;
};

class GradientOperandFaceHaloExchanger {
public:
    GradientOperandFaceHaloExchanger(
        const MpiRuntime& mpi,
        const GradientOperandFaceHaloPlan& plan)
        : mpi_(mpi), plan_(plan)
    {
    }

    void exchange(const GradientOperandFieldRegistry& fields) const;

private:
    const MpiRuntime& mpi_;
    const GradientOperandFaceHaloPlan& plan_;
};

struct GradientExchangeDescriptor {
    ConnectionId connection = invalid_connection_id;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId receiver_rank = invalid_rank_id;
    RankId donor_rank = invalid_rank_id;
    PeriodicTransform periodic {};
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    std::uint64_t version = 0;
    int dimension = 0;
    std::vector<HaloCellPair> pairs;

    [[nodiscard]] int message_tag(int tag_base = 20480) const;
};

[[nodiscard]] PrimitiveGradients transform_primitive_gradients_for_receiver(
    const PrimitiveGradients& donor,
    const GradientExchangeDescriptor& descriptor);

class GradientHaloPlan {
public:
    [[nodiscard]] static GradientHaloPlan build(
        const StructuredMesh& mesh,
        const DistributedTopology& topology,
        const AlgorithmProfile& profile,
        std::uint64_t version);

    [[nodiscard]] const std::vector<GradientExchangeDescriptor>& exchanges() const noexcept
    {
        return exchanges_;
    }

private:
    std::vector<GradientExchangeDescriptor> exchanges_;
};

class GradientFieldRegistry {
public:
    void add(BlockId block, PrimitiveGradientField& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] PrimitiveGradientField& field(BlockId block) const;

private:
    std::unordered_map<BlockId, PrimitiveGradientField*> fields_;
};

class GradientHaloExchanger {
public:
    GradientHaloExchanger(const MpiRuntime& mpi, const GradientHaloPlan& plan)
        : mpi_(mpi), plan_(plan)
    {
    }

    void exchange(const GradientFieldRegistry& fields) const;

private:
    const MpiRuntime& mpi_;
    const GradientHaloPlan& plan_;
};

} // namespace wcns
