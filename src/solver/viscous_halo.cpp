#include <wcns/solver/viscous_halo.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace wcns {
namespace {

constexpr std::uint64_t maximum_exact_message_version = 9007199254740992ULL;
constexpr int operand_tag_base = 16384;

std::array<std::array<Real, 3>, 3> transpose(
    const std::array<std::array<Real, 3>, 3>& matrix)
{
    std::array<std::array<Real, 3>, 3> result {};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]
                = matrix[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
        }
    }
    return result;
}

std::array<std::array<Real, 3>, 3> transform_tensor(
    const std::array<std::array<Real, 3>, 3>& donor,
    const PeriodicTransform& periodic)
{
    const auto& q = periodic.rotation;
    const auto qt = transpose(q);
    std::array<std::array<Real, 3>, 3> temporary {};
    std::array<std::array<Real, 3>, 3> result {};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                temporary[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]
                    += qt[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]
                    * donor[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                result[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]
                    += temporary[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]
                    * q[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
            }
        }
    }
    return result;
}

void require_version(std::uint64_t version, const char* label)
{
    if (version == 0 || version > maximum_exact_message_version) {
        throw std::invalid_argument(std::string(label) + " version is invalid");
    }
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("viscous halo message exceeds MPI int count");
    }
    return static_cast<int>(count);
}
#endif

GradientOperandState load_operand(
    const GradientOperandFaceField& field, Axis axis, Index3 index)
{
    GradientOperandState result {};
    const auto& values = field.field(axis);
    for (int component = 0; component < gradient_operand_components; ++component) {
        result[static_cast<std::size_t>(component)]
            = values(index.i, index.j, index.k, component);
    }
    return result;
}

void store_operand(
    GradientOperandFaceField& field, Axis axis, Index3 index,
    const GradientOperandState& state)
{
    auto& values = field.field(axis);
    for (int component = 0; component < gradient_operand_components; ++component) {
        values(index.i, index.j, index.k, component)
            = state[static_cast<std::size_t>(component)];
    }
}

void validate_operand_field(
    const GradientOperandFaceField& field,
    const FaceFluxExchangeDescriptor& descriptor)
{
    if (field.profile() != descriptor.profile || field.version() != descriptor.version) {
        throw std::invalid_argument("gradient operand field metadata mismatch");
    }
}

PrimitiveGradients load_gradients(
    const PrimitiveGradientField& field, Index3 index)
{
    PrimitiveGradients result {};
    for (int variable = 0; variable < viscous_primitive_components; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            result[static_cast<std::size_t>(variable)][static_cast<std::size_t>(direction)]
                = field(index, static_cast<ViscousPrimitive>(variable), direction);
        }
    }
    return result;
}

void store_gradients(
    PrimitiveGradientField& field, Index3 index,
    const PrimitiveGradients& gradients)
{
    for (int variable = 0; variable < viscous_primitive_components; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            field(index, static_cast<ViscousPrimitive>(variable), direction)
                = gradients[static_cast<std::size_t>(variable)]
                    [static_cast<std::size_t>(direction)];
        }
    }
}

void validate_gradient_field(
    const PrimitiveGradientField& field,
    const GradientExchangeDescriptor& descriptor)
{
    if (field.profile() != descriptor.profile || field.version() != descriptor.version
        || field.dimension() != descriptor.dimension) {
        throw std::invalid_argument("primitive gradient field metadata mismatch");
    }
}

const ConnectivityPatch& find_connection(
    const StructuredMesh& mesh,
    const DirectedExchange& exchange)
{
    const auto& block = mesh.block(exchange.halo.receiver_block);
    const auto iterator = std::find_if(
        block.connectivities.begin(), block.connectivities.end(),
        [&](const ConnectivityPatch& connection) {
            return connection.name == exchange.halo.connectivity_name
                && connection.donor_block == exchange.halo.donor_block;
        });
    if (iterator == block.connectivities.end()) {
        throw TopologyError("gradient halo cannot find its connectivity descriptor");
    }
    return *iterator;
}

int ghost_layer(Index3 index, Extent3 extent)
{
    int result = 0;
    for (int axis = 0; axis < 3; ++axis) {
        const int value = index[static_cast<std::size_t>(axis)];
        const int count = extent[static_cast<std::size_t>(axis)];
        if (value < 0) result = std::max(result, -value);
        if (value >= count) result = std::max(result, value - count + 1);
    }
    return result;
}

} // namespace

GradientOperandState transform_gradient_operand_for_receiver(
    const GradientOperandState& donor,
    const FaceFluxExchangeDescriptor& descriptor,
    int dimension)
{
    std::array<std::array<Real, 3>, 3> velocity_tensor {};
    for (int variable = 0; variable < 3; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            velocity_tensor[static_cast<std::size_t>(variable)]
                [static_cast<std::size_t>(direction)]
                = donor[static_cast<std::size_t>(variable * 3 + direction)];
        }
    }
    const auto transformed = transform_tensor(velocity_tensor, descriptor.periodic);
    const auto transformed_temperature = descriptor.periodic.inverse().apply_vector({{
        donor[9], donor[10], donor[11],
    }});
    GradientOperandState result {};
    for (int variable = 0; variable < 3; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            result[static_cast<std::size_t>(variable * 3 + direction)]
                = descriptor.orientation
                * transformed[static_cast<std::size_t>(variable)]
                    [static_cast<std::size_t>(direction)];
        }
    }
    for (int direction = 0; direction < 3; ++direction) {
        result[static_cast<std::size_t>(9 + direction)]
            = descriptor.orientation
            * transformed_temperature[static_cast<std::size_t>(direction)];
    }
    if (dimension == 2) {
        for (int direction = 0; direction < 3; ++direction) {
            result[static_cast<std::size_t>(6 + direction)] = 0.0;
        }
        for (int variable = 0; variable < viscous_primitive_components; ++variable) {
            result[static_cast<std::size_t>(variable * 3 + 2)] = 0.0;
        }
    }
    return result;
}

GradientOperandFaceHaloPlan GradientOperandFaceHaloPlan::build(
    const StructuredMesh& mesh,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    GradientOperandFaceHaloPlan result;
    const auto base = FaceFluxHaloPlan::build(mesh, profile, version);
    result.exchanges_ = base.exchanges();
    return result;
}

void GradientOperandFieldRegistry::add(
    BlockId block, GradientOperandFaceField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("gradient operand registry has an invalid block");
    }
}

bool GradientOperandFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

GradientOperandFaceField& GradientOperandFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("gradient operand field is not registered");
    }
    return *iterator->second;
}

void GradientOperandFaceHaloExchanger::exchange(
    const GradientOperandFieldRegistry& fields) const
{
    struct Pending {
        const FaceFluxExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };
    const RankId rank = mpi_.rank();
    std::vector<Pending> receives;
    std::vector<Pending> sends;
    for (const auto& descriptor : plan_.exchanges()) {
        const std::size_t count = 1 + descriptor.pairs.size()
            * static_cast<std::size_t>(gradient_operand_components);
        if (descriptor.receiver_rank == rank && descriptor.donor_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            const auto& donor = fields.field(descriptor.donor_block);
            validate_operand_field(receiver, descriptor);
            validate_operand_field(donor, descriptor);
            for (const auto& pair : descriptor.pairs) {
                store_operand(receiver, descriptor.receiver_axis, pair.receiver,
                    transform_gradient_operand_for_receiver(
                        load_operand(donor, descriptor.donor_axis, pair.donor),
                        descriptor, receiver.dimension()));
            }
        } else if (descriptor.receiver_rank == rank) {
            validate_operand_field(fields.field(descriptor.receiver_block), descriptor);
            receives.push_back({&descriptor, std::vector<Real>(count)});
        } else if (descriptor.donor_rank == rank) {
            const auto& donor = fields.field(descriptor.donor_block);
            validate_operand_field(donor, descriptor);
            Pending pending {&descriptor, std::vector<Real>(count)};
            pending.values[0] = static_cast<Real>(descriptor.version);
            std::size_t offset = 1;
            for (const auto& pair : descriptor.pairs) {
                const auto state = load_operand(
                    donor, descriptor.donor_axis, pair.donor);
                for (const Real value : state) pending.values[offset++] = value;
            }
            sends.push_back(std::move(pending));
        }
    }
#if WCNS_HAS_MPI
    std::vector<MPI_Request> requests(receives.size() + sends.size(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives) {
        check_mpi(MPI_Irecv(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->donor_rank,
            pending.descriptor->message_tag(operand_tag_base), mpi_.communicator(),
            &requests[request++]), "MPI_Irecv gradient operand");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->receiver_rank,
            pending.descriptor->message_tag(operand_tag_base), mpi_.communicator(),
            &requests[request++]), "MPI_Isend gradient operand");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall gradient operand");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote gradient operand exchange requires MPI");
    }
#endif
    for (const auto& pending : receives) {
        if (pending.values[0] != static_cast<Real>(pending.descriptor->version)) {
            throw MpiError("gradient operand message version mismatch");
        }
        auto& receiver = fields.field(pending.descriptor->receiver_block);
        std::size_t offset = 1;
        for (const auto& pair : pending.descriptor->pairs) {
            GradientOperandState donor {};
            for (auto& value : donor) value = pending.values[offset++];
            store_operand(receiver, pending.descriptor->receiver_axis, pair.receiver,
                transform_gradient_operand_for_receiver(
                    donor, *pending.descriptor, receiver.dimension()));
        }
    }
}

int GradientExchangeDescriptor::message_tag(int tag_base) const
{
    if (tag_base < 0 || connection < 0 || receiver_block < 0 || donor_block < 0) {
        throw TopologyError("gradient message tag inputs are invalid");
    }
    const int direction = receiver_block < donor_block ? 0 : 1;
    const long long tag = static_cast<long long>(tag_base)
        + 8LL * connection + 2LL * static_cast<int>(profile) + direction;
    if (tag > std::numeric_limits<int>::max()) {
        throw TopologyError("gradient message tag exceeds int range");
    }
    return static_cast<int>(tag);
}

PrimitiveGradients transform_primitive_gradients_for_receiver(
    const PrimitiveGradients& donor,
    const GradientExchangeDescriptor& descriptor)
{
    std::array<std::array<Real, 3>, 3> velocity {};
    for (int variable = 0; variable < 3; ++variable) {
        velocity[static_cast<std::size_t>(variable)]
            = donor[static_cast<std::size_t>(variable)];
    }
    const auto transformed_velocity = transform_tensor(velocity, descriptor.periodic);
    const auto transformed_temperature = descriptor.periodic.inverse().apply_vector(
        donor[static_cast<int>(ViscousPrimitive::Temperature)]);
    PrimitiveGradients result {};
    for (int variable = 0; variable < 3; ++variable) {
        result[static_cast<std::size_t>(variable)]
            = transformed_velocity[static_cast<std::size_t>(variable)];
    }
    result[static_cast<int>(ViscousPrimitive::Temperature)]
        = transformed_temperature;
    if (descriptor.dimension == 2) {
        result[static_cast<int>(ViscousPrimitive::VelocityZ)] = {{0.0, 0.0, 0.0}};
        for (auto& gradient : result) gradient[2] = 0.0;
    }
    return result;
}

GradientHaloPlan GradientHaloPlan::build(
    const StructuredMesh& mesh,
    const DistributedTopology& topology,
    const AlgorithmProfile& profile,
    std::uint64_t version)
{
    require_version(version, "gradient halo");
    GradientHaloPlan result;
    const int width = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 2 : 3;
    for (const auto& exchange : topology.exchanges()) {
        const auto& connection = find_connection(mesh, exchange);
        GradientExchangeDescriptor descriptor;
        descriptor.connection = exchange.connection;
        descriptor.receiver_block = exchange.halo.receiver_block;
        descriptor.donor_block = exchange.halo.donor_block;
        descriptor.receiver_rank = exchange.receiver_rank;
        descriptor.donor_rank = exchange.donor_rank;
        descriptor.periodic = connection.periodic;
        descriptor.profile = profile.kind();
        descriptor.version = version;
        descriptor.dimension = mesh.block(descriptor.receiver_block).cell_dimension();
        const auto receiver_extent = mesh.block(descriptor.receiver_block).cell_extent();
        for (const auto& pair : exchange.halo.cell_pairs) {
            if (ghost_layer(pair.receiver_ghost, receiver_extent) <= width) {
                descriptor.pairs.push_back(pair);
            }
        }
        result.exchanges_.push_back(std::move(descriptor));
    }
    return result;
}

void GradientFieldRegistry::add(BlockId block, PrimitiveGradientField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("gradient registry has an invalid block");
    }
}

bool GradientFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

PrimitiveGradientField& GradientFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("primitive gradient field is not registered");
    }
    return *iterator->second;
}

void GradientHaloExchanger::exchange(const GradientFieldRegistry& fields) const
{
    struct Pending {
        const GradientExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };
    const RankId rank = mpi_.rank();
    std::vector<Pending> receives;
    std::vector<Pending> sends;
    for (const auto& descriptor : plan_.exchanges()) {
        const std::size_t count = 1 + descriptor.pairs.size()
            * static_cast<std::size_t>(gradient_operand_components);
        if (descriptor.receiver_rank == rank && descriptor.donor_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            const auto& donor = fields.field(descriptor.donor_block);
            validate_gradient_field(receiver, descriptor);
            validate_gradient_field(donor, descriptor);
            for (const auto& pair : descriptor.pairs) {
                store_gradients(receiver, pair.receiver_ghost,
                    transform_primitive_gradients_for_receiver(
                        load_gradients(donor, pair.donor_interior), descriptor));
            }
        } else if (descriptor.receiver_rank == rank) {
            validate_gradient_field(fields.field(descriptor.receiver_block), descriptor);
            receives.push_back({&descriptor, std::vector<Real>(count)});
        } else if (descriptor.donor_rank == rank) {
            const auto& donor = fields.field(descriptor.donor_block);
            validate_gradient_field(donor, descriptor);
            Pending pending {&descriptor, std::vector<Real>(count)};
            pending.values[0] = static_cast<Real>(descriptor.version);
            std::size_t offset = 1;
            for (const auto& pair : descriptor.pairs) {
                const auto gradients = load_gradients(donor, pair.donor_interior);
                for (const auto& gradient : gradients) {
                    for (const Real value : gradient) pending.values[offset++] = value;
                }
            }
            sends.push_back(std::move(pending));
        }
    }
#if WCNS_HAS_MPI
    std::vector<MPI_Request> requests(receives.size() + sends.size(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives) {
        check_mpi(MPI_Irecv(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->donor_rank, pending.descriptor->message_tag(),
            mpi_.communicator(), &requests[request++]), "MPI_Irecv gradient");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(
            pending.values.data(), mpi_count(pending.values.size()), MPI_DOUBLE,
            pending.descriptor->receiver_rank, pending.descriptor->message_tag(),
            mpi_.communicator(), &requests[request++]), "MPI_Isend gradient");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall gradient");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote gradient exchange requires MPI");
    }
#endif
    for (const auto& pending : receives) {
        if (pending.values[0] != static_cast<Real>(pending.descriptor->version)) {
            throw MpiError("gradient message version mismatch");
        }
        auto& receiver = fields.field(pending.descriptor->receiver_block);
        std::size_t offset = 1;
        for (const auto& pair : pending.descriptor->pairs) {
            PrimitiveGradients donor {};
            for (auto& gradient : donor) {
                for (auto& value : gradient) value = pending.values[offset++];
            }
            store_gradients(receiver, pair.receiver_ghost,
                transform_primitive_gradients_for_receiver(
                    donor, *pending.descriptor));
        }
    }
}

} // namespace wcns
