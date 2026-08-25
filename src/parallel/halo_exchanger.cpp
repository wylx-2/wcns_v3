#include <wcns/parallel/halo_exchanger.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace wcns {
namespace {

struct PendingBuffer {
    const DirectedExchange* exchange = nullptr;
    std::vector<Real> values;
};

std::size_t value_count(const DirectedExchange& exchange, int components)
{
    const auto cells = exchange.halo.cell_pairs.size();
    const auto component_count = static_cast<std::size_t>(components);
    if (cells > std::numeric_limits<std::size_t>::max() / component_count) {
        throw std::overflow_error("halo message value count exceeds size_t range");
    }
    return cells * component_count;
}

void validate_field_for_exchange(
    const Field<Real>& field,
    const DirectedExchange& exchange,
    bool receiver)
{
    const auto extent = field.interior_extent();
    const int ghost = field.ghost_width();
    for (const auto& pair : exchange.halo.cell_pairs) {
        const auto index = receiver ? pair.receiver_ghost : pair.donor_interior;
        const bool valid = receiver
            ? index.i >= -ghost && index.i < extent.ni + ghost && index.j >= -ghost
                && index.j < extent.nj + ghost && index.k >= -ghost
                && index.k < extent.nk + ghost
            : index.i >= 0 && index.i < extent.ni && index.j >= 0
                && index.j < extent.nj && index.k >= 0 && index.k < extent.nk;
        if (!valid) {
            throw std::invalid_argument(
                receiver ? "halo receiver index is outside field storage"
                         : "halo donor index is outside field interior");
        }
    }
}

void copy_local(
    const DirectedExchange& exchange,
    const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(exchange.halo.receiver_block);
    const auto& donor = fields.field(exchange.halo.donor_block);
    validate_field_for_exchange(receiver, exchange, true);
    validate_field_for_exchange(donor, exchange, false);
    for (const auto& pair : exchange.halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            receiver(
                pair.receiver_ghost.i,
                pair.receiver_ghost.j,
                pair.receiver_ghost.k,
                component)
                = donor(
                    pair.donor_interior.i,
                    pair.donor_interior.j,
                    pair.donor_interior.k,
                    component);
        }
    }
}

void pack_send(
    PendingBuffer& pending,
    const BlockFieldRegistry& fields)
{
    const auto& donor = fields.field(pending.exchange->halo.donor_block);
    validate_field_for_exchange(donor, *pending.exchange, false);
    std::size_t output = 0;
    for (const auto& pair : pending.exchange->halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            pending.values[output++] = donor(
                pair.donor_interior.i,
                pair.donor_interior.j,
                pair.donor_interior.k,
                component);
        }
    }
}

void unpack_receive(
    const PendingBuffer& pending,
    const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(pending.exchange->halo.receiver_block);
    validate_field_for_exchange(receiver, *pending.exchange, true);
    std::size_t input = 0;
    for (const auto& pair : pending.exchange->halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            receiver(
                pair.receiver_ghost.i,
                pair.receiver_ghost.j,
                pair.receiver_ghost.k,
                component)
                = pending.values[input++];
        }
    }
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t values)
{
    if (values > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("halo message exceeds MPI int count range");
    }
    return static_cast<int>(values);
}
#endif

} // namespace

BlockFieldRegistry::BlockFieldRegistry(int components)
    : components_(components)
{
    if (components <= 0) {
        throw std::invalid_argument("field registry component count must be positive");
    }
}

void BlockFieldRegistry::add(BlockId block, Field<Real>& field)
{
    if (block < 0 || field.components() != components_) {
        throw std::invalid_argument("registered block field has incompatible metadata");
    }
    if (!fields_.emplace(block, &field).second) {
        throw std::invalid_argument("a block field is already registered");
    }
}

bool BlockFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

Field<Real>& BlockFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("block field is not registered on this rank");
    }
    return *iterator->second;
}

HaloExchanger::HaloExchanger(
    const MpiRuntime& mpi,
    const DistributedTopology& topology,
    int distribution_rank_count)
    : mpi_(mpi)
    , topology_(topology)
{
    if (distribution_rank_count != mpi.size()) {
        throw std::invalid_argument("MPI size differs from the block distribution rank count");
    }
}

void HaloExchanger::exchange(const BlockFieldRegistry& fields) const
{
    const auto rank = mpi_.rank();
    const auto local = topology_.local_copies(rank);
    const auto receives = topology_.receives(rank);
    const auto sends = topology_.sends(rank);

    for (const auto* exchange : local) {
        copy_local(*exchange, fields);
    }

    std::vector<PendingBuffer> receive_buffers;
    receive_buffers.reserve(receives.size());
    for (const auto* exchange : receives) {
        if (!fields.contains(exchange->halo.receiver_block)) {
            throw std::invalid_argument("receiver field is missing on its owner rank");
        }
        receive_buffers.push_back(
            {exchange, std::vector<Real>(value_count(*exchange, fields.components()))});
    }
    std::vector<PendingBuffer> send_buffers;
    send_buffers.reserve(sends.size());
    for (const auto* exchange : sends) {
        if (!fields.contains(exchange->halo.donor_block)) {
            throw std::invalid_argument("donor field is missing on its owner rank");
        }
        send_buffers.push_back(
            {exchange, std::vector<Real>(value_count(*exchange, fields.components()))});
        pack_send(send_buffers.back(), fields);
    }

#if WCNS_HAS_MPI
    int* tag_upper_bound = nullptr;
    int has_tag_upper_bound = 0;
    check_mpi(
        MPI_Comm_get_attr(
            mpi_.communicator(), MPI_TAG_UB, &tag_upper_bound, &has_tag_upper_bound),
        "MPI_Comm_get_attr MPI_TAG_UB");
    if (has_tag_upper_bound == 0 || tag_upper_bound == nullptr) {
        throw MpiError("MPI_TAG_UB is unavailable");
    }

    std::vector<MPI_Request> requests(
        receive_buffers.size() + send_buffers.size(), MPI_REQUEST_NULL);
    std::size_t request_index = 0;
    for (auto& pending : receive_buffers) {
        const int tag = pending.exchange->message_tag();
        if (tag > *tag_upper_bound) {
            throw MpiError("halo receive tag exceeds MPI_TAG_UB");
        }
        check_mpi(
            MPI_Irecv(
                pending.values.data(),
                mpi_count(pending.values.size()),
                MPI_DOUBLE,
                pending.exchange->donor_rank,
                tag,
                mpi_.communicator(),
                &requests[request_index++]),
            "MPI_Irecv halo");
    }
    for (auto& pending : send_buffers) {
        const int tag = pending.exchange->message_tag();
        if (tag > *tag_upper_bound) {
            throw MpiError("halo send tag exceeds MPI_TAG_UB");
        }
        check_mpi(
            MPI_Isend(
                pending.values.data(),
                mpi_count(pending.values.size()),
                MPI_DOUBLE,
                pending.exchange->receiver_rank,
                tag,
                mpi_.communicator(),
                &requests[request_index++]),
            "MPI_Isend halo");
    }
    if (!requests.empty()) {
        check_mpi(
            MPI_Waitall(
                static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall halo");
    }
#else
    if (!receive_buffers.empty() || !send_buffers.empty()) {
        throw MpiError("remote halo exchange requires WCNS_ENABLE_MPI");
    }
#endif

    for (const auto& pending : receive_buffers) {
        unpack_receive(pending, fields);
    }
}

} // namespace wcns
