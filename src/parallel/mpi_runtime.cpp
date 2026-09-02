#include <wcns/parallel/mpi_runtime.hpp>

#include <array>
#include <limits>

namespace wcns {

#if WCNS_HAS_MPI

void check_mpi(int status, const char* operation)
{
    if (status == MPI_SUCCESS) {
        return;
    }
    std::array<char, MPI_MAX_ERROR_STRING> message {};
    int length = 0;
    MPI_Error_string(status, message.data(), &length);
    throw MpiError(
        std::string(operation) + ": "
        + std::string(message.data(), static_cast<std::size_t>(length)));
}

#endif

MpiRuntime::MpiRuntime(int& argc, char**& argv)
{
#if WCNS_HAS_MPI
    int finalized = 0;
    check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized != 0) {
        throw MpiError("MPI has already been finalized");
    }

    int initialized = 0;
    check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");
    if (initialized == 0) {
        check_mpi(
            MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level_),
            "MPI_Init_thread");
        owns_mpi_ = true;
    } else {
        check_mpi(MPI_Query_thread(&thread_level_), "MPI_Query_thread");
    }
    if (thread_level_ < MPI_THREAD_FUNNELED) {
        if (owns_mpi_) {
            MPI_Finalize();
            owns_mpi_ = false;
        }
        throw MpiError("MPI does not provide MPI_THREAD_FUNNELED");
    }

    int raw_rank = 0;
    check_mpi(MPI_Comm_rank(MPI_COMM_WORLD, &raw_rank), "MPI_Comm_rank");
    check_mpi(MPI_Comm_size(MPI_COMM_WORLD, &size_), "MPI_Comm_size");
    check_mpi(
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN),
        "MPI_Comm_set_errhandler");
    rank_ = static_cast<RankId>(raw_rank);
    enabled_ = true;
#else
    static_cast<void>(argc);
    static_cast<void>(argv);
#endif
}

MpiRuntime::~MpiRuntime()
{
#if WCNS_HAS_MPI
    if (!owns_mpi_) {
        return;
    }
    int finalized = 0;
    if (MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0) {
        MPI_Finalize();
    }
#endif
}

void MpiRuntime::barrier() const
{
#if WCNS_HAS_MPI
    check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
#endif
}

Real MpiRuntime::sum(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(
        MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD),
        "MPI_Allreduce sum");
    return global_value;
#else
    return local_value;
#endif
}

Real MpiRuntime::min(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(
        MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD),
        "MPI_Allreduce min");
    return global_value;
#else
    return local_value;
#endif
}

Real MpiRuntime::max(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(
        MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD),
        "MPI_Allreduce max");
    return global_value;
#else
    return local_value;
#endif
}

bool MpiRuntime::all_true(bool local_value) const
{
#if WCNS_HAS_MPI
    const int local = local_value ? 1 : 0;
    int global = 0;
    check_mpi(
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD),
        "MPI_Allreduce logical and");
    return global != 0;
#else
    return local_value;
#endif
}

bool MpiRuntime::all_equal(std::uint64_t local_value) const
{
#if WCNS_HAS_MPI
    unsigned long long local = static_cast<unsigned long long>(local_value);
    unsigned long long minimum = 0;
    unsigned long long maximum = 0;
    check_mpi(
        MPI_Allreduce(
            &local,
            &minimum,
            1,
            MPI_UNSIGNED_LONG_LONG,
            MPI_MIN,
            MPI_COMM_WORLD),
        "MPI_Allreduce uint64 minimum");
    check_mpi(
        MPI_Allreduce(
            &local,
            &maximum,
            1,
            MPI_UNSIGNED_LONG_LONG,
            MPI_MAX,
            MPI_COMM_WORLD),
        "MPI_Allreduce uint64 maximum");
    return minimum == maximum;
#else
    static_cast<void>(local_value);
    return true;
#endif
}

std::string MpiRuntime::broadcast_string(
    std::string value,
    RankId root) const
{
    if (root < 0 || root >= size_) {
        throw MpiError("broadcast root is outside MPI rank range");
    }
#if WCNS_HAS_MPI
    unsigned long long size = rank_ == root
        ? static_cast<unsigned long long>(value.size()) : 0;
    check_mpi(
        MPI_Bcast(
            &size,
            1,
            MPI_UNSIGNED_LONG_LONG,
            static_cast<int>(root),
            MPI_COMM_WORLD),
        "MPI_Bcast string size");
    if (size > static_cast<unsigned long long>(std::numeric_limits<int>::max())) {
        throw MpiError("broadcast string exceeds MPI int count range");
    }
    if (rank_ != root) value.resize(static_cast<std::size_t>(size));
    if (size > 0) {
        check_mpi(
            MPI_Bcast(
                value.data(),
                static_cast<int>(size),
                MPI_CHAR,
                static_cast<int>(root),
                MPI_COMM_WORLD),
            "MPI_Bcast string data");
    }
#else
    static_cast<void>(root);
#endif
    return value;
}

} // namespace wcns
