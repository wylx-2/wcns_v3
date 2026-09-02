#include "test_support.hpp"

#include <wcns/parallel/mpi_runtime.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>

// 验收 MPI 生命周期、rank/size 和全局归约封装。
int main(int argc, char** argv)
{
    try {
        wcns::MpiRuntime mpi(argc, argv);
        WCNS_REQUIRE(mpi.rank() >= 0);
        WCNS_REQUIRE(mpi.rank() < mpi.size());
        WCNS_REQUIRE(mpi.size() >= 1);
        const auto rank_value = static_cast<wcns::Real>(mpi.rank() + 1);
        const auto expected_sum
            = static_cast<wcns::Real>(mpi.size() * (mpi.size() + 1) / 2);
        WCNS_REQUIRE(mpi.sum(rank_value) == expected_sum);
        WCNS_REQUIRE(mpi.min(rank_value) == 1.0);
        WCNS_REQUIRE(mpi.max(rank_value) == static_cast<wcns::Real>(mpi.size()));
        WCNS_REQUIRE(mpi.all_true(mpi.rank() < mpi.size()));
        WCNS_REQUIRE(mpi.all_equal(42));
        WCNS_REQUIRE(
            mpi.size() == 1
            || !mpi.all_equal(static_cast<std::uint64_t>(mpi.rank())));
        const auto message = mpi.broadcast_string(
            mpi.rank() == 0 ? std::string("wcns-config") : std::string(),
            0);
        WCNS_REQUIRE(message == "wcns-config");
        const auto gathered = mpi.gather_reals({
            static_cast<wcns::Real>(mpi.rank()), 2.0});
        if (mpi.rank() == 0) {
            WCNS_REQUIRE(gathered.size()
                == static_cast<std::size_t>(2 * mpi.size()));
            for (int rank = 0; rank < mpi.size(); ++rank) {
                WCNS_REQUIRE_NEAR(
                    gathered[static_cast<std::size_t>(2 * rank)],
                    static_cast<wcns::Real>(rank),
                    1.0e-15);
            }
        } else {
            WCNS_REQUIRE(gathered.empty());
        }
        std::vector<wcns::Real> scatter_values;
        std::vector<std::size_t> scatter_counts;
        if (mpi.rank() == 0) {
            scatter_counts.assign(static_cast<std::size_t>(mpi.size()), 2);
            for (int rank = 0; rank < mpi.size(); ++rank) {
                scatter_values.push_back(static_cast<wcns::Real>(rank));
                scatter_values.push_back(3.0);
            }
        }
        const auto scattered = mpi.scatter_reals(
            scatter_values, scatter_counts);
        WCNS_REQUIRE(scattered.size() == 2);
        WCNS_REQUIRE_NEAR(
            scattered[0], static_cast<wcns::Real>(mpi.rank()), 1.0e-15);
        WCNS_REQUIRE_NEAR(scattered[1], 3.0, 1.0e-15);
        mpi.barrier();
        if (mpi.rank() == 0) {
            std::cout << "MPI runtime test passed with " << mpi.size() << " ranks\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
