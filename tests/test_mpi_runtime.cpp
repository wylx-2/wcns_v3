#include "test_support.hpp"

#include <wcns/parallel/mpi_runtime.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>

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

