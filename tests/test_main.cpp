#include <wcns/core/types.hpp>

#include <cstdlib>
#include <iostream>
#include <type_traits>

void test_index();

int main()
{
    try {
        static_assert(std::is_same_v<wcns::Real, double>);
        static_assert(std::is_signed_v<wcns::BlockId>);
        static_assert(std::is_signed_v<wcns::RankId>);

        if (wcns::invalid_block_id >= 0 || wcns::invalid_rank_id >= 0) {
            std::cerr << "invalid identifiers must be negative\n";
            return EXIT_FAILURE;
        }

        test_index();

        std::cout << "stage A tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
