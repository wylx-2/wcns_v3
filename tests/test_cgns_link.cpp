#include "test_support.hpp"

#include <cgnslib.h>

void test_cgns_link()
{
    WCNS_REQUIRE(CGNS_VERSION == 4400);
    WCNS_REQUIRE(cg_get_error() != nullptr);
}

