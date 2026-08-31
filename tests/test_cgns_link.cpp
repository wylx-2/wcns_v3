#include "test_support.hpp"

#include <cgnslib.h>

// 验收固定 CGNS 库能够被当前构建正确链接和调用。
void test_cgns_link()
{
    WCNS_REQUIRE(CGNS_VERSION == 4400);
    WCNS_REQUIRE(cg_get_error() != nullptr);
}
