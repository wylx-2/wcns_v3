add_executable(wcns_stage_a_tests
    test_main.cpp
    test_index.cpp
    test_array3d.cpp
)

target_link_libraries(wcns_stage_a_tests PRIVATE wcns::core)

add_test(NAME wcns.stage_a COMMAND wcns_stage_a_tests)
