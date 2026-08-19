add_executable(wcns_unit_tests
    test_main.cpp
    test_index.cpp
    test_array3d.cpp
    test_field.cpp
    test_structured_block.cpp
    test_cgns_link.cpp
)

target_link_libraries(wcns_unit_tests PRIVATE wcns::core)
if(WCNS_ENABLE_CGNS)
    target_link_libraries(wcns_unit_tests PRIVATE wcns::cgns)

    add_executable(wcns_generate_test_cgns
        ${CMAKE_CURRENT_SOURCE_DIR}/../tools/generate_test_cgns.cpp
    )
    target_link_libraries(wcns_generate_test_cgns PRIVATE wcns::cgns)

    add_test(
        NAME wcns.generate_test_cgns
        COMMAND wcns_generate_test_cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
    )
endif()

add_test(NAME wcns.unit COMMAND wcns_unit_tests)
