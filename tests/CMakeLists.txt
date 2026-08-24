add_executable(wcns_unit_tests
    test_main.cpp
    test_index.cpp
    test_array3d.cpp
    test_field.cpp
    test_structured_block.cpp
    test_cgns_link.cpp
    test_metrics.cpp
    test_topology.cpp
)

target_link_libraries(wcns_unit_tests PRIVATE wcns::core wcns::mesh)
if(WCNS_ENABLE_CGNS)
    target_link_libraries(wcns_unit_tests PRIVATE wcns::cgns)

    add_executable(wcns_generate_test_cgns
        ${CMAKE_CURRENT_SOURCE_DIR}/../tools/generate_test_cgns.cpp
    )
    target_link_libraries(wcns_generate_test_cgns PRIVATE wcns::core wcns::cgns)

    add_test(
        NAME wcns.generate_test_cgns
        COMMAND wcns_generate_test_cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_invalid_2d.cgns
    )
    set_tests_properties(wcns.generate_test_cgns PROPERTIES FIXTURES_SETUP cgns_files)

    add_executable(wcns_cgns_reader_tests
        test_cgns_reader.cpp
    )
    target_link_libraries(wcns_cgns_reader_tests PRIVATE wcns::io wcns::mesh)
    add_test(
        NAME wcns.cgns_reader
        COMMAND wcns_cgns_reader_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_invalid_2d.cgns
    )
    set_tests_properties(wcns.cgns_reader PROPERTIES FIXTURES_REQUIRED cgns_files)
endif()

add_test(NAME wcns.unit COMMAND wcns_unit_tests)
