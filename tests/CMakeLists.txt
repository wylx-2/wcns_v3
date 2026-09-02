add_executable(wcns_unit_tests
    test_main.cpp
    test_index.cpp
    test_array3d.cpp
    test_field.cpp
    test_geometry_operators.cpp
    test_geometry_halo.cpp
    test_high_order_metrics.cpp
    test_topology_field.cpp
    test_structured_block.cpp
    test_cgns_link.cpp
    test_conservation_weights.cpp
    test_metrics.cpp
    test_topology.cpp
    test_distribution.cpp
    test_euler.cpp
    test_thermodynamics.cpp
    test_wcns.cpp
    test_spatial_operator.cpp
    test_source_terms.cpp
    test_physical_boundary.cpp
    test_inviscid_reconstruction.cpp
    test_inviscid_flux.cpp
    test_stage_j_sources.cpp
    test_viscous_flux.cpp
    test_viscous_gradient.cpp
    test_viscous_boundary.cpp
    test_case_config.cpp
    test_structured_partition.cpp
    test_flow_initializer.cpp
)

target_link_libraries(
    wcns_unit_tests
    PRIVATE wcns::core wcns::mesh wcns::parallel wcns::solver wcns::runtime
)
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
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_one_sided_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_unknown_donor_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_partition_2d.cgns
    )
    set_tests_properties(wcns.generate_test_cgns PROPERTIES FIXTURES_SETUP cgns_files)

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
        @ONLY
    )
    add_test(
        NAME wcns.run.freestream.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
    )
    set_tests_properties(
        wcns.run.freestream.serial
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )

    add_executable(wcns_cgns_reader_tests
        test_cgns_reader.cpp
    )
    target_link_libraries(
        wcns_cgns_reader_tests
        PRIVATE wcns::io wcns::mesh wcns::parallel
    )
    add_test(
        NAME wcns.cgns_reader
        COMMAND wcns_cgns_reader_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_invalid_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_one_sided_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_unknown_donor_2d.cgns
    )
    set_tests_properties(wcns.cgns_reader PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_halo_exchange_tests test_halo_exchange.cpp)
    target_link_libraries(
        wcns_halo_exchange_tests
        PRIVATE wcns::io wcns::parallel
    )
    add_test(
        NAME wcns.halo_exchange.serial
        COMMAND wcns_halo_exchange_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.halo_exchange.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_euler_solver_tests test_euler_solver.cpp)
    target_link_libraries(
        wcns_euler_solver_tests
        PRIVATE wcns::io wcns::parallel wcns::solver
    )
    add_test(
        NAME wcns.euler_solver.serial
        COMMAND wcns_euler_solver_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.euler_solver.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_geometry_plan_tests test_geometry_plan_mpi.cpp)
    target_link_libraries(
        wcns_geometry_plan_tests
        PRIVATE wcns::io wcns::parallel wcns::mesh wcns::solver
    )
    add_test(
        NAME wcns.geometry_plan.serial
        COMMAND wcns_geometry_plan_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.geometry_plan.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    if(WCNS_ENABLE_MPI)
        add_test(
            NAME wcns.run.freestream.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
        )
        set_tests_properties(
            wcns.run.freestream.2
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
        foreach(rank_count IN ITEMS 4 8)
            add_test(
                NAME wcns.run.freestream.${rank_count}
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
            )
            set_tests_properties(
                wcns.run.freestream.${rank_count}
                PROPERTIES FIXTURES_REQUIRED cgns_files
            )
        endforeach()
        add_test(
            NAME wcns.halo_exchange.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_halo_exchange_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.halo_exchange.2 PROPERTIES FIXTURES_REQUIRED cgns_files)

        add_test(
            NAME wcns.euler_solver.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_euler_solver_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.euler_solver.2 PROPERTIES FIXTURES_REQUIRED cgns_files)

        add_test(
            NAME wcns.geometry_plan.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_geometry_plan_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.geometry_plan.2 PROPERTIES FIXTURES_REQUIRED cgns_files)
    endif()
endif()

add_test(NAME wcns.unit COMMAND wcns_unit_tests)

add_executable(wcns_mpi_runtime_tests test_mpi_runtime.cpp)
target_link_libraries(wcns_mpi_runtime_tests PRIVATE wcns::parallel)
add_test(NAME wcns.mpi_runtime.serial COMMAND wcns_mpi_runtime_tests)

add_executable(wcns_stage_j_solver_tests test_inviscid_wcns_solver_mpi.cpp)
target_link_libraries(
    wcns_stage_j_solver_tests
    PRIVATE wcns::parallel wcns::solver wcns::mesh
)
add_test(NAME wcns.stage_j_solver.serial COMMAND wcns_stage_j_solver_tests)

add_executable(wcns_stage_k_solver_tests test_viscous_wcns_solver_mpi.cpp)
target_link_libraries(
    wcns_stage_k_solver_tests
    PRIVATE wcns::parallel wcns::solver wcns::mesh
)
add_test(NAME wcns.stage_k_solver.serial COMMAND wcns_stage_k_solver_tests)

if(WCNS_ENABLE_MPI)
    if(NOT WCNS_MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "MPI tests require WCNS_MPIEXEC_EXECUTABLE")
    endif()
    add_test(
        NAME wcns.mpi_runtime.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2 $<TARGET_FILE:wcns_mpi_runtime_tests>
    )
    add_test(
        NAME wcns.stage_j_solver.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
            $<TARGET_FILE:wcns_stage_j_solver_tests>
    )
    add_test(
        NAME wcns.stage_k_solver.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
            $<TARGET_FILE:wcns_stage_k_solver_tests>
    )
endif()
