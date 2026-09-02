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
    test_simulation_driver.cpp
    test_stop_controller.cpp
    test_structured_partition.cpp
    test_flow_initializer.cpp
    test_output_manager.cpp
    test_quantity_registry.cpp
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
    set(WCNS_OUTPUT_SUBDIR output-validation-serial)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_output_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_serial.wcns
        @ONLY
    )
    set(WCNS_OUTPUT_SUBDIR output-validation-mpi)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_output_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_mpi.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_OUTPUT checkpoint-validation)
    set(WCNS_CHECKPOINT_ENABLED true)
    set(WCNS_RESTART_FIELD false)
    set(WCNS_RESTART_LINE "# fresh start")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint_source.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_ENABLED false)
    set(WCNS_RESTART_FIELD true)
    set(WCNS_RESTART_LINE
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/checkpoint-validation/checkpoint-freestream.checkpoint.latest.cgns")
    set(WCNS_CHECKPOINT_OUTPUT restart-validation-serial)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_serial.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_OUTPUT restart-validation-mpi)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_mpi.wcns
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
    add_test(
        NAME wcns.run.output.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_serial.wcns
    )
    set_tests_properties(
        wcns.run.output.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP output_serial_files
            RESOURCE_LOCK wcns_output_validation
    )

    add_executable(wcns_field_output_tests test_field_output.cpp)
    target_link_libraries(wcns_field_output_tests PRIVATE wcns::core wcns::cgns)
    add_test(
        NAME wcns.output_files.serial
        COMMAND wcns_field_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/output-validation-serial
            1
    )
    set_tests_properties(
        wcns.output_files.serial
        PROPERTIES FIXTURES_REQUIRED output_serial_files
    )
    add_executable(wcns_checkpoint_file_tests test_checkpoint_file.cpp)
    target_link_libraries(wcns_checkpoint_file_tests PRIVATE wcns::core wcns::cgns)
    add_test(
        NAME wcns.run.checkpoint_source
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/checkpoint_source.wcns
    )
    set_tests_properties(
        wcns.run.checkpoint_source
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP checkpoint_file
            RESOURCE_LOCK wcns_checkpoint_validation
    )
    add_test(
        NAME wcns.checkpoint_file
        COMMAND wcns_checkpoint_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/checkpoint-validation/checkpoint-freestream.checkpoint.latest.cgns
    )
    set_tests_properties(
        wcns.checkpoint_file
        PROPERTIES FIXTURES_REQUIRED checkpoint_file
    )
    add_test(
        NAME wcns.run.restart.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/restart_serial.wcns
    )
    set_tests_properties(
        wcns.run.restart.serial
        PROPERTIES
            FIXTURES_REQUIRED checkpoint_file
            RESOURCE_LOCK wcns_checkpoint_validation
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
        add_test(
            NAME wcns.run.output.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_mpi.wcns
        )
        set_tests_properties(
            wcns.run.output.2
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP output_mpi_files
                RESOURCE_LOCK wcns_output_validation
        )
        add_test(
            NAME wcns.output_files.2
            COMMAND wcns_field_output_tests
                ${CMAKE_CURRENT_BINARY_DIR}/output-validation-mpi
                2
        )
        set_tests_properties(
            wcns.output_files.2
            PROPERTIES FIXTURES_REQUIRED output_mpi_files
        )
        add_test(
            NAME wcns.run.restart.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/restart_mpi.wcns
        )
        set_tests_properties(
            wcns.run.restart.2
            PROPERTIES
                FIXTURES_REQUIRED checkpoint_file
                RESOURCE_LOCK wcns_checkpoint_validation
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
