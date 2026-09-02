if(NOT DEFINED PROGRAM OR NOT DEFINED CONFIG OR NOT DEFINED EXPECTED_EXIT)
    message(FATAL_ERROR "PROGRAM, CONFIG and EXPECTED_EXIT are required")
endif()

execute_process(
    COMMAND "${PROGRAM}" --config "${CONFIG}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
)
if(NOT result EQUAL EXPECTED_EXIT)
    message(FATAL_ERROR
        "expected exit ${EXPECTED_EXIT}, got ${result}\nstdout:\n${output}\nstderr:\n${error}")
endif()
