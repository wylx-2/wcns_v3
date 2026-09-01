cmake_minimum_required(VERSION 3.20)

project(wcns VERSION 0.1.0 LANGUAGES C CXX)

option(WCNS_BUILD_TESTS "Build WCNS unit tests" ON)
option(WCNS_ENABLE_CGNS "Enable CGNS mesh I/O" ON)
option(WCNS_ENABLE_MPI "Enable MPI parallel execution" OFF)

add_library(wcns_core INTERFACE)
add_library(wcns::core ALIAS wcns_core)

target_include_directories(
    wcns_core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(wcns_core INTERFACE cxx_std_20)

if(MSVC)
    target_compile_options(wcns_core INTERFACE /W4 /permissive-)
else()
    target_compile_options(wcns_core INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_library(wcns_physics STATIC
    src/physics/source_terms.cpp
    src/physics/thermodynamics.cpp
)
add_library(wcns::physics ALIAS wcns_physics)
target_link_libraries(wcns_physics PUBLIC wcns::core)

add_library(wcns_mesh STATIC
    src/mesh/algorithm_profile.cpp
    src/mesh/conservation_weights.cpp
    src/mesh/geometry_halo.cpp
    src/mesh/halo_exchange.cpp
    src/mesh/high_order_metrics.cpp
    src/mesh/linear_operators.cpp
    src/mesh/metrics.cpp
    src/mesh/structured_mesh.cpp
    src/mesh/topology.cpp
)
add_library(wcns::mesh ALIAS wcns_mesh)
target_link_libraries(wcns_mesh PUBLIC wcns::core wcns::physics)

set(WCNS_MPI_TARGET "")
if(WCNS_ENABLE_MPI)
    if(WIN32 AND MINGW)
        find_path(
            WCNS_INTEL_MPI_INCLUDE_DIR
            mpi.h
            HINTS
                "$ENV{I_MPI_ROOT}/include"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/include"
        )
        find_file(
            WCNS_INTEL_MPI_LIBRARY
            impi.lib
            HINTS
                "$ENV{I_MPI_ROOT}/lib/release"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/lib/release"
        )
        find_program(
            WCNS_MPIEXEC_EXECUTABLE
            NAMES mpiexec mpiexec.exe
            HINTS
                "$ENV{I_MPI_ROOT}/bin"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin"
        )
        if(NOT WCNS_INTEL_MPI_INCLUDE_DIR OR NOT WCNS_INTEL_MPI_LIBRARY)
            message(FATAL_ERROR "WCNS_ENABLE_MPI requires an Intel MPI installation")
        endif()
        add_library(wcns_intel_mpi UNKNOWN IMPORTED)
        set_target_properties(
            wcns_intel_mpi
            PROPERTIES
                IMPORTED_LOCATION "${WCNS_INTEL_MPI_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${WCNS_INTEL_MPI_INCLUDE_DIR}"
        )
        set(WCNS_MPI_TARGET wcns_intel_mpi)
    else()
        find_package(MPI REQUIRED COMPONENTS CXX)
        set(WCNS_MPI_TARGET MPI::MPI_CXX)
        set(WCNS_MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}")
    endif()
endif()

add_library(wcns_parallel STATIC
    src/parallel/block_distribution.cpp
    src/parallel/distributed_topology.cpp
    src/parallel/halo_exchanger.cpp
    src/parallel/mpi_runtime.cpp
)
add_library(wcns::parallel ALIAS wcns_parallel)
target_link_libraries(wcns_parallel PUBLIC wcns::mesh)
if(WCNS_ENABLE_MPI)
    target_link_libraries(wcns_parallel PUBLIC ${WCNS_MPI_TARGET})
    target_compile_definitions(wcns_parallel PUBLIC WCNS_HAS_MPI=1)
    if(MINGW)
        target_compile_definitions(wcns_parallel PUBLIC "__int64=long long")
    endif()
else()
    target_compile_definitions(wcns_parallel PUBLIC WCNS_HAS_MPI=0)
endif()

if(WCNS_ENABLE_CGNS)
    include(FetchContent)
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
    endif()

    set(CGNS_ENABLE_HDF5 OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_64BIT OFF CACHE BOOL "" FORCE)
    set(CGNS_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(CGNS_BUILD_TESTING OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        wcns_cgns_source
        URL ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cgns/CGNS-4.4.0.zip
        URL_HASH SHA256=72AD499936089102C1D77E3CDD78A623D61D1C6AD94C4DF3B23F36DB9DDFAD6A
    )
    FetchContent_MakeAvailable(wcns_cgns_source)

    add_library(wcns_cgns INTERFACE)
    add_library(wcns::cgns ALIAS wcns_cgns)
    target_link_libraries(wcns_cgns INTERFACE CGNS::cgns-static)
    if(MINGW)
        target_compile_definitions(wcns_cgns INTERFACE "__int64=long long")
    endif()

    add_library(wcns_io STATIC
        src/io/cgns_reader.cpp
    )
    add_library(wcns::io ALIAS wcns_io)
    target_link_libraries(wcns_io PUBLIC wcns::mesh PRIVATE wcns::cgns)
endif()

add_library(wcns_solver STATIC
    src/solver/boundary_conditions.cpp
    src/solver/euler.cpp
    src/solver/euler_solver.cpp
    src/solver/inviscid_flux.cpp
    src/solver/inviscid_wcns_solver.cpp
    src/solver/physical_boundary.cpp
    src/solver/riemann_solver.cpp
    src/solver/source_operator.cpp
    src/solver/spatial_operator.cpp
    src/solver/time_integrator.cpp
    src/solver/transport_model.cpp
    src/solver/viscous_flux.cpp
    src/solver/wcns_reconstruction.cpp
)
add_library(wcns::solver ALIAS wcns_solver)
target_link_libraries(wcns_solver PUBLIC wcns::parallel wcns::physics)

if(WCNS_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()
