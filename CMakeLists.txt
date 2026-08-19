cmake_minimum_required(VERSION 3.20)

project(wcns VERSION 0.1.0 LANGUAGES C CXX)

option(WCNS_BUILD_TESTS "Build WCNS unit tests" ON)
option(WCNS_ENABLE_CGNS "Enable CGNS mesh I/O" ON)

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

add_library(wcns_mesh STATIC
    src/mesh/metrics.cpp
)
add_library(wcns::mesh ALIAS wcns_mesh)
target_link_libraries(wcns_mesh PUBLIC wcns::core)

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
    target_link_libraries(wcns_io PUBLIC wcns::core PRIVATE wcns::cgns)
endif()

if(WCNS_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()
