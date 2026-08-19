cmake_minimum_required(VERSION 3.20)

project(wcns VERSION 0.1.0 LANGUAGES CXX)

option(WCNS_BUILD_TESTS "Build WCNS unit tests" ON)

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

if(WCNS_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()

