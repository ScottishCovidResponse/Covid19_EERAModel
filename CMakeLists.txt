cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME Covid19EERAModel)

project(${PROJECT_NAME} VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)
set (CMAKE_CXX_FLAGS "-std=gnu++11 -Wall -Wextra -DROOT_DIR=\\\"${CMAKE_SOURCE_DIR}\\\"")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include(CodeCoverage)
append_coverage_compiler_flags()

setup_target_for_coverage_lcov(NAME coverage EXECUTABLE ${CMAKE_BINARY_DIR}/Covid19EERAModel)

add_subdirectory(src)
add_subdirectory(test)
