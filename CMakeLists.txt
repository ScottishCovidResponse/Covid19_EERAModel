cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME Covid19EERAModel)

project(${PROJECT_NAME} VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(src)
add_subdirectory(test)