file(GLOB_RECURSE src_files "*.h" "*.cpp")

find_package(GSL REQUIRED)

add_executable(${PROJECT_NAME} ${src_files})
add_compile_options(${PROJECT_NAME} ${BUILD_TYPE_COMPILE_FLAGS})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_FLAGS "-Wall -Wextra")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas)
