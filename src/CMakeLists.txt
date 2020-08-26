file(GLOB_RECURSE src_files "*.h" "*.cpp")
list(REMOVE_ITEM src_files ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(GSL REQUIRED)

set(PRE_CONFIGURE_FILE "Git.cpp.in")
set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/Git.cpp")
include(${CMAKE_SOURCE_DIR}/cmake/git_watcher.cmake)
add_library(git STATIC ${POST_CONFIGURE_FILE})
target_include_directories(git PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(git check_git)

set(PRECONFIGURE_DEPENDENCY_FILE "Dependencies.cpp.in")
set(POSTCONFIGURE_DEPENDENCY_FILE "${CMAKE_CURRENT_BINARY_DIR}/Dependencies.cpp")
configure_file(${PRECONFIGURE_DEPENDENCY_FILE} ${POSTCONFIGURE_DEPENDENCY_FILE} @ONLY)
list(APPEND src_files ${POSTCONFIGURE_DEPENDENCY_FILE})

set (PROJECT_LIBS ${PROJECT_NAME}-lib)
add_library(${PROJECT_LIBS} STATIC ${src_files})
target_link_libraries(${PROJECT_LIBS} PUBLIC GSL::gsl GSL::gslcblas)
target_include_directories(${PROJECT_LIBS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_LIBS} PUBLIC ${TCLAP})
target_link_libraries(${PROJECT_LIBS} PUBLIC data_pipeline_lib)
target_link_libraries(${PROJECT_LIBS} PUBLIC Python3::Python)
target_link_libraries(${PROJECT_LIBS} PUBLIC git)


add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_LIBS})
