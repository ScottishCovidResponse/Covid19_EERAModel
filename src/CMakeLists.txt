file(GLOB_RECURSE src_files "*.h" "*.cpp")
list(REMOVE_ITEM src_files ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(GSL REQUIRED)

set (PROJECT_LIBS ${PROJECT_NAME}-lib)
add_library(${PROJECT_LIBS} SHARED ${src_files})
target_link_libraries(${PROJECT_LIBS} PUBLIC GSL::gsl GSL::gslcblas)
target_include_directories(${PROJECT_LIBS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_LIBS})
