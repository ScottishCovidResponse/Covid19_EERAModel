file(GLOB_RECURSE src_files "*.h" "*.cpp")
list(REMOVE_ITEM src_files ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(GSL REQUIRED)

# Fetch the TCLAP library from the git repository and place it
# into the source folder under src/tclap
IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tclap )
    execute_process( COMMAND git archive --remote=git://git.code.sf.net/p/tclap/code HEAD:include -o tclap.zip )
    execute_process( COMMAND unzip -d ${CMAKE_CURRENT_SOURCE_DIR} tclap.zip )
ENDIF()

set (PROJECT_LIBS ${PROJECT_NAME}-lib)
add_library(${PROJECT_LIBS} SHARED ${src_files})
target_link_libraries(${PROJECT_LIBS} PUBLIC GSL::gsl GSL::gslcblas)
target_include_directories(${PROJECT_LIBS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_LIBS})
