cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME Covid19EERAModel)

project(${PROJECT_NAME} VERSION 0.10.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-DROOT_DIR=\\\"${CMAKE_SOURCE_DIR}\\\" -DVERSION=\\\"${PROJECT_VERSION}\\\" ")

# Enable or Disable Clang Tidy in Build
option(CLANG_TIDY "Enable Clang Tidy in Build" OFF)

# Enable or Disable Coverage with LCov During Compilation
option(CODE_COVERAGE "Enable code coverage with LCov" OFF)

# Allow user to override default of running build in Debug Mode
if(CMAKE_BUILD_TYPE STREQUAL "")
	set(CMAKE_BUILD_TYPE Debug)
endif()

message("----- EERA Model Build Options -----")
message(STATUS "Clang Tidy: ${CLANG_TIDY}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Code Coverage: ${CODE_COVERAGE}")
message("------------------------------------")

if(CLANG_TIDY)
	if(APPLE)
    # Disable Clang-Tidy for macOS (as not available by default)
		message(WARNING "Clang Tidy option not currently supported on macOS, ignoring option")
	else()
		set(CLANG_TIDY_FLAGS -checks=*,-*default-arguments-calls*,-llvm-include-order*,-modernize-use-trailing-return-type*,-readability-isolate-declaration*,-google-runtime-references*,-fuchsia-*,-llvm-header-guard* -header-filter=${CMAKE_SOURCE_DIR}/src,-${CMAKE_CURRENT_SOURCE_DIR}/src/tclap)
		set(CMAKE_CXX_CLANG_TIDY clang-tidy ${CLANG_TIDY_FLAGS})
	endif(APPLE)
endif()
if(CODE_COVERAGE)
	message(STATUS "Including Code Coverage")
	set(COVERAGE_COMPILER_FLAGS "-g -fprofile-arcs -ftest-coverage")
	set (CMAKE_CXX_FLAGS "-std=gnu++11 -Wall -Wextra ${COVERAGE_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
endif()
	

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include(external/tclap.cmake)
add_subdirectory(src)
add_subdirectory(test)

if(CODE_COVERAGE)
	find_program( LCOV_PATH  NAMES lcov lcov.bat lcov.exe lcov.perl)
	find_program( GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat )
	find_program( GCOV_PATH gcov )
	
	set(LCOV_EXCLUDES "build/*" "/usr*")
	
	add_custom_target(coverage
	
	        # Cleanup lcov
	        COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} -directory . -b ${PROJECT_SOURCE_DIR} --zerocounters
	
	        # Create baseline to make sure untouched files show up in the report
		COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} -c -i -d . -b ${PROJECT_SOURCE_DIR} -o ${PROJECT_NAME}.base
	
	        # Run executables
		COMMAND ${CMAKE_SOURCE_DIR}/build/bin/Covid19EERAModel ${COV_TEST_FLAGS} -s original  -m inference
		COMMAND ${CMAKE_SOURCE_DIR}/build/bin/Covid19EERAModel-unit_tests
	
	        # Capturing lcov counters and generating report
		COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} --directory . -b ${PROJECT_SOURCE_DIR} --capture --output-file ${PROJECT_NAME}.capture
	
	        # add baseline counters
		COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} -a ${PROJECT_NAME}.base -a ${PROJECT_NAME}.capture --output-file ${PROJECT_NAME}.total
	
	        # filter collected data to final coverage report and merge outputs
		COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} --remove ${PROJECT_NAME}.total ${LCOV_EXCLUDES} --output-file ${PROJECT_NAME}.info
	
	        # Generate HTML output
		
		COMMAND ${GENHTML_PATH} -o ${PROJECT_NAME}_coverage ${PROJECT_NAME}.info
	
	        # Set output files as GENERATED (will be removed on 'make clean')
	        BYPRODUCTS
	            ${PROJECT_NAME}.base
	            ${PROJECT_NAME}.capture
	            ${PROJECT_NAME}.total
	            ${PROJECT_NAME}.info
	            ${PROJECT_NAME}  # report directory
	
	        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	        VERBATIM # Protect arguments to commands
	        COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
	)
endif()
