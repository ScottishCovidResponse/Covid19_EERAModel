cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(covid19-cxx-data-management)

enable_language(CXX)

option(ENABLE_HDF5 "Enable HDF5 support" ON)
option(ENABLE_NETWORKING "Enable HDF5 support" ON)

################################ Poco ###############################
if(ENABLE_NETWORKING)
    find_package(Poco REQUIRED COMPONENTS Foundation Util Net NetSSL Crypto JSON)
    include_directories(${Poco_INCLUDE_DIRS})

    # openssl-dev must be installed, but not needed to find the package
    #find_package(OpenSSL REQUIRED)

    if (UNIX)
        if(ON) # centos, rhel7 may have problem in Poco_LIBRARIES
            set(Poco_LIBRARIES "-lPocoNet -lPocoNetSSL -lPocoFoundation -lPocoJSON -lPocoUtil")
            message(STATUS " manually set Poco_LIBRARIES = ${Poco_LIBRARIES}")
            # to make it compatible with Poco_LIBRARIES of FindPoco.cmake
        else()    # system provided Poco::JSON
            message(STATUS "found poco by cmake: Poco_LIBRARIES = ${Poco_LIBRARIES}")
        endif()
    endif()

    if(Poco_FOUND)  # add executable
        add_executable(data_client_example  "data_client_example.cpp")
        target_link_libraries(data_client_example 
            PUBLIC ${Poco_LIBRARIES}
        )
        target_compile_options(data_client_example PRIVATE -Wall -Wextra -Wno-unused)

    else()
        message(FATAL_ERROR "Poco components is not found")
    endif()

endif()


# adapted from https://github.com/bast/hdf5-cmake-example/blob/master/CMakeLists.txt
#  https://cmake.org/cmake/help/v3.6/module/FindHDF5.html
if(ENABLE_HDF5)
    # HL: high-level API: Table, Image data structure
    # NOTE: CXX not CPP for C++ in cmake lang naming
    # prefer serial if both serial and parallel versions are found
    # setting the variable HDF5_PREFER_PARALLEL to true, before FIND_PACKAGE()
    FIND_PACKAGE(HDF5 COMPONENTS C CXX HL REQUIRED)
    if(HDF5_FOUND)
        message(STATUS "Found HDF5 version ${HDF5_VERSION}, ${HDF5_CXX_LIBRARIES}, ${HDF5_C_LIBRARIES}")
        include_directories(${HDF5_INCLUDE_DIR})
        link_directories(${HDF5_INCLUDE_DIRS})  # it is not installed to default lib folder
        set(_hdf5_libs ${HDF5_C_LIBRARIES})
        list(APPEND _hdf5_libs ${HDF5_CXX_LIBRARIES})
        list(APPEND _hdf5_libs ${HDF5_CXX_HL_LIBRARIES})
    else()
        # we did not find it so we fetch it from the web, configure, compile, and link
        # not tested!!!
        include(ExternalProject)
        set(ExternalProjectCMakeArgs
            -DHDF5_BUILD_CPP_LIB=ON
            )
        # need to check the download link to find out the current/latest version
        set(_hdf5_version "1.10.5")  # Ubuntu 18.04 -> version 1.10.0.1
        # for the download and local built, must link to the static lib, not shared
        file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/hdf5_local/src)
        ExternalProject_Add(hdf5_local
            DOWNLOAD_COMMAND curl http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-${_hdf5_version}.tar.gz | tar xvz -C ${PROJECT_BINARY_DIR}/hdf5_local/src
            CMAKE_ARGS ${ExternalProjectCMakeArgs}
            SOURCE_DIR ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}
            BINARY_DIR ${PROJECT_BINARY_DIR}/hdf5_local/build
            STAMP_DIR  ${PROJECT_BINARY_DIR}/hdf5_local/stamp
            TMP_DIR    ${PROJECT_BINARY_DIR}/hdf5_local/tmp
            INSTALL_COMMAND ""
            )
        include_directories(
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}/src
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}/c++/src
            ${PROJECT_BINARY_DIR}/hdf5_local/build
            )
        set(_hdf5_libs
            ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5_cpp.a
            ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5.a
            -ldl
            )
    endif()
endif()

if(ENABLE_HDF5)
    add_executable(data_hdf5_example
        data_hdf5_example.cpp
    )
    if(NOT HDF5_FOUND)
        add_dependencies(data_hdf5_example hdf5_local)
    endif()
    target_link_libraries(
        data_hdf5_example
        ${_hdf5_libs}
        )
endif()