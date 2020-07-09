include(ExternalProject)
ExternalProject_Add(tclap_src
    GIT_REPOSITORY https://git.code.sf.net/p/tclap/code
    GIT_TAG 07472612522e4da6b96f589ccd1d20607c28c2b8
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(tclap_src source_dir)
set(TCLAP_INCLUDE_DIRS ${source_dir}/include)

set(TCLAP tclap)
add_library(${TCLAP} UNKNOWN IMPORTED)
set_target_properties(${TCLAP} PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${TCLAP_INCLUDE_DIRS})
