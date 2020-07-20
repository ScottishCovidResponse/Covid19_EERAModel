if (NOT DATA_PIPELINE)
  message (FATAL_ERROR "Variable DATA_PIPELINE must be defined!")
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

execute_process(
    COMMAND bash "-c" "python3 -c 'import os; import pybind11; print(os.path.dirname(pybind11.__file__)+\"/include\")'"
    OUTPUT_VARIABLE PYBIND11_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (PYBIND11_INCLUDE_DIRS)
    message("Found pybind11 include directories at ${PYBIND11_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "Failed to find pybind11 include directories!")
endif()

add_library(data_pipeline_lib SHARED IMPORTED)
set_target_properties(data_pipeline_lib PROPERTIES
  IMPORTED_LOCATION "${DATA_PIPELINE}/bindings/cpp/build/libdatapipeline.a"
  INTERFACE_INCLUDE_DIRECTORIES "${DATA_PIPELINE}/bindings/cpp;${PYBIND11_INCLUDE_DIRS};${Python3_INCLUDE_DIRS}"
)