# Covid19_EERAModel
![Covid19EERAModel](https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=dev)

Simple COVID-19 simulation model with ABC-smc inference

## Project structure
 * `data/` Working input data directory: used to store inputs to a model run
 * `doc/` Documentation associated with the model
 * `outputs/` Working outputs data directory: used to store outputs from a model run
 * `scripts/` Utility scripts for setting up and executing model runs
 * `src/` Model source code
 * `test/` Model test code

## Dependencies
 * CMake (>= 3.10)
 * GSL (GNU Scientific Library)
 * PkgConfig
 * Threads
 
## Build
The build follows the normal CMake procedure. From the root project directory:
```
mkdir build
cd build
cmake ..
make
```

## Run the model
Following build, the model executable is `build/bin/Covid19EERAModel`. It presently requires that it 
is run from the top-level project directory, where the working input data files are located at 
relative path `./data`.

## Code Documentation Site

Code documentation generated using Doxygen and Code Coverage reports can be found [here](https://scottishcovidresponse.github.io/Covid19_EERAModel/).

## Automated Code Formatting

As part of GitHub actions `clang-format-10` is run on the source code, this ensures consistency between code files without each developer having to worry
about following a convention. Settings are given in the `.clang-format` file.

## Tests

### Regression tests
The regression tests can be found in `test/regression`. Each run uses a fixed seed value, fixed inputs, and a reference set of output data files. A regression test consists of:

* Copying the input data from `test/regression/runN/data` to the working data directory `data`
* Running the model executable `./build/bin/Covid19EERAModel`
* Compare the model outputs in `outputs` with the reference outputs in `test/regression/runN/outputs`

There are two sets of regression tests: one set which use the original model structure, and another 
set which use the Irish epidemiological structure. The former are regression tests 1-6; the latter
are tests 7-12.

The regression tests can be run automatically by running the script `scripts/RunRegressionTests.sh` 
from the top-level roject directory. Each test will be run consecutively, and the script will provide
a summary of success/failures. The script takes the first and last tests to run as arguments i.e.
to run tests 4 through 9:
```
$ ./scripts/RunRegressionTests 4 9
```

**Note:** The regression tests are an aid to refactoring with confidence: they should not be considered
confirmation of the code's correctness. The reference outputs are updates periodically based on 
changes in the core model logic.

### Unit tests
The unit tests can be found in `test/unit`. They are built using the Google Test unit-testing framework. CMake automatically downloads and builds GTest as an external project, so it is not required to have GTest installed on the build system.

Following build, the unit test executable is `build/bin/Covid19EERAModel-unit_tests`.

### Code Coverage
Code coverage is now checked by `lcov` as part of the GitHub actions Ubuntu GCC workflow, a summary of the coverage percentage being given in the output. In addition percentage coverage for each source file are depicted graphically within the file `coverage-output.pdf` which is generated as a downloadable [artifact](https://help.github.com/en/actions/configuring-and-managing-workflows/persisting-workflow-data-using-artifacts) available from within the GitHub action workflow window.

### Check with CppCheck

As part of the validation procedure source and header files are checked with CppCheck. It is recommended you run this on your code before
pushing to the remote repository, from the repository root directory run:

`cppcheck --language=c++ --std=c++11 <address-of-code-file(s)>`


