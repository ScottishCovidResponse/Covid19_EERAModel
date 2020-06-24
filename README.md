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

## Build dependencies
 * CMake (>= 3.10)
 * GSL (GNU Scientific Library)
 * PkgConfig
 * Threads

## Third-party libraries
The project uses Google Test (https://github.com/google/googletest) as it unit testing framework. It
is downloaded and built automatically as part of the project build process. It is not required to be
installed on the host system.

The project uses [TCLAP](http://tclap.sourceforge.net/) for parsing of command line arguments. It is
downloaded and built automatically as part of the project build process. It is not required to be 
installed on the host system. If fetching of the repository by CMake fails, a copy of TCLAP will be 
needed. The contents of the `include/tclap` folder should be placed into `src/tclap` within this repository.

## Build
The build follows the normal CMake procedure. To do an out-of-source build, from the root project
directory:
```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Run the model
Following build, the model executable is `build/bin/Covid19EERAModel`. Its usage is:
```
Brief USAGE: 
   ./Covid19_EERAModel/build/bin/Covid19EERAModel  [-d <string>] -s
                                        <original|irish> [-l <string>] [-m
                                        <inference|prediction>] [--]
                                        [--version] [-h]

For complete USAGE and HELP type: 
   ./Covid19_EERAModel/build/bin/Covid19EERAModel --help
```

If the options for `structure` or `mode` are not set they will have the value `Default` in the 
`[Arguments]` log section, this means that the value will either be obtained from the included 
parameters file (if specified in that file) or given the default value of `inference` for mode and
`original` for structure.

At the present time, the `d` and `l` options are unused and can be omitted.

## Code Documentation Site

Code documentation generated using Doxygen and Code Coverage reports can be found [here](https://scottishcovidresponse.github.io/Covid19_EERAModel/).

## Automated Code Formatting

As part of GitHub actions `clang-format-10` is run on the source code, this ensures consistency between code files without each developer having to worry
about following a convention. Settings are given in the `.clang-format` file.

## Tests

### Regression tests
The regression tests can be found in `test/regression`. Each run uses a fixed seed value, fixed inputs,
and a reference set of output data files. A regression test consists of:

* Copying the input data from `test/regression/runN/data` to the working data directory `data`
* Running the model executable `./build/bin/Covid19EERAModel`
* Compare the model outputs in `outputs` with the reference outputs in `test/regression/runN/outputs`

There are two sets of regression tests: one set which use the original model structure, and another 
set which use the Irish epidemiological structure. The former are regression tests 1-6; the latter
are tests 7-12.

The regression tests can be run automatically by running the script `scripts/RunRegressionTests.sh` 
from the top-level roject directory. Each test will be run consecutively, and on completion the 
script will provide a summary of successes and failures. The script takes the first and last tests
to run as arguments i.e. to run tests 4 through 9, execute the command:
```
$ ./scripts/RunRegressionTests 4 9
```

**Note:** The regression tests are an aid to refactoring with confidence: they should not be considered
confirmation of the code's correctness. The reference outputs are updated periodically based on 
changes in the core model logic.

### Unit tests
The unit tests can be found in `test/unit`. They are built using the Google Test unit-testing framework.
CMake automatically downloads and builds Google Test as an external project, so it is not required to have
Google Test installed on the build system.

Following build, the unit test executable is `build/bin/Covid19EERAModel-unit_tests`.

### Code Coverage
Code coverage is now checked by `lcov` as part of the GitHub actions Ubuntu GCC workflow, a summary of the coverage percentage being given in the output. In addition percentage coverage for each source file are depicted graphically within the file `coverage-output.pdf` which is generated as a downloadable [artifact](https://help.github.com/en/actions/configuring-and-managing-workflows/persisting-workflow-data-using-artifacts) available from within the GitHub action workflow window.

### Check with CppCheck

As part of the validation procedure source and header files are checked with CppCheck. It is recommended you run this on your code before
pushing to the remote repository, from the repository root directory run:

`cppcheck --language=c++ --std=c++11 <address-of-code-file(s)>`
