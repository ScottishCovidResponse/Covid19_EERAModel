# Covid19_EERAModel
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
 * CppCheck
 
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

## Tests

### Regression tests
The regression tests can be found in `test/regression`. Each run uses a fixed seed value, fixed inputs, and a reference set of output data files. A regression test consists of:

* Copying the input data from `test/regression/runN/data` to the working data directory `data`
* Running the model executable `./build/bin/Covid19EERAModel`
* Compare the model outputs in `outputs` with the reference outputs in `test/regression/runN/outputs`

The regression tests can be run automatically by running the script `scripts/RunRegressionTests.sh` from the top-level roject directory. Each test will be run consecutively, and the script will provide a summary of success/failures.

**Note:** The regression tests are an aid to refactoring with confidence: they should not be considered confirmation of the code's correctness. The reference outputs are based on the last baseline model code (Version 0.3.2.4 at the time of writing).

### Unit tests
The unit tests can be found in `test/unit`. They are built using the Google Test unit-testing framework. CMake automatically downloads and builds GTest as an external project, so it is not required to have GTest installed on the build system.

Following build, the unit test executable is `build/bin/Covid19EERAModel-unit_tests`.

### Check with CppCheck

As part of the validation procedure source and header files are checked with CppCheck. It is recommended you run this on your code before
pushing to the remote repository, from the repository root directory run:

`cppcheck --language=c++ --std=c++11 <address-of-code-file(s)>`


