# Covid19_EERAModel
![Covid19EERAModel](https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=dev)

Simple COVID-19 simulation model with ABC-smc inference

## Project structure
 * `data/` Working input data directory: used to store inputs to a model run
 * `doc/` Documentation associated with the model
 * `external/` External third-party libraries
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
The project uses [Google Test](https://github.com/google/googletest) as its unit testing framework. It
is downloaded and built automatically as part of the project build process. It is not required to be
installed on the host system.

The project uses [TCLAP](http://tclap.sourceforge.net/) for parsing of command line arguments. It is
downloaded and built automatically as part of the project build process. It is not required to be 
installed on the host system. 

## Build
The build follows the normal CMake procedure. To do an out-of-source build, from the root project
directory:
```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Running the model
Following build, the model executable is `build/bin/Covid19EERAModel`. Its usage is:
```

$ Covid19EERAModel  -m <inference|prediction>]
                    -s <original|irish|irish2>
                    [-d <string>]
                    [-l <string>] [--]
                    [--version] [-h]

```
The two mandatory options are "-s" for the model structure, and "-m" for the run mode. Omission of
either of these options will cause the run to terminate with an error message.

At the present time, the `-d` and `-l` options are unused by the code and can be omitted.

The command line options supplied are logged in the output log file for future reference.

### Prediction mode
The model can be run in a prediction mode, where a fixed set of parameters is supplied to the model,
and the model is run for a fixed number of simulation steps.

To run the model in prediction mode, set the `-m` switch to prediction:
```
$ .build/bin/Covid19EERAModel -m prediction ...
```
To configure the prediction run, two main pieces of configuration are required: a posterior parameters
file, and a `Prediction Config` category in the `parameters.ini` file.

The posterior parameters file should be named `posterior_parameters.csv` and should be placed in the
`data` directory. Its format should be:
```
Index,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,intro
0,0.153532,0.60916,37.9059,0.525139,0.313957,0.787278,0.516736,8.50135E-07
1,0.12602,0.429026,43.9404,0.277644,0.722916,0.470001,3.41006,6.27917E-07
...
```
Each row in the file contains 9 entries: the first is the index of the row, and the remaining 8 are 
the sets of posterior parameters that can be used by the predication run.

The `parameters.ini` file must contain a category with the configuration of the prediction run, as
below
```
[Prediction Configuration]
posterior_parameter_index=2
n_sim_steps=100000
```
The setting `n_sim_steps` is the number of iterations the model should be run for. The setting 
`posterior_parameter_index` is the index of the selected posterior parameter set from the 
`posterior_parameters.csv` file.

When the model is run in prediction mode, all of the above configurations are logged to the terminal
and the log file.

### Inference mode
To run the model in inference mode, set the `-m` switch to inference:
```
$ .build/bin/Covid19EERAModel -m inference ...
```
To configure the inference run, two main pieces of configuration must be supplied:
  * A set of data files containing observations and population/age-group descriptions
  * A `parameters.ini` file

Examples of these files can be found in the regression tests directories: `test/regression/run<N>/data`.

## Code Documentation Site

Code documentation generated using Doxygen and Code Coverage reports can be found [here](https://scottishcovidresponse.github.io/Covid19_EERAModel/).

## Automated Code Formatting

As part of GitHub actions `clang-format-10` is run on the source code, this ensures consistency between code files without each developer having to worry about following a convention. Settings are given in the `.clang-format` file.

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
