# Covid19_EERAModel
Simple COVID-19 simulation model with ABC-smc inference

## Project structure
 * `doc/` Documentation associated with the model
 * `data/` Test data files
 * `src/` Model source code
 * `test/` Model test code

## Dependencies
 * CMake (>= 3.10)
 * GSL (GNU Scientific Library)
 * PkgConfig
 * Threads
 * GoogleTest

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
is run from the top-level project directory, where the test data and parameters are located at 
relative paths `./data` and `./data/parameters.ini`, respectively.

## Run the tests
Following build, the unit test executable is `build/bin/Covid19EERAModel-unit_tests`.

## Known issues
    * Issues in initiating simulations (can require multiple restarts)
    * Multithreading via OpenMP not working on macOS
    * Because of memory allocation limits for structures, the model is not able to output the total
      number of deaths (community + hospital).
