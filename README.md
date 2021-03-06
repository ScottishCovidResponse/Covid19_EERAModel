# Covid19_EERAModel
![Covid19EERAModel](https://github.com/ScottishCovidResponse/Covid19_EERAModel/workflows/Covid19EERAModel/badge.svg?branch=dev)

Simple COVID-19 simulation model with ABC-smc inference

## Branches and releases
The most recently released code is always located on the `master` branch. The `dev` is the main development branch, containing the most recently developed pre-release features. 

Releases are versioned using the [semantic versioning](www.semver.org) scheme. Releases are tagged with the release number, in the form MAJOR.MINOR.PATCH.

## Project structure
 * `cmake/` Project CMake scripts
 * `data/` Working input data directory: used to store inputs to a model run
 * `doc/` Documentation associated with the model
 * `external/` External third-party libraries
 * `outputs/` Working outputs data directory: used to store outputs from a model run
 * `scripts/` Utility scripts for setting up and executing model runs
 * `src/` Model source code
 * `test/` Model test code

## Build dependencies
The project has the following dependencies, which must be satisfied by installing the relevant packages on the host system. Where the version is prefixed by '~', this version has been tested, but it is likely that earlier versions will also work.

 * CMake (>= 3.10)
 * GNU Scientific Library (~2.4)
 * PkgConfig (~0.29.0)
 * System threads package (pthread or other)

## Third-party libraries
The project uses [Google Test](https://github.com/google/googletest) as its unit testing framework. It is downloaded and built automatically as part of the project build process. It is not required to be installed on the host system.

The project uses [TCLAP](http://tclap.sourceforge.net/) for parsing of command line arguments. It is downloaded and built automatically as part of the project build process. It is not required to be installed on the host system.

The project makes use of a [CMake-based version-tracking tool](https://github.com/andrew-hardin/cmake-git-version-tracking) for encoding Git repository information in the project binary files. It was developed by Andrew Hardin (https://github.com/andrew-hardin), and is licensed uder an MIT license.

The project makes use of an [INI file parser](https://www.codeproject.com/Articles/8342/CIniFile-Class-for-C-A-robust-cross-platform-INI-f) to read values from parameter files. It is used under the terms of the license [here](https://www.codeproject.com/info/cpol10.aspx).

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
Following build, the model executable is located at the path`build/bin/Covid19EERAModel` relative to the root of the project directory. Its usage is:
```

$ ./bin/build/Covid19EERAModel  -m <inference|prediction>]
                    -s <original|irish|irish2>
                    [-i <integer>] 
                    [-d <string>]
                    [-l <string>] [--]
                    [--version] [-h]

```
The two mandatory options are `-s` for the model structure, and `-m` for the run mode. Omission of either of these options will cause the run to terminate with an error message. The index option specifies the parameter set that will be used in prediction mode: it is unused in inference mode.

At the present time, the `-d` and `-l` options are unused by the code and can be omitted.

The command line options used on a given model run are logged in the output log file (see below) for traceability and future reference.

### Inputs
The model requires a number of input files to run, in addition to the command line arguments. Input files must be placed in the `data` directory, and must be named according to the table below. The required contents of each file are described in more detail underneath the table. 

**Examples of each of the files below can be found in the `data/example` directory.**

| File name        | Description           | Usage|
| ------------- |:-------------:|:-------------:|
| parameters.ini      | General parameters configuring model runs | All |
| cfr_byage.csv | Case Fatality Ratios, by age group      | All |
| scot_age.csv | Proportion of health board populations in each age group      | All |
| scot_data.csv | Timeseries of observed disease cases, by health board      | Inference only|
| scot_deaths.csv | Timeseries of observed disease deaths, by health board      | Inference only |
| scot_frail.csv | Probability of frailty, by age group      | All |
| waifw_home.csv | Age Mixing Matrix (Home)| All |
| waifw_norm.csv | Age Mixing Matrix (All contact included)| All |
| waifw_sdist.csv |  Age Mixing Matrix (Social Distancing)| All |
|posterior_parameters.csv| Posterior model parameters| Prediction only|

#### parameters.ini
This file contains general model parameters, in `.ini` format. Parameters are grouped into sections. The required sections and parameters are as follows. This file is used in both inference and prediction runs.

| Section name        | Parameter name           | Parameter type        | Description           |
| ------------- |:-------------:|:-------------:|:-------------:|
| Settings      | shb\_id           | Integer        | Identifier for selected health board (1-15)           |
| Settings       | tau           | Float        |  Time step scale factor           |
| Settings       | nHealthBoards | Int        |  Number of Health Boards (used to validate the input files - not used in a model run)           |
| Settings       | nAgeGroups | Int        |  Number of Age Groups (used to validate the input files - not used in a model run)           |
| Settings       | nCfrCategories | Int        | Number of Case State Categories (used to validate the input files - not used in a model run)          |
| Settings       | nCasesDays | Int        | Number of days that observations are recorded for (used to validate the input files - not used in a model run)           |
| Seed settings       | seedmethod           | String        | Seeding method ("background" or "random")           |
| Seed settings       | nseed           | Integer        | Population seeding number <br>(Random seeding only)           |
| Seed settings       | hrp           | Integer        | High Risk Period in days <br>(Background seeding only)        |
| Seed settings       | use\_fixed\_seed           | Integer        | If 1, use fixed randomiser seed value<br>If 0, use time-based seed           |
| Seed settings       | seed\_value          | Integer        | If use\_fixed\_seed is 1, the randomiser seed value to use|
| Fit settings       | nsteps          | Integer        | Number of inference steps to run|
| Fit settings       | nParticLimit          | Integer        | Maximum number of inference particles to accept in an inference step|
| Fit settings       | nSim          | Integer        | Maximum number of model runs to execute per inference step|
| Fit settings       | kernelFactor          | Float        | Scale factor for inference parameter kernel window|
| Tolerance settings       | Key1..10          | Float        | Tolerance factor for accepting inference particles|
| Fixed parameters        | totN\_hcw          | Integer        | Total number of health care workers in Scotland|
| Fixed parameters        | day\_shut          | Integer        | Time at which lockdown began <br>(days with respect to time series start)|
| Fixed parameters        | T\_lat          | Float        | Mean latent period (days)|
| Fixed parameters        | juvp\_s          | Float        | Probability of juvenile developing symptoms|
| Fixed parameters        | T\_inf          | Float        | Mean asymptomatic period (days)|
| Fixed parameters        | T\_rec          | Float        | Mean time to recovery if symptomatic (days)|
| Fixed parameters        | T\_sym          | Float        | Mean symptomatic period prior to hospitalisation (days)|
| Fixed parameters        | T\_hos          | Float        | Mean hospitalisation stay (days)|
| Fixed parameters        | K          | Integer        | Hospital bed capacity|
| Fixed parameters        | inf\_asym          | Float        | Reduction factor of infectiousness for asymptomatic infectious individuals|
| Priors settings        | prior\_pinf\_shape1          | Float        | Probability of Infection <br>Beta distribution shapre parameter 1|
| Priors settings        | prior\_pinf\_shape2          | Float        | Probability of Infection <br>Beta distribution shape parameter 2 |
| Priors settings        | prior\_phcw\_shape1          | Float        | Probability of Infection (HCW) <br>Beta distribution shapre parameter 1|
| Priors settings        | prior\_phcw\_shape2          | Float        |  Probability of Infection (HCW) <br>Beta distribution shapre parameter 2|
| Priors settings        | prior\_chcw\_mean          | Float        | Mean number of HCW contacts per day<br>Poisson distribution mean|
| Priors settings        | prior\_d\_shape1          | Float        | Proportion of population observing social distancing<br> Beta distribution shape parameter 1|
| Priors settings        | prior\_d\_shape2          | Float        | Proportion of population observing social distancing<br> Beta distribution shape parameter 2|
| Priors settings        | prior\_q\_shape1          | Float        | Proportion of normal contact made by people self-isolating<br> Beta distribution shape parameter 1|
| Priors settings        | prior\_q\_shape2          | Float        | Proportion of normal contact made by people self-isolating<br> Beta distribution shape parameter 2|
| Priors settings        | prior\_ps\_shape1          | Float        | Age-dependent probability of developing symptoms<br> Beta distribution shape parameter 1|
| Priors settings        | prior\_ps\_shape2          | Float        | Age-dependent probability of developing symptoms<br> Beta distribution shape parameter 2|
| Priors settings        | prior\_rrd\_shape1          | Float        | Risk of death if not hospitalised<br> Gamma distribution shape parameter 1|
| Priors settings        | prior\_rrd\_shape2          | Float        |  Risk of death if not hospitalised<br> Gamma distribution shape parameter 2|
| Priors settings        | prior\_lambda\_shape1          | Float        | Background transmission rate<br> Uniform distribution shape parameter 1|
| Priors settings        | prior\_lambda\_shape2          | Float        |  Background transmission rate<br> Uniform distribution shape parameter 1|
| Prediction Configuration        | n\_sim\_steps          | Float        | Number of model iterations <br>(prediction mode only)|

#### cfr_byage.csv
CSV file containing the Case Fatality Ratio, by age group. Each row is a different age group. The four columns are:

* **Column 0**: Probability of hospitalisation
* **Column 1**: Case Fatality Ratio
* **Column 2**: Probability of Death, given hospitalisation
* **Column 3**: Unused

#### scot_age.csv
CSV file containing the proportion of people in each age group, per health board population. Each row corresponds to a different health board, whicle each column is an age group. This file does not include HCW as a distinct age group. The proprotion of HCW in a population is estimated at run time.

#### scot_data.csv, scot_deaths.csv
CSV file containing the timeseries of cases and deaths, per health board. Each row corresponds to a different health board, while ach column is a day in the time series. The first column is the toal population of the health board.

#### scot_frail.csv
CSV file containing the probabilities of frailty for each age group, by health board. Each column is an age group. Each row is a health board, with the exception of the last row, which is for the whole of Scotland.

#### waifw_home.csv, waifw_norm.csv, waifw_sdist.csv
CSV files containing the age mixing matrices for people (1) isolating at home, (2) behaving normally, and (3) socially distancing. 

#### posterior_parameters.csv
CSV file containing batched parameter sets (fixed and inferrred parameters). This file is only used in **prediction mode**.  Its format is:
```
Index,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,lambda, T_lat, juvp_s, T_inf, T_rec, T_sym, T_hos, K, inf_asym
0,0.153532,0.60916,37.9059,0.525139,0.313957,0.787278,0.516736,8.50135E-07,4,0.1, 1.5,11,7,5,2000,1
...
```
Each row in the file contains 17 entries: the first is the index of the row; the following 8 are the inferred posterior parameters; and the remaining 8 are model fixed parameters. The row selected for use in the prediction run will be that specified by the index argument on the command line (see Prediction Mode discussion below).

### Prediction mode
The model can be run in a prediction mode, where a fixed set of parameters is supplied to the model,
and the model is run for a fixed number of simulation steps.

To run the model in prediction mode, set the `-m` switch to prediction:
```
$ .build/bin/Covid19EERAModel -m prediction [-i <integer>]...
```
To configure the prediction run, three main pieces of configuration are required: a posterior parameters
file (described above); a `Prediction Config` settings category in the `parameters.ini` file; and an index as a command line argument. The index should be provided by using the `-i` option on the command line. If it is omitted, it will default to 0.

The `parameters.ini` file must contain a category with the configuration of the prediction run, as below
```
[Prediction Configuration]
n_iterations=100
n_sim_steps=100000

```
The setting `n_iterations` sets the number of model runs which should occur. The setting `n_sim_steps` sets the number of days each individual model run should should be run for. 

When the model is run in prediction mode, all of the above configuration is logged to the terminal and the log file.

### Inference mode
To run the model in inference mode, set the `-m` command line switch to inference:
```
$ .build/bin/Covid19EERAModel -m inference ...
```

## Outputs
The model generates a number of output files on each run. Output files are stored in the `outputs` directory.

Different output files are produced depending on whether the model is run in inference or prediction mode. These files are described in the section below. 

Regardless of the mode in which the model is run, a *log file* is always produced. This file records the same information printed to the terminal while the model is running. The contents of the log file include command line arguments, versioning information, and copies of most of the significant input parameters. The log file is stored in the `outputs/logs` directory.

The log file includes a section listing Git repository version information, of the form
```
[Git Versioning]
    Commit SHA: xxxxxxx
    Commit Date: 
    Tag: 
    Uncommitted changes:
```
Listed are the SHA of the `HEAD` commit, the corresponding commit date, the tag (if any), and a message to say if there are any uncommitted changes in the repository at the time of the last build.

### Inference mode outputs
In inference mode, three types of file are produced in the `outputs` directory. These files have names of the form:

* `output_abc-smc_ends_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt`
* `output_abc-smc_particles_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt`
* `output_abc-smc_simu_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt`

In each case, the placeholder `<n>` will be the inference step number (a separate copy of each file is produced for each inference step), and `<m>` is the ID of the Health Board for which inference is being performed (this corresponds to the `shb_id` parameter in the input `parameters.ini` file).

#### End State File
The `output_abc-smc_ends_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt` (hereafter referred to as the end-state file) records the state of the simulated population at the end of a simulation run. The format is: 
```
iterID,age_group,comparts,value
0, 0, 0, 750169
...
```
where `iterID` is the simulation number, `age_group` is the age\_group number, `comparts` is the epidemiological compartment number, and `value` is the population of that compartment at the end of the simulation run.

The number of each age group is defined in the table below:

|**Age group number**|**Description**|
|----|----|
| 0  | Under 20	|
| 1  | 20-29	|
| 2  | 30-39	|
| 3  | 40-49	|
| 4  | 50-59	|
| 5  | 60-69	|
| 6  | 70+	|
| 7  | Health Care Workers	|

(Note that Health Care Workers are assumed in the model to have behaviour similar to the average of groups between the age of 20 and 59.)

The number of each epidemiological compartment number is given in the table below:

|**ID**|**Compartment number**|**Description**|
|----|----|----|
| S   | 0  | Number of susceptible individuals (not infected).	|
| E   | 1  | Number of infected individuals but not yet infectious (exposed).	|
| E_t | 2  | Number of exposed individuals and tested positive.	|
| I_p | 3  | Number of infected and infectious symptomatic individuals but at pre-clinical stage (show yet no symptoms).	|
| I_t | 4  | Number of tested positive individuals that infectious.	|
| I1  | 5  | Number of infected and infectious asymptomatic individuals: first stage.	|
| I2  | 6  | Number of infected and infectious asymptomatic individuals: second stage. 	|
| I3  | 7  | Number of infected and infectious asymptomatic individuals: third stage. 	|
| I4  | 8  | Number of infected and infectious asymptomatic individuals: last stage. 	|
| I_s1| 9  | Number of infected and infectious symptomatic individuals: first stage.	|
| I_s2| 10 | Number of infected and infectious symptomatic individuals: second stage. 	|
| I_s3| 11 | Number of infected and infectious symptomatic individuals: thrid stage.	|
| I_s4| 12 | Number of infected and infectious symptomatic individuals: last stage. 	|
| H   | 13 | Number of infected individuals that are hospitalised. 	|
| R   | 14 | Number of infected individuals that are recovered from infection.  	|
| D   | 15 | Number of dead individuals due to disease.|

#### Particles file

The `output_abc-smc_particles_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt` (hereafter referred to as the particles file) records the contents of each particle accepted as part of the inference process. A particle consists of a collection of inferred parameter values, and a weighting.

The format is: 
```
iterID,nsse_cases,nsse_deaths,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,intro,weight
0, 0.0866585, 0.126525, 0.157731, 0.425757, 37, 0.433313, 0.431202, 0.783405, 0.796279, 8.78239e-07, 1
...
```
where `iterID` is the simulation number for which the particle was considered, as before. The `weight` is the weighting of the particle. The remaining fields are described in the table below.

|**Particle**|**Description**|
|----|----|
| nsse_cases | Normalised sum of square error for the number of cases. |
| nsse_deaths | Normalised sum of square error for the number of deaths. |
| p_inf | Probability of Infection |
| p_hcw | Probability of Infection (Healthcare Worker) |
| c_hcw | Mean number of Healthcare Worker contacts per day |
| d | Proportion of population observing social distancing |
| q | Proportion of normal contact made by people self-isolating |
| p_s | Age-dependent probability of developing symptoms |
| rrd | Risk of death if not hospitalised |
| lambda | Background transmission rate |

#### Simulation file
The `output_abc-smc_simu_step<n>_shb<m>_dd-mm-yyyy_hh-mm-ss.txt` (hereafter referred to as the simulation file) records the disease incidence at each day in each simulation run. Its format is: 
```
iterID,day,inc_case,inc_death_hospital,inc_death
0, 0, 0, 0, 0
...
```
where `iterID` is the simulation number, `day` is the day number, `inc_cases` is the simulated incidence of cases occurring on that day, `inc_death_hospital` is the simulated incidence of hospital deaths on that day, and `inc_deaths` is the simulated incidence of non-hospital deaths occurring on that day.

### Prediction mode outputs
In prediction mode, two types of file are produced in the `outputs` directory. These files have names of the form:

* `output_prediction_full_dd-mm-yyyy_hh-mm-ss.txt`
* `output_prediction_simu_dd-mm-yyyy_hh-mm-ss.txt`

#### Full file
The `output_prediction_full_dd-mm-yyyy_hh-mm-ss.txt` (hereafter referred to as the full file) records the state of the simulated population on each day in each simulation run. The format is: 
```
iter, day, age_group, S, E, E_t, I_p, I_t, I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D
0, 0, 0, 1266783, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
...
```
where `iter` is the simulation number, `day` is the day number, `age_group` is the age\_group number, and the remaining fields are the population of the relevant compartment on that day in the given simulation run. The compartment names are the same as those described for the inference mode end-state file above. The age group numbers are likewise the same.
 
#### Simulation file
The `output_prediction_simu_dd-mm-yyyy_hh-mm-ss.txt` file (hereafter referred to as the simulation file) records the disease incidence at each day in each simulation run. Its format is identical to that of the inference mode simulation file described above.

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

There are multiple sets of regression tests, which exercise different model structures in both inference and forward prediction modes. The table below lists the configuration of each test:

| Test numbers        | Model           | Mode|
| ------------- |:-------------:|:-------------:|
| 1 - 6      | Original | Inference |
| 7 - 12      | Irish | Inference |
| 13 - 18      | Original | Forward Prediction |
| 19 - 24      | Irish | Forward Prediction |

The regression tests can be run automatically by running the script `scripts/RunRegressionTests.sh` from the top-level roject directory. Each test will be run consecutively, and on completion the script will provide a summary of successes and failures. The script takes the first and last tests to run as arguments i.e. to run tests 4 through 9, execute the command:
```
$ ./scripts/RunRegressionTests 4 9
```
The regression test script automatically configures each run in line with the table above: the user does not need to do this.

**Note:** The regression tests are an aid to refactoring with confidence: they should not be considered confirmation of the code's correctness. The reference outputs are updated periodically based on changes in the core model logic.

### Unit tests
The unit tests can be found in `test/unit`. They are built using the Google Test unit-testing framework. CMake automatically downloads and builds Google Test as an external project, so it is not required to have Google Test installed on the build system.

Following build, the unit test executable is `build/bin/Covid19EERAModel-unit_tests`.

### Code Coverage
Code coverage is now checked by `lcov` as part of the GitHub actions Ubuntu GCC workflow, a summary of the coverage percentage being given in the output. In addition percentage coverage for each source file are depicted graphically within the file `coverage-output.pdf` which is generated as a downloadable [artifact](https://help.github.com/en/actions/configuring-and-managing-workflows/persisting-workflow-data-using-artifacts) available from within the GitHub action workflow window.

### Check with CppCheck

As part of the validation procedure source and header files are checked with CppCheck. It is recommended you run this on your code before pushing to the remote repository, from the repository root directory run:

`cppcheck --language=c++ --std=c++11 <address-of-code-file(s)>`
