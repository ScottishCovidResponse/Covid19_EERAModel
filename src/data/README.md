# C++ data management library for Covid19 RAMP

## Introduction RAMP

Unified data management and workflow
Zulip link:


## Third-party library selection
+ poco for networking: the only lib I have used
+ json.hpp header only C++ json: extensible to decode/encode user types
+ HDF5: result packing

### Install dependency on Ubuntu 18.04

` sudo apt install libopenssl-dev, libpoco-dev, libhdf5-dev`

### RHEL 7
Centos versions are older than Ubuntu, sometimes versions cause trouble.
poco-devel-1.6.1-3.el7.x86_64.rpm
hdf5-devel-1.8.12-11.el7.x86_64.rpm
`yum install epel-release -y && yum install poco-devel openssl-devel cmake3 hdf5-devel -y`

This project's cmake can download and build latest HDF5 lib in local build folder. 

#### MacOS
homebrew maybe is the way to go

## Data flow 

### model data flow
-> download json parameters/dataset (by restful GET API) from central store

-> parse parameters into model-specific config file (for the time being)/C++ data structure as model input (later) 
-> prepare dataset into model-specific input file (for the time being)/C++ data structure as model input (later)  

-> model-specific simulation to generate log and result files
-> model-specific post-processing into images, report files

-> packed all into one HDF5 file
-> upload to central store

It is assume model needs 2 kinds of inputs
+ small textual config file to tweak model behavior
+ large dataset files

Model author will focus on model simulation code and judge the quality of result

### EERA specific

2 Inputs adn 1 Output:  IO.h
 `EERAModel::ModelInputParameters ReadParametersFromFile()`
 `EERAModel::InputObservations ReadObservationsFromFiles();`
 `WriteOutputsToFiles()`
Will have their DataFlow version, a macro switch can select. currently IO interface with local files will be maintained as requested by model author.

example inputs can be found in "test/regression/data/"

My understanding is "parameter.ini" does not needs any paramter from remote data store. While the other "Observations" could be retrieved from remote data store.

Input dataset `Observations` can be several 2D datasets in HDF5 or TABLEs in database
All data types are defined in ``ModelTypes.h
```c++
struct InputObservations {
	std::vector<std::vector<int>> cases;
	std::vector<std::vector<int>> deaths;
	std::vector<std::vector<double>> age_pop;
	std::vector<std::vector<double>> waifw_norm;
	std::vector<std::vector<double>> waifw_home;
	std::vector<std::vector<double>> waifw_sdist;
	std::vector<std::vector<double>> cfr_byage;
	std::vector<std::vector<double> > pf_pop;
};
// each field map to a CSV files
```
Currently, `Observations` are read from local csv files in a model root_dir
`std::string(ROOT_DIR)+"/data/waifw_home.csv"`

## Code structure

###  build with cmake

Currently, it is not clear where this source folder should locate.
Only 2 demo cpp to test out platform cmake configuration. 




