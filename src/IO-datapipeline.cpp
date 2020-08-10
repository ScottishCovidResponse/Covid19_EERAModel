#include "IO-datapipeline.h"
#include "IO.h"
#include "ModelCommon.h"
#include "Utilities.h"
#include "Git.h"
#include "array.hh"

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>

namespace EERAModel {
namespace IO {

IOdatapipeline::IOdatapipeline(string params_path, string model_config, Utilities::logging_stream::Sptr log_stream, string dpconfig_path)
{
    ParamsPath = params_path;
    ModelConfigDir = model_config; // Shouldn't need this if data pipeline active
    log = log_stream;

    datapipelineActive = false;

    std::cout << "ParamsPath: " << ParamsPath << "\n";
    std::cout << "ModelConfig: " << ModelConfigDir << "\n";
    std::cout << "ConfigPath: " << dpconfig_path << "\n";

    if (dpconfig_path != "")
    {
        std::string uri = GitMetadata::URL(); // "https://whatever"; I'm guessing this is the repo for the model
        std::string git_sha = GitMetadata::CommitSHA1(); // And this the version ID, need to find these both

        std::cout << "URI: " << uri << "\n";
        std::cout << "SHA: " << git_sha << "\n";

        // If something goes wrong in the opening of the data pipeline an exception will get thrown
        dp.reset(new DataPipeline(dpconfig_path, uri.c_str(), git_sha.c_str()));
        datapipelineActive = true;
    }
}

CommonModelInputParameters IOdatapipeline::ReadCommonParameters()
{
    CommonModelInputParameters commonParameters;

    if (datapipelineActive)
    {
        commonParameters.paramlist  = ReadFixedModelParameters();
        commonParameters.totN_hcw   = dp->read_estimate("fixed-parameters/total_hcw", "total_hcw");
        commonParameters.day_shut = dp->read_estimate("fixed-parameters/day_shut", "day_shut");
    }
    else
    {
        commonParameters.paramlist  = IO::ReadFixedModelParameters(ParamsPath);
        commonParameters.totN_hcw   = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", ParamsPath);
        commonParameters.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", ParamsPath);
    }
    
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);

    return commonParameters;
}

params IOdatapipeline::ReadFixedModelParameters()
{
    std::cout << "(Datapipeline): ReadFixedModelParameters\n";
    params paramlist;

    paramlist.T_lat = dp->read_estimate("fixed-parameters/T_lat", "T_lat");
    paramlist.juvp_s = dp->read_estimate("fixed-parameters/juvp_s", "juvp_s");
    paramlist.T_inf = dp->read_estimate("fixed-parameters/T_inf", "T_inf");
    paramlist.T_rec = dp->read_estimate("fixed-parameters/T_rec", "T_rec");
    paramlist.T_sym = dp->read_estimate("fixed-parameters/T_sym", "T_sym");
    paramlist.T_hos = dp->read_estimate("fixed-parameters/T_hos", "T_hos");
    paramlist.K = dp->read_estimate("fixed-parameters/K", "K");
    paramlist.inf_asym = dp->read_estimate("fixed-parameters/inf_asym", "inf_asym");

    return paramlist;
}

namespace {
    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        os << "[";
        for (auto &ve : v) {
            os << " \"" << ve << "\"";
        }
        os << " ]";
        return os;
    }
}

ObservationsForModels IOdatapipeline::ReadModelObservations()
{
    if (!datapipelineActive)
    {
        return IO::ReadModelObservations(ModelConfigDir, log);
    }
    else
    {
        ObservationsForModels observations;

        (*log) << "Observations For Models:" << std::endl;

        const std::string scot_data_file = ModelConfigDir + "/scot_data.csv";
        if (!Utilities::fileExists(scot_data_file)) throw IOException(scot_data_file + ": File not found!");
        
        const std::string scot_ages_file = ModelConfigDir + "/scot_age.csv";
        if (!Utilities::fileExists(scot_ages_file)) throw IOException(scot_ages_file + ": File not found!");

        const std::string waifw_norm_file = ModelConfigDir + "/waifw_norm.csv";
        if (!Utilities::fileExists(waifw_norm_file)) throw IOException(waifw_norm_file + ": File not found!");

        const std::string waifw_home_file = ModelConfigDir + "/waifw_home.csv";
        if (!Utilities::fileExists(waifw_home_file)) throw IOException(waifw_home_file + ": File not found!");

        const std::string waifw_sdist_file = ModelConfigDir + "/waifw_sdist.csv";
        if (!Utilities::fileExists(waifw_sdist_file)) throw IOException(waifw_sdist_file + ": File not found!");

        const std::string cfr_byage_file = ModelConfigDir + "/cfr_byage.csv";
        if (!Utilities::fileExists(cfr_byage_file)) throw IOException(cfr_byage_file + ": File not found!");

        const std::string scot_frail_file = ModelConfigDir + "/scot_frail.csv";
        if (!Utilities::fileExists(scot_frail_file)) throw IOException(scot_frail_file + ": File not found!");


        //Uploading observed disease data
        //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
        //rows from 1 are indivudual health board
        //last row is for all of scotland

        // vector or vectors returned
        
        (*log) << "\t- " << scot_data_file << std::endl;
        //population-data/data_for_scotland data
        // Array<double> cases = dp->read_array("population-data/data_for_scotland", "data");
        // std::vector<int> cases_size = cases.size();
        // std::cout << "cases.size() = " << cases_size << "\n";

        observations.cases = Utilities::read_csv<int>(scot_data_file, ',');

        //Uploading population per age group
        //columns are for each individual Health Borad
        //last column is for Scotland
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        (*log) << "\t- " << scot_ages_file << std::endl;
        //population-data/data_for_scotland age
        observations.age_pop = Utilities::read_csv<double>(scot_ages_file, ',');

        //mean number of daily contacts per age group (overall)	
        dparray_to_csv("contact-data/who_acquired_infection_from_whom", "norm", &observations.waifw_norm);

        //mean number of daily contacts per age group (home only)
        (*log) << "\t- " << waifw_home_file << std::endl;
        //contact-data/who_acquired_infection_from_whom/data_for_scotland home
        observations.waifw_home = Utilities::read_csv<double>(waifw_home_file, ',');

        //mean number of daily contacts per age group (not school, not work)
        (*log) << "\t- " << waifw_sdist_file << std::endl;
        //contact-data/who_acquired_infection_from_whom/data_for_scotland sdist
        observations.waifw_sdist = Utilities::read_csv<double>(waifw_sdist_file, ',');

        //Upload cfr by age group
        //col0: p_h: probability of hospitalisation
        //col1: cfr: case fatality ratio
        //col2: p_d: probability of death, given hospitalisation
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        (*log) << "\t- " << cfr_byage_file << std::endl;
        // case-fatality-ratios/by_age/data_for_scotland ? (just a copy of contact-data file currently)
        observations.cfr_byage = Utilities::read_csv<double>(cfr_byage_file, ',');

        //Upload frailty probability p_f by age group
        //columns are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        //rows are for each individual Health Borad
        //last row is for Scotland
        (*log) << "\t- " << scot_frail_file << std::endl;
        // Not available currently
        observations.pf_pop = Utilities::read_csv<double>(scot_frail_file, ',');
        return observations;

    }
}

void IOdatapipeline::dparray_to_csv(
    const std::string& data_product, const std::string& component, std::vector<std::vector<double>> *result) {

    (*log) << "\t- (data pipeline) \"" << data_product << "\", \"" << component << "\"" << std::endl;

    Array<double> input = dp->read_array(data_product, component);
    std::vector<int> array_sizes = input.size();

    if (array_sizes.size() != 2) {
        // Should complain about this... and probably should check the dimensions as matching what is expected... if that matters.
    }

    result->resize(0);

    for (int j = 0; j < array_sizes[1]; ++j) {
        // Construct a new element of isize size
        result->emplace_back(array_sizes[0]);
        auto& row = result->back();

        // Copy the data row
        for (int i = 0; i < array_sizes[0]; ++i) {
            row[i] = input(i, j);
        }
    }

    // Some checking
    std::cout << "Original: \"" << data_product << "\" \"" << component << "\"\n";
    std::cout << "    size() = " << array_sizes << "\n";
    std::cout << "    (0,0)=" << input(0,0) << " (1,0)=" << input(1,0) << "\n";
    std::cout << "    (0,1)=" << input(0,1) << " (1,1)=" << input(1,1) << "\n";

    std::cout << "VoV:\n";
    std::cout << "    size() = [" << result->size() << ", " << (*result)[0].size() << "]\n";
    std::cout << "    (0,0)=" << (*result)[0][0] << " (1,0)=" << (*result)[0][1] << "\n";
    std::cout << "    (0,1)=" << (*result)[1][0] << " (1,1)=" << (*result)[1][1] << "\n";
}

} // namespace IO
} // namespace EERAModel