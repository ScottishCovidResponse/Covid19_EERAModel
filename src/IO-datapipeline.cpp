#include "IO-datapipeline.h"
#include "IO.h"
#include "ModelCommon.h"
#include "Utilities.h"
#include "Git.h"

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

        const std::string cfr_byage_file = ModelConfigDir + "/cfr_byage.csv";
        if (!Utilities::fileExists(cfr_byage_file)) throw IOException(cfr_byage_file + ": File not found!");

        const std::string scot_frail_file = ModelConfigDir + "/scot_frail.csv";
        if (!Utilities::fileExists(scot_frail_file)) throw IOException(scot_frail_file + ": File not found!");


        //Uploading observed disease data
        //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
        //rows from 1 are indivudual health board
        //last row is for all of scotland
        dparray_to_csv<int>("population-data/data_for_scotland", "data", &observations.cases);
        // auto cases = Utilities::read_csv<int>(scot_data_file, ',');

        // TODO: this is indexed by herd_id, and the data file has a titles row that the data pipeline doesn't
        // so need to do something about that. Fix this properly, but for now...

        observations.cases.insert(observations.cases.begin(), std::vector<int>(observations.cases.front().size()));

        int v = -1;
        for (auto& el : observations.cases.front()) {
            el = v;
            v++;
        }

        //Uploading population per age group
        //columns are for each individual Health Borad
        //last column is for Scotland
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        dparray_to_csv<double>("population-data/data_for_scotland", "age", &observations.age_pop);
        // auto age_pop = Utilities::read_csv<double>(scot_ages_file, ',');

        // std::cout << observations.age_pop << "\n\n";
        // std::cout << age_pop << "\n";

        //mean number of daily contacts per age group (overall)	
        //contact-data/who_acquired_infection_from_whom/data_for_scotland norm
        dparray_to_csv<double>("contact-data/who_acquired_infection_from_whom", "norm", &observations.waifw_norm);

        //mean number of daily contacts per age group (home only)
        //contact-data/who_acquired_infection_from_whom/data_for_scotland home
        dparray_to_csv<double>("contact-data/who_acquired_infection_from_whom", "home", &observations.waifw_home);
        // auto waifw_home = Utilities::read_csv<double>(waifw_home_file, ',');

        //mean number of daily contacts per age group (not school, not work)
        //contact-data/who_acquired_infection_from_whom/data_for_scotland sdist
        dparray_to_csv<double>("contact-data/who_acquired_infection_from_whom", "sdist", &observations.waifw_sdist);
        // auto waifw_sdist = Utilities::read_csv<double>(waifw_sdist_file, ',');

        //Upload cfr by age group
        //col0: p_h: probability of hospitalisation
        //col1: cfr: case fatality ratio
        //col2: p_d: probability of death, given hospitalisation
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        // prob_hosp_and_cfr/data_for_scotland cfr_byage
        dptable_to_csv<double>("prob_hosp_and_cfr/data_for_scotland", "cfr_byage", &observations.cfr_byage);
        //observations.cfr_byage = Utilities::read_csv<double>(cfr_byage_file, ',');

        //std::cout << observations.cfr_byage << "\n\n";
        //std::cout << cfr_byage << "\n";

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

} // namespace IO
} // namespace EERAModel