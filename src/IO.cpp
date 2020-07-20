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

class IOException: public std::exception
{
public:
    IOException(const std::string& message) : message_(message) {}
    
    const char* what() const throw() 
    {
        return message_.c_str();
    }

private:
    std::string message_;
};

seed ReadSeedSettings(const std::string& ParamsPath, Utilities::logging_stream::Sptr log)
{
    seed seedSettings;
    
    std::string SettingName("seedmethod");
    std::string SettingsCategory("Seed settings");
    seedSettings.seedmethod = CIniFile::GetValue(SettingName, SettingsCategory, ParamsPath); 
    if (seedSettings.seedmethod == "background") 
    {
        seedSettings.hrp = ReadNumberFromFile<int>("hrp", SettingsCategory, ParamsPath);
    } 
    else
    {
        if (seedSettings.seedmethod != "random") 
        {
            (*log) << "WARNING: unknown seeding method, defaulting to random" << std::endl;
        }

        seedSettings.nseed = ReadNumberFromFile<int>("nseed", SettingsCategory, ParamsPath);
    }

    seedSettings.use_fixed_seed = static_cast<bool>(
        ReadNumberFromFile<int>("use_fixed_seed", SettingsCategory, ParamsPath)
    );
    seedSettings.seed_value = ReadNumberFromFile<int>(
        "seed_value", SettingsCategory, ParamsPath
    );

    return seedSettings;
}

params ReadFixedModelParameters(const std::string& ParamsPath)
{
    params paramlist;
    
    std::string SettingsCategory("Fixed parameters");
    paramlist.T_lat = ReadNumberFromFile<double>("T_lat", SettingsCategory, ParamsPath);
    paramlist.juvp_s = ReadNumberFromFile<double>("juvp_s", SettingsCategory, ParamsPath);
    paramlist.T_inf = ReadNumberFromFile<double>("T_inf", SettingsCategory, ParamsPath);
    paramlist.T_rec = ReadNumberFromFile<double>("T_rec", SettingsCategory, ParamsPath);
    paramlist.T_sym = ReadNumberFromFile<double>("T_sym", SettingsCategory, ParamsPath);
    paramlist.T_hos = ReadNumberFromFile<double>("T_hos", SettingsCategory, ParamsPath);
    paramlist.K = ReadNumberFromFile<int>("K", SettingsCategory, ParamsPath);
    paramlist.inf_asym = ReadNumberFromFile<double>("inf_asym", SettingsCategory, ParamsPath);

    return paramlist;
}

SupplementaryInputParameters ReadSupplementaryParameters(const std::string& ParamsPath,
    Utilities::logging_stream::Sptr log)
{
    SupplementaryInputParameters supplementaryParameters;

    supplementaryParameters.seedlist = ReadSeedSettings(ParamsPath, log);

    return supplementaryParameters;
}

CommonModelInputParameters ReadCommonParameters(const std::string& ParamsPath)
{
    CommonModelInputParameters commonParameters;

    commonParameters.paramlist  = ReadFixedModelParameters(ParamsPath);
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);
    commonParameters.totN_hcw   = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", ParamsPath);

    return commonParameters;
}

InferenceConfig ReadInferenceConfig(const std::string& configDir, Utilities::logging_stream::Sptr log) 
{
    std::string ParamsPath(configDir + "/parameters.ini");
    if (!Utilities::fileExists(ParamsPath)) throw IOException(ParamsPath + ": File not found!");;

    InferenceConfig inferenceConfig;
    
    inferenceConfig.seedlist = ReadSeedSettings(ParamsPath, log);
    inferenceConfig.paramlist = ReadFixedModelParameters(ParamsPath);

    inferenceConfig.herd_id = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);
    inferenceConfig.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", ParamsPath);
    inferenceConfig.tau = ReadNumberFromFile<double>("tau", "Settings", ParamsPath);

    inferenceConfig.prior_pinf_shape1 = ReadNumberFromFile<double>("prior_pinf_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_pinf_shape2 = ReadNumberFromFile<double>("prior_pinf_shape2", "Priors settings", ParamsPath);
    inferenceConfig.prior_phcw_shape1 = ReadNumberFromFile<double>("prior_phcw_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_phcw_shape2 = ReadNumberFromFile<double>("prior_phcw_shape2", "Priors settings", ParamsPath);
    inferenceConfig.prior_chcw_mean = ReadNumberFromFile<double>("prior_chcw_mean", "Priors settings", ParamsPath);
    inferenceConfig.prior_d_shape1 = ReadNumberFromFile<double>("prior_d_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_d_shape2 = ReadNumberFromFile<double>("prior_d_shape2", "Priors settings", ParamsPath);
    inferenceConfig.prior_q_shape1 = ReadNumberFromFile<double>("prior_q_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_q_shape2 = ReadNumberFromFile<double>("prior_q_shape2", "Priors settings", ParamsPath);
    inferenceConfig.prior_lambda_shape1 = ReadNumberFromFile<double>("prior_lambda_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_lambda_shape2 = ReadNumberFromFile<double>("prior_lambda_shape2", "Priors settings", ParamsPath);
    
    inferenceConfig.prior_ps_shape1 = ReadNumberFromFile<double>("prior_ps_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_ps_shape2 = ReadNumberFromFile<double>("prior_ps_shape2", "Priors settings", ParamsPath);
    inferenceConfig.prior_rrd_shape1 = ReadNumberFromFile<double>("prior_rrd_shape1", "Priors settings", ParamsPath);
    inferenceConfig.prior_rrd_shape2 = ReadNumberFromFile<double>("prior_rrd_shape2", "Priors settings", ParamsPath);

    inferenceConfig.nsteps = ReadNumberFromFile<int>("nsteps", "Fit settings", ParamsPath);
    inferenceConfig.kernelFactor = ReadNumberFromFile<double>("kernelFactor", "Fit settings", ParamsPath);
    inferenceConfig.nSim = ReadNumberFromFile<int>("nSim", "Fit settings", ParamsPath);
    inferenceConfig.nParticleLimit = ReadNumberFromFile<int>("nParticLimit", "Fit settings", ParamsPath);

    for (int ii = 1; ii <= inferenceConfig.nsteps; ii++) {
        inferenceConfig.toleranceLimit.push_back(0.0);
    }

    for (int ii = 0; ii < inferenceConfig.nsteps; ii++) {
        std::stringstream KeyName;
        KeyName << "Key" << (ii + 1);
        inferenceConfig.toleranceLimit[ii] = ReadNumberFromFile<double>(KeyName.str(), "Tolerance settings", ParamsPath);
    }

    inferenceConfig.observations = ReadInferenceObservations(configDir, log);

    return inferenceConfig;
}

PredictionConfig ReadPredictionConfig(const std::string& configDir, int index, Utilities::logging_stream::Sptr log)
{
    std::string filePath(configDir + "/parameters.ini");
    if (!Utilities::fileExists(filePath)) throw IOException(filePath + ": File not found!");

    PredictionConfig predictionConfig;
    
    predictionConfig.seedlist = ReadSeedSettings(filePath, log);

    predictionConfig.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", filePath);

    std::string sectionId("Prediction Configuration");    
    predictionConfig.n_sim_steps = ReadNumberFromFile<int>("n_sim_steps",
        sectionId, filePath);

    predictionConfig.index = index;

    std::string parametersFile(configDir + "/posterior_parameters.csv");
    if (!Utilities::fileExists(parametersFile)) throw IOException(parametersFile + ": File not found!");
    
    predictionConfig.posterior_parameters = ReadPosteriorParametersFromFile(parametersFile,
        predictionConfig.index);

    return predictionConfig;
}

ObservationsForInference ReadInferenceObservations(const std::string& configDir, Utilities::logging_stream::Sptr log)
{
    ObservationsForInference observations;
    
    (*log) << "Observations For Inference Config:" << std::endl;

    const std::string scot_data_file = configDir + "/scot_data.csv";
    if (!Utilities::fileExists(scot_data_file)) throw IOException(scot_data_file + ": File not found!");
    
    const std::string scot_deaths_file = configDir + "/scot_deaths.csv";
    if (!Utilities::fileExists(scot_deaths_file)) throw IOException(scot_deaths_file + ": File not found!");

    //Uploading observed disease data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log) << "\t- " << scot_data_file << std::endl;
    observations.cases = Utilities::read_csv<int>(scot_data_file, ',');

    //Uploading observed death data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log) << "\t- " << scot_deaths_file << std::endl;
    observations.deaths = Utilities::read_csv<int>(scot_deaths_file, ',');

    return observations;
}

ObservationsForModels ReadModelObservations(const std::string& configDir, Utilities::logging_stream::Sptr log)
{
    ObservationsForModels observations;

    (*log) << "Observations For Models:" << std::endl;

    const std::string scot_data_file = configDir + "/scot_data.csv";
    if (!Utilities::fileExists(scot_data_file)) throw IOException(scot_data_file + ": File not found!");
    
    const std::string scot_ages_file = configDir + "/scot_age.csv";
    if (!Utilities::fileExists(scot_ages_file)) throw IOException(scot_ages_file + ": File not found!");

    const std::string waifw_norm_file = configDir + "/waifw_norm.csv";
    if (!Utilities::fileExists(waifw_norm_file)) throw IOException(waifw_norm_file + ": File not found!");

    const std::string waifw_home_file = configDir + "/waifw_home.csv";
    if (!Utilities::fileExists(waifw_home_file)) throw IOException(waifw_home_file + ": File not found!");

    const std::string waifw_sdist_file = configDir + "/waifw_sdist.csv";
    if (!Utilities::fileExists(waifw_sdist_file)) throw IOException(waifw_sdist_file + ": File not found!");

    const std::string cfr_byage_file = configDir + "/cfr_byage.csv";
    if (!Utilities::fileExists(cfr_byage_file)) throw IOException(cfr_byage_file + ": File not found!");

    const std::string scot_frail_file = configDir + "/scot_frail.csv";
    if (!Utilities::fileExists(scot_frail_file)) throw IOException(scot_frail_file + ": File not found!");


    //Uploading observed disease data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log) << "\t- " << scot_data_file << std::endl;
    observations.cases = Utilities::read_csv<int>(scot_data_file, ',');

    //Uploading population per age group
    //columns are for each individual Health Borad
    //last column is for Scotland
    //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
    (*log) << "\t- " << scot_ages_file << std::endl;
    observations.age_pop = Utilities::read_csv<double>(scot_ages_file, ',');

    //mean number of daily contacts per age group (overall)	
    (*log) << "\t- " << waifw_norm_file << std::endl;
    observations.waifw_norm = Utilities::read_csv<double>(waifw_norm_file, ',');

    //mean number of daily contacts per age group (home only)
    (*log) << "\t- " << waifw_home_file << std::endl;
    observations.waifw_home = Utilities::read_csv<double>(waifw_home_file, ',');

    //mean number of daily contacts per age group (not school, not work)
    (*log) << "\t- " << waifw_sdist_file << std::endl;
    observations.waifw_sdist = Utilities::read_csv<double>(waifw_sdist_file, ',');

    //Upload cfr by age group
    //col0: p_h: probability of hospitalisation
    //col1: cfr: case fatality ratio
    //col2: p_d: probability of death, given hospitalisation
    //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
    (*log) << "\t- " << cfr_byage_file << std::endl;
    observations.cfr_byage = Utilities::read_csv<double>(cfr_byage_file, ',');

    //Upload frailty probability p_f by age group
    //columns are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
    //rows are for each individual Health Borad
    //last row is for Scotland
    (*log) << "\t- " << scot_frail_file << std::endl;
    observations.pf_pop = Utilities::read_csv<double>(scot_frail_file, ',');

    return observations;
}

std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, int set_selection)
{
    // Temporary matrix to hold data from input file
    std::vector<std::vector<double>> lines;
    char delimiter = ',';
    
    lines = Utilities::read_csv<double>(filePath, delimiter);

    // Select line from input file and store result in another temporary vector
    if (set_selection >= static_cast<int>(lines.size())){
        std::stringstream SetSelectError;
        SetSelectError << "Parameter set selection out of bounds! Please select between 0-" << (lines.size() - 1) << "..." << std::endl;
        throw std::overflow_error(SetSelectError.str());
    }
    std::vector<double> line_select = lines[set_selection];

    auto first = line_select.cbegin() + 1;
    auto last = line_select.cend();
    /** @todo Replace with constant from inference parameter description class */
    const int nPar = 8;
    if ((last - first) != nPar) {
        std::stringstream PosteriorFileFormatError;
        PosteriorFileFormatError << "Please check formatting of posterior parameter input file, 8 parameter values are needed..." << std::endl;
        throw std::runtime_error(PosteriorFileFormatError.str());
    }

    std::vector<double> parameter_sets(first, last);

    return parameter_sets;
}

void WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, 
    const std::vector<particle>& particleList, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
    std::stringstream namefile, namefile_simu, namefile_ends;
    namefile << (outDirPath + "/output_abc-smc_particles_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
    namefile_simu << (outDirPath + "/output_abc-smc_simu_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
    namefile_ends << (outDirPath + "/output_abc-smc_ends_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";		

    std::ofstream output_step (namefile.str().c_str());
    std::ofstream output_simu (namefile_simu.str().c_str());
    std::ofstream output_ends (namefile_ends.str().c_str());
    
    //add the column names for each output list of particles
    output_step << "iterID,nsse_cases,nsse_deaths,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,intro,weight" << std::endl;

    //add the column names for each output list of chosen simulations
    output_simu << "iterID" << "," << "day" << "," << "inc_case" << "," << "inc_death_hospital" << "," << "inc_death" << std::endl;
    
    //add the column names for each output list of the compartment values of the last day of the chosen simulations
    output_ends << "iterID" << "," << "age_group" << "," << "comparts" << "," << "value" << std::endl;		

    // outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
    for (int kk = 0; kk < Nparticle; ++kk) {
        output_step << particleList[kk].iter << ", " << particleList[kk].nsse_cases <<  ", " << particleList[kk].nsse_deaths <<  ", " ;

        for (int var = 0; var < nPar; ++var) {
            output_step << particleList[kk].parameter_set[var] << ", ";
        }
        output_step	<< particleList[kk].weight<< '\n';
        
        for (unsigned int var = 0; var < particleList[kk].simu_outs.size(); ++var) {
            output_simu << particleList[kk].iter << ", " << var << ", " <<  particleList[kk].simu_outs[var] << ", " \
                        <<  particleList[kk].hospital_death_outs[var] << ", " << particleList[kk].death_outs[var] << '\n';
        }
        
        for (unsigned int age = 0; age < particleList[kk].end_comps.size(); ++age) {
            for (unsigned int var = 0; var < particleList[kk].end_comps[0].size(); ++var) {
                output_ends << particleList[kk].iter << ", " << age << ", " << var << ", " <<  particleList[kk].end_comps[age][var] << '\n';
            }
        }
    }
    
    output_step.close();
    output_simu.close();
    output_ends.close();
}

void WritePredictionsToFiles(Status status, std::vector<std::vector<int>>& end_comps, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
    std::stringstream namefile_simu, namefile_ends, namefile_full;
    namefile_simu << (outDirPath + "/output_prediction_simu") << "_" << log->getLoggerTime() << ".txt";
    namefile_ends << (outDirPath + "/output_prediction_ends") << "_" << log->getLoggerTime() << ".txt";
    namefile_full << (outDirPath + "/output_prediction_full") << "_" << log->getLoggerTime() << ".txt";

    std::ofstream output_simu (namefile_simu.str().c_str());
    std::ofstream output_ends (namefile_ends.str().c_str());
    std::ofstream output_full (namefile_full.str().c_str());

    output_simu << "day" << "," << "inc_case" << "," << "inc_death_hospital" << "," << "inc_death" << std::endl;

    output_ends << "age_group" << "," << "comparts" << "," << "value" << std::endl;

    for (unsigned int var = 0; var < status.simulation.size(); ++var) {
        output_simu << var << ", " << status.simulation[var] << ", "\
                    << status.hospital_deaths[var] << ", " << status.deaths[var] << '\n';
    }
    
    for (unsigned int age = 0; age < end_comps.size(); ++age) {
        for (unsigned int var = 0; var < end_comps[0].size(); ++var) {
            output_ends << age << ", " << var << ", " << end_comps[age][var] << '\n';
        }
    }

    /** The first two lines are just markers for age group and compartments
     * similar to the first and second columns of output_ends. These can be 
     * removed if it'll make parsing easier.
     */
    for (unsigned int age = 0; age < end_comps.size(); ++age) {
        for (unsigned int comp = 0; comp < end_comps[age].size(); ++comp) {
            output_full << age;
            if (comp < end_comps[0].size() - 1) { output_full << ", "; }
        }
        output_full << "\t";
    }
    output_full << std::endl;

    for (unsigned int age = 0; age < end_comps.size(); ++age) {
        for (unsigned int comp = 0; comp < end_comps[age].size(); ++comp) {
            output_full << comp;
            if (comp < end_comps[0].size() - 1) { output_full << ", "; }
        }
        output_full << "\t";
    }

    output_full << std::endl;
    

    for (unsigned int var = 0; var < status.pop_array.size(); ++var) {
        std::vector<std::vector<int>> pop_array_compartment_to_vector = Model::compartments_to_vector(status.pop_array[var]);
        for (auto & age: pop_array_compartment_to_vector) {
            for (unsigned int comp = 0; comp < age.size(); ++comp) {
                output_full << age[comp];
                if (comp < age.size() - 1) { output_full << ", "; }
            }
            output_full << "\t";
        }
        if (var < status.pop_array.size() - 1) { output_full << '\n'; }
    }

    output_simu.close();
    output_ends.close();
    output_full.close();
}

void LogFixedParameters(const CommonModelInputParameters& params, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Fixed parameter values]:\n";
    (*log) << "    latent period (theta_l): " << params.paramlist.T_lat <<std::endl;
    (*log) << "    pre-clinical period (theta_i): " << params.paramlist.T_inf <<std::endl;
    (*log) << "    asymptomatic period (theta_r): " << params.paramlist.T_rec <<std::endl;
    (*log) << "    symptomatic period (theta_s): " << params.paramlist.T_sym <<std::endl;
    (*log) << "    hospitalisation stay (theta_h): " << params.paramlist.T_hos <<std::endl;
    (*log) << "    pre-adult probability of symptoms devt (p_s[0]): " << params.paramlist.juvp_s <<std::endl;
    (*log) << "    bed capacity at hospital (K): " << params.paramlist.K <<std::endl;
    (*log) << "    relative infectiousness of asymptomatic (u): " << params.paramlist.inf_asym <<std::endl;
}

void LogRandomiserSettings(const SupplementaryInputParameters& params, unsigned long randomiser_seed, 
    Utilities::logging_stream::Sptr log)
{
    (*log) << "[Randomisation Settings]:\n";
    (*log) << "    Seed type: ";
    (*log) << ((params.seedlist.use_fixed_seed) ? "Fixed" : "Time based") << std::endl;
    (*log) << "    Seed value: " << randomiser_seed << std::endl;
}

void LogSeedSettings(const seed& params, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Disease seeding settings]:" << std::endl;
    (*log) << "    seeding method: "<< params.seedmethod <<  std::endl;
    if (params.seedmethod == "random"){
        (*log) << "    number of seed: " << params.nseed << std::endl;
    } else if(params.seedmethod == "background"){
        (*log) << "    duration of the high risk period (hrp): " << params.hrp << std::endl;
    }
}

void LogPredictionConfig(const PredictionConfig& config, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Prediction Configuration]" << std::endl;
    (*log) << "    n_sim_steps: "       << config.n_sim_steps << std::endl;
    (*log) << "    parameter index: "   << config.index << std::endl;
    (*log) << "    p_inf: "             << config.posterior_parameters[0] << std::endl;
    (*log) << "    p_hcw: "             << config.posterior_parameters[1] << std::endl;
    (*log) << "    c_hcw: "             << config.posterior_parameters[2] << std::endl;
    (*log) << "    d: "                 << config.posterior_parameters[3] << std::endl;
    (*log) << "    q: "                 << config.posterior_parameters[4] << std::endl;
    (*log) << "    p_s: "               << config.posterior_parameters[5] << std::endl;
    (*log) << "    rrd: "               << config.posterior_parameters[6] << std::endl;
    (*log) << "    intro: "             << config.posterior_parameters[7] << std::endl;
}

void LogGitVersionInfo(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Git Versioning]"    << std::endl;
    (*log) << "    Commit SHA: "    << GitMetadata::CommitSHA1() << std::endl;
    (*log) << "    Commit Date: "   << GitMetadata::CommitDate() << std::endl;
    (*log) << "    Tag: "           << GitMetadata::Tag() << std::endl;
    (*log) << "    Uncommitted changes: " << (GitMetadata::AnyUncommittedChanges() ? "Yes" : "No") << std::endl;
}

} // namespace IO
} // namespace EERAModel
