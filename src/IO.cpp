#include "IO.h"
#include "ModelCommon.h"
#include "Utilities.h"
#include "Git.h"
#include "Dependencies.h"

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>

namespace EERAModel {
namespace IO {

void ImportConsistencyCheck(const std::string& filePath, const unsigned int& axisLength, const unsigned int& expectedValue, const std::string& axisID)
{
    if (axisLength != expectedValue) {
        std::stringstream IOMessage;
        IOMessage << "Error in : " << filePath <<
        "\n Number of " << axisID << ": " << axisLength <<
        "\n Expected number of " << axisID << ": " << expectedValue << std::endl;
        throw IOException(IOMessage.str());
    }

}

ValidationParameters ImportValidationParameters(const std::string& filePath)
{
    ValidationParameters parameters;

    parameters.nHealthBoards = ReadNumberFromFile<int>("nHealthBoards", "Settings", filePath);
    parameters.nAgeGroups = ReadNumberFromFile<int>("nAgeGroups", "Settings", filePath);
    parameters.nCfrCategories = ReadNumberFromFile<int>("nCfrCategories", "Settings", filePath);
    parameters.nCasesDays = ReadNumberFromFile<int>("nCasesDays", "Settings", filePath);

    return parameters;
}

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

InferenceConfig ReadInferenceConfig(const std::string& configDir, Utilities::logging_stream::Sptr log, const CommonModelInputParameters& commonParameters) 
{
    std::string ParamsPath(configDir + "/parameters.ini");

    InferenceConfig inferenceConfig;
    
    ReadLocalInferenceConfig(ParamsPath, log, commonParameters, &inferenceConfig);    

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

    inferenceConfig.observations = ReadInferenceObservations(configDir, log);

    return inferenceConfig;
}


void ReadLocalInferenceConfig(
    const std::string& ParamsPath, Utilities::logging_stream::Sptr log, const CommonModelInputParameters& commonParameters,
    InferenceConfig *inferenceConfig)
{
    inferenceConfig->seedlist = ReadSeedSettings(ParamsPath, log);
    //inferenceConfig->paramlist = ReadFixedModelParameters(ParamsPath);
    inferenceConfig->paramlist = commonParameters.paramlist;
    inferenceConfig->herd_id = commonParameters.herd_id;
    inferenceConfig->day_shut = commonParameters.day_shut;
    inferenceConfig->tau = ReadNumberFromFile<double>("tau", "Settings", ParamsPath);

    inferenceConfig->nsteps = ReadNumberFromFile<int>("nsteps", "Fit settings", ParamsPath);
    inferenceConfig->kernelFactor = ReadNumberFromFile<double>("kernelFactor", "Fit settings", ParamsPath);
    inferenceConfig->nSim = ReadNumberFromFile<int>("nSim", "Fit settings", ParamsPath);
    inferenceConfig->nParticleLimit = ReadNumberFromFile<int>("nParticLimit", "Fit settings", ParamsPath);

    for (int ii = 1; ii <= inferenceConfig->nsteps; ii++) {
        inferenceConfig->toleranceLimit.push_back(0.0);
    }

    for (int ii = 0; ii < inferenceConfig->nsteps; ii++) {
        std::stringstream KeyName;
        KeyName << "Key" << (ii + 1);
        inferenceConfig->toleranceLimit[ii] = ReadNumberFromFile<double>(KeyName.str(), "Tolerance settings", ParamsPath);
    }
}

PredictionConfig ReadPredictionConfig(const std::string& configDir, int index, Utilities::logging_stream::Sptr log, const CommonModelInputParameters& commonParameters)
{
    std::string filePath(configDir + "/parameters.ini");

    PredictionConfig predictionConfig;

    ReadLocalPredictionConfig(filePath, index, log, commonParameters, &predictionConfig);

    std::string parametersFile(configDir + "/posterior_parameters.csv");
    if (!Utilities::fileExists(parametersFile)) {
        std::stringstream error_message;
        error_message << "Cannot locate posterior parameters file at " << parametersFile << std::endl;
        throw std::runtime_error(error_message.str());
    }

    // The posterior parameters are columns 0 to 7; the fixed parameters are columns 8 to 15
    std::vector<double> modelParameters = ReadPredictionParametersFromFile(parametersFile, predictionConfig.index);
    predictionConfig.posterior_parameters = std::vector<double>(modelParameters.begin(), modelParameters.begin() + 8);
    predictionConfig.fixedParameters.T_lat      = modelParameters[8];
    predictionConfig.fixedParameters.juvp_s     = modelParameters[9];
    predictionConfig.fixedParameters.T_inf      = modelParameters[10];
    predictionConfig.fixedParameters.T_rec      = modelParameters[11];
    predictionConfig.fixedParameters.T_sym      = modelParameters[12];
    predictionConfig.fixedParameters.T_hos      = modelParameters[13];
    predictionConfig.fixedParameters.K          = static_cast<int>(modelParameters[14]);
    predictionConfig.fixedParameters.inf_asym   = modelParameters[15];

    return predictionConfig;
}

void ReadLocalPredictionConfig(
    const std::string& filePath, int index, Utilities::logging_stream::Sptr log,
    const CommonModelInputParameters& commonParameters, PredictionConfig *predictionConfig)
{
    predictionConfig->seedlist = ReadSeedSettings(filePath, log);

    predictionConfig->day_shut = commonParameters.day_shut;

    std::string sectionId("Prediction Configuration");    
    predictionConfig->n_sim_steps = ReadNumberFromFile<int>("n_sim_steps",
        sectionId, filePath);
    predictionConfig->n_iterations = ReadNumberFromFile<int>("n_iterations",
        sectionId, filePath);

    predictionConfig->index = index;
}

ObservationsForInference ReadInferenceObservations(const std::string& configDir, Utilities::logging_stream::Sptr log)
{
    ObservationsForInference observations;
    
    (*log) << "Observations For Inference Config:" << std::endl;

    const std::string scot_data_file = configDir + "/scot_data.csv";
    const std::string scot_deaths_file = configDir + "/scot_deaths.csv";

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

    const std::string settings_file = configDir + "/parameters.ini";

    ValidationParameters validationParams = ImportValidationParameters(settings_file);
    int nHealthBoards = validationParams.nHealthBoards;
    int nCasesDays = validationParams.nCasesDays;

    unsigned int cases_rows = observations.cases.size();
    unsigned int cases_cols = observations.cases[0].size();

    ImportConsistencyCheck(scot_data_file, cases_rows, (nHealthBoards + 1), "rows");
    ImportConsistencyCheck(scot_data_file, cases_cols, nCasesDays, "columns");

    unsigned int deaths_rows = observations.deaths.size();
    unsigned int deaths_cols = observations.deaths[0].size();

    ImportConsistencyCheck(scot_deaths_file, deaths_rows, (nHealthBoards + 1), "rows");
    ImportConsistencyCheck(scot_deaths_file, deaths_cols, nCasesDays, "columns");

    return observations;
}

ObservationsForModels ReadModelObservations(const std::string& configDir, Utilities::logging_stream::Sptr log)
{
    ObservationsForModels observations;

    const std::string scot_data_file = configDir + "/scot_data.csv";
    const std::string scot_ages_file = configDir + "/scot_age.csv";
    const std::string waifw_norm_file = configDir + "/waifw_norm.csv";
    const std::string waifw_home_file = configDir + "/waifw_home.csv";
    const std::string waifw_sdist_file = configDir + "/waifw_sdist.csv";
    const std::string cfr_byage_file = configDir + "/cfr_byage.csv";
    const std::string settings_file = configDir + "/parameters.ini";

    ValidationParameters validationParameters = ImportValidationParameters(settings_file);
    int nHealthBoards = validationParameters.nHealthBoards;
    int nAgeGroups = validationParameters.nAgeGroups;
    int nCfrCategories = validationParameters.nCfrCategories;
    int nCasesDays = validationParameters.nCasesDays;

    //Uploading observed disease data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log) << "\t- " << scot_data_file << std::endl;
    observations.cases = Utilities::read_csv<int>(scot_data_file, ',');

    unsigned int cases_rows = observations.cases.size();
    unsigned int cases_cols = observations.cases[0].size();

    ImportConsistencyCheck(scot_data_file, cases_rows, (nHealthBoards + 1), "rows");
    ImportConsistencyCheck(scot_data_file, cases_cols, nCasesDays, "columns");

    //Uploading population per age group
    //columns are for each individual Health Borad
    //last column is for Scotland
    //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
    (*log) << "\t- " << scot_ages_file << std::endl;
    observations.age_pop = Utilities::read_csv<double>(scot_ages_file, ',');

    unsigned int age_pop_rows = observations.age_pop.size();
    unsigned int age_pop_cols = observations.age_pop[0].size();

    ImportConsistencyCheck(scot_ages_file, age_pop_rows, nHealthBoards, "rows");
    ImportConsistencyCheck(scot_ages_file, age_pop_cols, (nAgeGroups - 1), "columns");


    //mean number of daily contacts per age group (overall)	
    (*log) << "\t- " << waifw_norm_file << std::endl;
    observations.waifw_norm = Utilities::read_csv<double>(waifw_norm_file, ',');

    unsigned int waifw_norm_rows = observations.waifw_norm.size();
    unsigned int waifw_norm_cols = observations.waifw_norm[0].size();

    ImportConsistencyCheck(waifw_norm_file, waifw_norm_rows, nAgeGroups, "rows");
    ImportConsistencyCheck(waifw_norm_file, waifw_norm_cols, nAgeGroups, "columns");

    //mean number of daily contacts per age group (home only)
    (*log) << "\t- " << waifw_home_file << std::endl;
    observations.waifw_home = Utilities::read_csv<double>(waifw_home_file, ',');

    unsigned int waifw_home_rows = observations.waifw_home.size();
    unsigned int waifw_home_cols = observations.waifw_home[0].size();

    ImportConsistencyCheck(waifw_home_file, waifw_home_rows, nAgeGroups, "rows");
    ImportConsistencyCheck(waifw_home_file, waifw_home_cols, nAgeGroups, "columns");

    //mean number of daily contacts per age group (not school, not work)
    (*log) << "\t- " << waifw_sdist_file << std::endl;
    observations.waifw_sdist = Utilities::read_csv<double>(waifw_sdist_file, ',');

    unsigned int waifw_sdist_rows = observations.waifw_sdist.size();
    unsigned int waifw_sdist_cols = observations.waifw_sdist[0].size();
    
    ImportConsistencyCheck(waifw_sdist_file, waifw_sdist_rows, nAgeGroups, "rows");
    ImportConsistencyCheck(waifw_sdist_file, waifw_sdist_cols, nAgeGroups, "columns");

    //Upload cfr by age group
    //col0: p_h: probability of hospitalisation
    //col1: cfr: case fatality ratio
    //col2: p_d: probability of death, given hospitalisation
    //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
    (*log) << "\t- " << cfr_byage_file << std::endl;
    observations.cfr_byage = Utilities::read_csv<double>(cfr_byage_file, ',');

    unsigned int cfr_rows = observations.cfr_byage.size();
    unsigned int cfr_cols = observations.cfr_byage[0].size();

    ImportConsistencyCheck(cfr_byage_file, cfr_rows, nAgeGroups, "rows");
    ImportConsistencyCheck(cfr_byage_file, cfr_cols, nCfrCategories, "columns");

    return observations;
}

std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, int set_selection)
{
    // Temporary matrix to hold data from input file
    std::vector<std::vector<double>> lines;
    char delimiter = ',';
    
    lines = Utilities::read_csv<double>(filePath, delimiter, true);

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

std::vector<double> ReadPredictionParametersFromFile(const std::string& filePath, int index)
{
    // Temporary matrix to hold data from input file
    std::vector<std::vector<double>> lines;
    char delimiter = ',';
    
    lines = Utilities::read_csv<double>(filePath, delimiter, true);

    // Select line from input file and store result in another temporary vector
    if (index >= static_cast<int>(lines.size())){
        std::stringstream SetSelectError;
        SetSelectError << "Parameter set selection out of bounds! Please select between 0-" << (lines.size() - 1) << "..." << std::endl;
        throw std::overflow_error(SetSelectError.str());
    }
    std::vector<double> line_select = lines[index];

    // Skip the index column
    auto first = line_select.cbegin() + 1;
    auto last = line_select.cend();
    
    // There are 16 columns in the file: 8 posterior parameters + 8 fixed parameters
    const int nPar = 16;
    if ((last - first) != nPar) {
        std::stringstream PosteriorFileFormatError;
        PosteriorFileFormatError << "Please check formatting of posterior parameter input file, 16 parameter values are needed..." << std::endl;
        throw std::runtime_error(PosteriorFileFormatError.str());
    }

    return std::vector<double>(first, last);
}

bool ReadBoolFromFile(std::string SettingName, std::string SettingCategory, const std::string& filePath) 
{
	std::string SettingValue = CIniFile::GetValue(SettingName, SettingCategory, filePath);

	char* endptr = nullptr;
	bool Value = false;

	/* Convert to upper case */
    SettingValue = Utilities::toUpper(SettingValue);

	if (SettingValue == "TRUE" || SettingValue == "T" || SettingValue == "1")
	{
		Value = true;
	}
	else if (SettingValue == "FALSE" || SettingValue == "F" || SettingValue == "0")
	{
		Value = false;
	}
	else
	{
		std::stringstream SettingParseError;
		SettingParseError << std::endl;
		SettingParseError << "Invalid value in Parameter File: " << filePath.c_str() <<  std::endl;
		SettingParseError << "Category: " << SettingCategory.c_str() << std::endl;
		SettingParseError << "Setting: " << SettingName.c_str() << std::endl;
		SettingParseError << "Value: " << SettingValue.c_str() << std::endl;
		throw std::runtime_error(SettingParseError.str());
	}
	return Value;
}

std::string ReadStringFromFile(std::string SettingName, std::string SettingCategory, const std::string& filePath) 
{
	return CIniFile::GetValue(SettingName, SettingCategory, filePath);
}

bool ExistsInFile(std::string SettingName, std::string SettingCategory, const std::string& filePath)
{
	return CIniFile::RecordExists(SettingName, SettingCategory, filePath);
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
    WriteInferenceParticlesHeader(output_step);

    //add the column names for each output list of chosen simulations
    WriteSimuHeader(output_simu);
    
    //add the column names for each output list of the compartment values of the last day of the chosen simulations
    WriteInferenceEndsHeader(output_ends);	

    // outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
    for (int kk = 0; kk < Nparticle; ++kk) {
        WriteInferenceParticlesRow(output_step, kk, particleList[kk]);
        
        for (unsigned int var = 0; var < particleList[kk].simu_outs.size(); ++var) {
            WriteSimuRow(output_simu, particleList[kk].iter, var , particleList[kk].simu_outs[var],
                        particleList[kk].hospital_death_outs[var], particleList[kk].death_outs[var]);
        }
        
        for (unsigned int age = 0; age < particleList[kk].end_comps.size(); ++age) {
            WriteInferenceEndsRow(output_ends, particleList[kk].iter, age, particleList[kk].end_comps[age]);
        }
    }
    
    output_step.close();
    output_simu.close();
    output_ends.close();
}

void WritePredictionsToFiles(std::vector<Status> statuses, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
    std::stringstream namefile_simu, namefile_full;
    namefile_simu << (outDirPath + "/output_prediction_simu") << "_" << log->getLoggerTime() << ".txt";
    namefile_full << (outDirPath + "/output_prediction_full") << "_" << log->getLoggerTime() << ".txt";

    std::ofstream output_simu(namefile_simu.str());
    std::ofstream output_full(namefile_full.str());

    WriteSimuHeader(output_simu);
    for (unsigned int iter = 0; iter < statuses.size(); ++iter) {
        const Status& status = statuses[iter];
        
        for (unsigned int day = 0; day < status.simulation.size(); ++day) {
            WriteSimuRow(output_simu, iter, day, status.simulation[day],
                status.hospital_deaths[day], status.deaths[day]);
        }
    }

    WritePredictionFullHeader(output_full);
    for (unsigned int iter = 0; iter < statuses.size(); ++iter) {
        const Status& status = statuses[iter];
        const auto& pop_array = status.pop_array;
        
        for (unsigned int day = 0; day < pop_array.size(); ++day) {
            const auto& age_groups = pop_array[day];
            
            for (unsigned int age = 0; age < age_groups.size(); age++) {
                const auto& comp = age_groups[age];
                WritePredictionFullRow(output_full, iter, day, age, comp);
            }
        }
    }
}

void WritePredictionFullHeader(std::ostream& os)
{
    os << "iter, day, age_group, S, E, E_t, I_p, I_t,"
        " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D" << std::endl;
}

void WriteInferenceEndsHeader(std::ostream& os)
{
    os << "iter, age_group, S, E, E_t, I_p, I_t,"
        " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D" << std::endl;
}

void WriteInferenceParticlesHeader(std::ostream& os)
{
    os << "iter, nsse_cases, nsse_deaths, p_inf, "
        "p_hcw, c_hcw, d, q, p_s, rrd, lambda, weight" << std::endl;
}

void WritePredictionFullRow(std::ostream& os, int iter, int day, int age_group, const Compartments& comp)
{
    os << iter          << ", ";
    os << day           << ", ";
    os << age_group     << ", ";
    os << comp.S        << ", ";
    os << comp.E        << ", ";
    os << comp.E_t      << ", ";
    os << comp.I_p      << ", ";
    os << comp.I_t      << ", ";
    os << comp.I1       << ", ";
    os << comp.I2       << ", ";
    os << comp.I3       << ", ";
    os << comp.I4       << ", ";
    os << comp.I_s1     << ", ";
    os << comp.I_s2     << ", ";
    os << comp.I_s3     << ", ";
    os << comp.I_s4     << ", ";
    os << comp.H        << ", ";
    os << comp.R        << ", ";
    os << comp.D        << std::endl;
}

void WriteInferenceEndsRow(std::ostream& os, int iter, int age_group, const Compartments& comp)
{
    os << iter          << ", ";
    os << age_group     << ", ";
    os << comp.S        << ", ";
    os << comp.E        << ", ";
    os << comp.E_t      << ", ";
    os << comp.I_p      << ", ";
    os << comp.I_t      << ", ";
    os << comp.I1       << ", ";
    os << comp.I2       << ", ";
    os << comp.I3       << ", ";
    os << comp.I4       << ", ";
    os << comp.I_s1     << ", ";
    os << comp.I_s2     << ", ";
    os << comp.I_s3     << ", ";
    os << comp.I_s4     << ", ";
    os << comp.H        << ", ";
    os << comp.R        << ", ";
    os << comp.D        << std::endl;
}

void WriteInferenceParticlesRow(std::ostream& os, int iter, const particle particle)
{
    os << iter                                                     << ", ";
    os << particle.nsse_cases                                      << ", ";
    os << particle.nsse_deaths                                     << ", ";
    os << particle.parameter_set[Model::ModelParameters::PINF]     << ", ";
    os << particle.parameter_set[Model::ModelParameters::PHCW]     << ", ";
    os << particle.parameter_set[Model::ModelParameters::CHCW]     << ", ";
    os << particle.parameter_set[Model::ModelParameters::D]        << ", ";
    os << particle.parameter_set[Model::ModelParameters::Q]        << ", ";
    os << particle.parameter_set[Model::ModelParameters::PS]       << ", ";
    os << particle.parameter_set[Model::ModelParameters::RRD]      << ", ";
    os << particle.parameter_set[Model::ModelParameters::LAMBDA]   << ", ";
    os << particle.weight                                          << std::endl;
}

void WriteSimuHeader(std::ostream& os)
{
    os << "iter, day, " << "inc_case, " << "inc_death_hospital, " << "inc_death" << std::endl;
}

void WriteSimuRow(std::ostream& os, int iter, int day, int inc_case, int inc_death_hospital,
    int inc_death) 
{
    os << iter  << ", ";
    os << day   << ", ";
    os << inc_case  << ", ";
    os << inc_death_hospital << ", ";
    os << inc_death << std::endl;
}

void LogFixedParameters(const params& paramlist, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Fixed parameter values]:\n";
    OutputFixedParameters(log, paramlist);
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
    (*log) << "    n_iterations: "       << config.n_iterations << std::endl;
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
    OutputFixedParameters(log, config.fixedParameters);
}

void LogGitVersionInfo(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Git Versioning]"    << std::endl;
    (*log) << "    Commit URL: "    << GitMetadata::URL() << std::endl;
    (*log) << "    Commit SHA: "    << GitMetadata::CommitSHA1() << std::endl;
    (*log) << "    Commit Date: "   << GitMetadata::CommitDate() << std::endl;
    (*log) << "    Tag: "           << GitMetadata::Tag() << std::endl;
    (*log) << "    Uncommitted changes: " << (GitMetadata::AnyUncommittedChanges() ? "Yes" : "No") << std::endl;
}

void LogDependencyVersionInfo(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Dependency Versioning]"    << std::endl;
    (*log) << "    Compiler id: "    << DependencyVersions::CompilerId() << std::endl;
    (*log) << "    Compiler version: "    << DependencyVersions::CompilerVersion() << std::endl;
    (*log) << "    CMake version: "    << DependencyVersions::CMakeVersion() << std::endl;
    (*log) << "    GSL version: " << DependencyVersions::GSLVersion() << std::endl;
}

void OutputFixedParameters(Utilities::logging_stream::Sptr& log, const params& paramlist)
{
    (*log) << "    latent period (theta_l): " << paramlist.T_lat << std::endl;
    (*log) << "    pre-adult probability of symptoms devt (p_s[0]): " << paramlist.juvp_s << std::endl;
    (*log) << "    pre-clinical period (theta_i): " << paramlist.T_inf <<std::endl;
    (*log) << "    asymptomatic period (theta_r): " << paramlist.T_rec <<std::endl;
    (*log) << "    symptomatic period (theta_s): " << paramlist.T_sym <<std::endl;
    (*log) << "    hospitalisation stay (theta_h): " << paramlist.T_hos <<std::endl;
    (*log) << "    bed capacity at hospital (K): " << paramlist.K <<std::endl;
    (*log) << "    relative infectiousness of asymptomatic (u): " << paramlist.inf_asym <<std::endl;
}

} // namespace IO
} // namespace EERAModel
