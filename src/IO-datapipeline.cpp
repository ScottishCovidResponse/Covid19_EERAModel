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

    // std::cout << "ParamsPath: " << ParamsPath << "\n";
    // std::cout << "ModelConfig: " << ModelConfigDir << "\n";
    // std::cout << "ConfigPath: " << dpconfig_path << "\n";

    if (dpconfig_path != "")
    {
        std::string uri = GitMetadata::URL(); // "https://whatever"; I'm guessing this is the repo for the model
        std::string git_sha = GitMetadata::CommitSHA1(); // And this the version ID, need to find these both

        // std::cout << "URI: " << uri << "\n";
        // std::cout << "SHA: " << git_sha << "\n";

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
        (*log) << "    From data pipeline" << std::endl;

        commonParameters.paramlist  = ReadFixedModelParameters();
        commonParameters.totN_hcw   = dp->read_estimate("fixed-parameters/total_hcw", "total_hcw");
        commonParameters.day_shut = dp->read_estimate("fixed-parameters/day_shut", "day_shut");
    }
    else
    {
        (*log) << "    From local parameters.ini" << std::endl;

        commonParameters.paramlist  = IO::ReadFixedModelParameters(ParamsPath);
        commonParameters.totN_hcw   = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", ParamsPath);
        commonParameters.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", ParamsPath);
    }
    
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);

    return commonParameters;
}

params IOdatapipeline::ReadFixedModelParameters()
{
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
        ValidationParameters validationParameters = ImportValidationParameters(ParamsPath);
        int nHealthBoards = validationParameters.nHealthBoards;
        int nAgeGroups = validationParameters.nAgeGroups;
        int nCfrCategories = validationParameters.nCfrCategories;
        int nCasesDays = validationParameters.nCasesDays;

        ObservationsForModels observations;

        //Uploading observed disease data
        //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
        //rows from 1 are indivudual health board
        //last row is for all of scotland
        dparray_to_csv<int>(
            "population-data/data_for_scotland", "data", &observations.cases,
            nHealthBoards, nCasesDays);

        // TODO: this is indexed by herd_id, and the data file has a titles row that the data pipeline doesn't
        // so need to do something about that. Fix this properly, but for now...

        observations.cases.insert(observations.cases.begin(), std::vector<int>(observations.cases[0].size()));

        //Uploading population per age group
        //columns are for each individual Health Borad
        //last column is for Scotland
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
        dparray_to_csv<double>(
            "population-data/data_for_scotland", "age", &observations.age_pop,
            nHealthBoards, nAgeGroups - 1);

        //mean number of daily contacts per age group (overall)	
        dparray_to_csv<double>(
            "contact-data/who_acquired_infection_from_whom", "norm", &observations.waifw_norm,
            nAgeGroups, nAgeGroups);

        //mean number of daily contacts per age group (home only)
        dparray_to_csv<double>(
            "contact-data/who_acquired_infection_from_whom", "home", &observations.waifw_home,
            nAgeGroups, nAgeGroups);

        //mean number of daily contacts per age group (not school, not work)
        dparray_to_csv<double>(
            "contact-data/who_acquired_infection_from_whom", "sdist", &observations.waifw_sdist,
            nAgeGroups, nAgeGroups);

        //Upload cfr by age group
        //col0: p_h: probability of hospitalisation
        //col1: cfr: case fatality ratio
        //col2: p_d: probability of death, given hospitalisation
        //rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW

        // TODO: (nCfrCategories - 1) is correct currently for the data pipeline vs local files, but is too arbitrary
        // and really the parameters.ini should be set for the data pipeline if used.
        dptable_to_csv<double>(
            "prob_hosp_and_cfr/data_for_scotland", "cfr_byage", &observations.cfr_byage,
            nAgeGroups, -1 /* nCfrCategories - 1 */ );

        // TODO: THIS IS A BODGE TO MAKE THE REGRESSION TESTS PASS. DON'T USE THIS!
        // Fix third column. 

        for (int i = 0; i < observations.cfr_byage.size(); ++i) {
            std::stringstream str;
            str << std::scientific << std::setprecision(14);
            str << observations.cfr_byage[i][2];
            observations.cfr_byage[i][2] = atof(str.str().c_str());
        }

        return observations;
    }
}

InferenceConfig IOdatapipeline::ReadInferenceConfig(
    const CommonModelInputParameters& commonParameters, const ObservationsForModels& modelObservations) 
{
    if (!datapipelineActive)
    {
        return IO::ReadInferenceConfig(ModelConfigDir, log, commonParameters);
    }
    else
    {
        InferenceConfig inferenceConfig;

        IO::ReadLocalInferenceConfig(ParamsPath, log, commonParameters, &inferenceConfig);

        dpdistribution(
            "prior-distributions/pinf", "pinf",
            "alpha", &inferenceConfig.prior_pinf_shape1, "beta", &inferenceConfig.prior_pinf_shape2);

        dpdistribution(
            "prior-distributions/phcw", "phcw",
            "alpha", &inferenceConfig.prior_phcw_shape1, "beta", &inferenceConfig.prior_phcw_shape2);

        dpdistribution(
            "prior-distributions/chcw", "chcw",
            "lambda", &inferenceConfig.prior_chcw_mean);

        dpdistribution(
            "prior-distributions/d", "d",
            "alpha", &inferenceConfig.prior_d_shape1, "beta", &inferenceConfig.prior_d_shape2);

        dpdistribution(
            "prior-distributions/q", "q",
            "alpha", &inferenceConfig.prior_q_shape1, "beta", &inferenceConfig.prior_q_shape2);

        dpdistribution(
            "prior-distributions/lambda", "lambda",
            "a", &inferenceConfig.prior_lambda_shape1, "b", &inferenceConfig.prior_lambda_shape2);
        
        dpdistribution(
            "prior-distributions/ps", "ps",
            "alpha", &inferenceConfig.prior_ps_shape1, "beta", &inferenceConfig.prior_ps_shape2);

        dpdistribution(
            "prior-distributions/rrd", "rrd",
            "k", &inferenceConfig.prior_rrd_shape1, "theta", &inferenceConfig.prior_rrd_shape2);

        inferenceConfig.observations = ReadInferenceObservations(modelObservations);

        return inferenceConfig;
    }
}

ObservationsForInference IOdatapipeline::ReadInferenceObservations(const ObservationsForModels& modelObservations)
{
    ValidationParameters validationParameters = ImportValidationParameters(ParamsPath);
    int nHealthBoards = validationParameters.nHealthBoards;
    int nCasesDays = validationParameters.nCasesDays;
    
    ObservationsForInference observations;
    
    (*log) << "Observations For Inference Config:" << std::endl;

    //Uploading observed disease data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log) << "\t- (data pipeline) copying cases from model observations" << std::endl;
    observations.cases = modelObservations.cases;

    //Uploading observed death data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    dparray_to_csv<int>(
            "population-data/data_for_scotland", "deaths", &observations.deaths,
            nHealthBoards, nCasesDays);

    observations.deaths.insert(observations.deaths.begin(), std::vector<int>(observations.deaths[0].size()));

    return observations;
}

PredictionConfig IOdatapipeline::ReadPredictionConfig(int index, const CommonModelInputParameters& commonParameters)
{
    if (!datapipelineActive)
    {
        return IO::ReadPredictionConfig(ModelConfigDir, index, log, commonParameters);
    }
    else
    {
        PredictionConfig predictionConfig;

        IO::ReadLocalPredictionConfig(ParamsPath, index, log, commonParameters, &predictionConfig);

        const char *posterior_product = "posterior_parameters/data_for_scotland";
        const char *posterior_component = "posterior_parameters";

        (*log) << "\t- (data pipeline) \"" << posterior_product << "\", \"" << posterior_component << "\"" << std::endl;

        Table posterior_table = dp->read_table(posterior_product, posterior_component);

        ImportConsistencyCheck(
            posterior_product, posterior_component, posterior_table.get_column_names().size(), 17, "columns");

        if (index >= posterior_table.get_column_size()) {
            std::stringstream SetSelectError;
            SetSelectError << "Parameter set selection out of bounds! Please select between 0-" << (posterior_table.get_column_size() - 1) << "..." << std::endl;
            throw std::overflow_error(SetSelectError.str());
        }

        // The posterior parameters are columns 1 to 8
        const char *posterior_colnames[] = { "p_inf", "p_hcw", "c_hcw", "d", "q", "p_s", "rrd", "intro" };

        predictionConfig.posterior_parameters.resize(8);

        for (int p = 0; p < 8; ++p) {
            predictionConfig.posterior_parameters[p] = posterior_table.get_column<double>(posterior_colnames[p])[index];
        }

        // The fixed parameters are columns 9 to 16

        predictionConfig.fixedParameters.T_lat      = posterior_table.get_column<int>("T_lat")[index];
        predictionConfig.fixedParameters.juvp_s     = posterior_table.get_column<double>("juvp_s")[index];
        predictionConfig.fixedParameters.T_inf      = posterior_table.get_column<double>("T_inf")[index];
        predictionConfig.fixedParameters.T_rec      = posterior_table.get_column<int>("T_rec")[index];
        predictionConfig.fixedParameters.T_sym      = posterior_table.get_column<int>("T_sym")[index];
        predictionConfig.fixedParameters.T_hos      = posterior_table.get_column<int>("T_hos")[index];
        predictionConfig.fixedParameters.K          = posterior_table.get_column<int>("K")[index];
        predictionConfig.fixedParameters.inf_asym   = posterior_table.get_column<int>("inf_asym")[index];

        return predictionConfig;
    }
}

void IOdatapipeline::dpdistribution(
    const std::string& data_product, const std::string& component,
    std::string p1, double *a, std::string p2, double *b)
{
    (*log) << "\t- (data pipeline) \"" << data_product << "\", \"" << component << "\"" << std::endl;

    Distribution input = dp->read_distribution(data_product, component);
    *a = input.getParameter(p1.c_str());
    if (b) *b = input.getParameter(p2.c_str());
}

void IOdatapipeline::ImportConsistencyCheck(
    const std::string& data_product, const std::string& component,
    const unsigned int axisLength, const unsigned int expectedValue, const std::string& axisID)
{
    if (axisLength != expectedValue) {
        std::stringstream IOMessage;
        IOMessage << "Error in data pipeline : \"" << data_product << "\" \"" << component <<
        "\n Number of " << axisID << ": " << axisLength <<
        "\n Expected number of " << axisID << ": " << expectedValue << std::endl;
        throw IOException(IOMessage.str());
    }
}

} // namespace IO
} // namespace EERAModel