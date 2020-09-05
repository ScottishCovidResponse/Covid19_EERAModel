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

IOdatapipeline::IOdatapipeline(
    string params_path, string model_config, string outdir_path,
    Utilities::logging_stream::Sptr log_stream, string dpconfig_path)
{
    paramsPath_ = params_path;
    modelConfigDir_ = model_config; // Shouldn't need this if data pipeline active
    outDirPath_ = outdir_path; // nor this
    log_ = log_stream;

    datapipelineActive_ = false;

    timeStamp_ = Utilities::timeString();

    // std::cout << "paramsPath_: " << paramsPath_ << "\n";
    // std::cout << "modelConfigDir_: " << modelConfigDir_ << "\n";
    // std::cout << "outDirPath_: " << outDirPath_ << "\n";
    // std::cout << "DPConfigPath: " << dpconfig_path << "\n";

    if (dpconfig_path != "")
    {
        std::string uri = GitMetadata::URL(); // "https://whatever"; I'm guessing this is the repo for the model
        std::string git_sha = GitMetadata::CommitSHA1(); // And this the version ID, need to find these both

        // std::cout << "URI: " << uri << "\n";
        // std::cout << "SHA: " << git_sha << "\n";

        // If something goes wrong in the opening of the data pipeline an exception will get thrown
        dp_.reset(new DataPipeline(dpconfig_path, uri.c_str(), git_sha.c_str()));
        datapipelineActive_ = true;
    }
}

CommonModelInputParameters IOdatapipeline::ReadCommonParameters()
{
    CommonModelInputParameters commonParameters;

    if (datapipelineActive_)
    {
        (*log_) << "    From data pipeline" << std::endl;

        commonParameters.paramlist = ReadFixedModelParameters();
        commonParameters.totN_hcw = dp_->read_estimate("fixed-parameters/total_hcw", "total_hcw");
        commonParameters.day_shut = dp_->read_estimate("fixed-parameters/day_shut", "day_shut");
    }
    else
    {
        (*log_) << "    From local parameters.ini" << std::endl;

        commonParameters.paramlist = IO::ReadFixedModelParameters(paramsPath_);
        commonParameters.totN_hcw = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", paramsPath_);
        commonParameters.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", paramsPath_);
    }
    
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", paramsPath_);

    return commonParameters;
}

params IOdatapipeline::ReadFixedModelParameters()
{
    params paramlist;

    paramlist.T_lat = dp_->read_estimate("fixed-parameters/T_lat", "T_lat");
    paramlist.juvp_s = dp_->read_estimate("fixed-parameters/juvp_s", "juvp_s");
    paramlist.T_inf = dp_->read_estimate("fixed-parameters/T_inf", "T_inf");
    paramlist.T_rec = dp_->read_estimate("fixed-parameters/T_rec", "T_rec");
    paramlist.T_sym = dp_->read_estimate("fixed-parameters/T_sym", "T_sym");
    paramlist.T_hos = dp_->read_estimate("fixed-parameters/T_hos", "T_hos");
    paramlist.K = dp_->read_estimate("fixed-parameters/K", "K");
    paramlist.inf_asym = dp_->read_estimate("fixed-parameters/inf_asym", "inf_asym");

    return paramlist;
}

ObservationsForModels IOdatapipeline::ReadModelObservations()
{
    if (!datapipelineActive_)
    {
        return IO::ReadModelObservations(modelConfigDir_, log_);
    }
    else
    {
        ValidationParameters validationParameters = ImportValidationParameters(paramsPath_);
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
    if (!datapipelineActive_)
    {
        return IO::ReadInferenceConfig(modelConfigDir_, log_, commonParameters);
    }
    else
    {
        InferenceConfig inferenceConfig;

        IO::ReadLocalInferenceConfig(paramsPath_, log_, commonParameters, &inferenceConfig);

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
    ValidationParameters validationParameters = ImportValidationParameters(paramsPath_);
    int nHealthBoards = validationParameters.nHealthBoards;
    int nCasesDays = validationParameters.nCasesDays;
    
    ObservationsForInference observations;
    
    (*log_) << "Observations For Inference Config:" << std::endl;

    //Uploading observed disease data
    //Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
    //rows from 1 are indivudual health board
    //last row is for all of scotland
    (*log_) << "\t- (data pipeline) copying cases from model observations" << std::endl;
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
    if (!datapipelineActive_)
    {
        return IO::ReadPredictionConfig(modelConfigDir_, index, log_, commonParameters);
    }
    else
    {
        PredictionConfig predictionConfig;

        IO::ReadLocalPredictionConfig(paramsPath_, index, log_, commonParameters, &predictionConfig);

        const char *posterior_product = "posterior_parameters/data_for_scotland";
        const char *posterior_component = "posterior_parameters";

        (*log_) << "\t- (data pipeline) \"" << posterior_product << "\", \"" << posterior_component << "\"" << std::endl;

        Table posterior_table = dp_->read_table(posterior_product, posterior_component);

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

void IOdatapipeline::WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, 
    const std::vector<particle>& particleList, const std::string &modelType)
{
    if (!datapipelineActive_)
    {
        IO::WriteOutputsToFiles(smc, herd_id, Nparticle, nPar, particleList, outDirPath_, log_);
    }
    else
    {
        std::string product = "outputs/" + modelType + "/inference/" + timeStamp_;
        std::string component_prefix = "steps/" + std::to_string(smc);

        std::cout << "Writing product: " << product << "\n";

        // write_table("original/inference/2020-08-29_10_11_00", "log_", log_table);
        // for step in 0:nsteps
        //   write_table("eera_outputs/original/inference/2020-08-29_10_11_00", "steps/step/ends", ends_step_n);
        //   write_table("eera_outputs/original/inference", "steps/step/simu", simu_step_n);
        //   write_table("eera_outputs/original/inference", "steps/step/particles", particles_step_n);
        // ends

        // Particles

        {
            std::string component_particles = component_prefix + "/particles";
            std::cout << "    " << component_particles << "\n";

            Table output_part;
            WriteInferenceParticlesTable(Nparticle, particleList, &output_part);
            dp_->write_table(product, component_particles, output_part);
        }

        // Simu

        {
            std::string component_simu = component_prefix + "/simu";
            std::cout << "    " << component_simu << "\n";

            Table output_simu;
            WriteInferenceSimuTable(Nparticle, particleList, &output_simu);
            dp_->write_table(product, component_simu, output_simu);
        }

        // Ends

        {
            std::string component_ends = component_prefix + "/ends";
            std::cout << "    " << component_ends << "\n";

            Table output_ends;
            WriteInferenceEndsTable(Nparticle, particleList, &output_ends);
            dp_->write_table(product, component_ends, output_ends);
        }
    }
}

void IOdatapipeline::WriteInferenceParticlesTable(
    int Nparticle, const std::vector<particle>& particleList, Table *table)
{
    std::vector<int> int_values(Nparticle);
    std::vector<double> double_values(Nparticle);

    auto intcolumn =
        [&] (const char *colname, std::function<int(const particle& particle)> element) {
            for (int kk = 0; kk < Nparticle; ++kk) {
                int_values[kk] = element(particleList[kk]);
            }

            table->add_column(colname, int_values);
        };

    auto doublecolumn =
        [&] (const char *colname, std::function<double(const particle& particle)> element) {
            for (int kk = 0; kk < Nparticle; ++kk) {
                double_values[kk] = element(particleList[kk]);
            }

            table->add_column(colname, double_values);
        };

    intcolumn("iter", [=] (const particle& particle) -> int {
        return particle.iter; });

    doublecolumn("nsse_cases", [=] (const particle& particle) -> double {
        return particle.nsse_cases; });

    doublecolumn("nsse_deaths", [=] (const particle& particle) -> double {
        return particle.nsse_deaths; });

    doublecolumn("p_inf", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PINF]; });

    doublecolumn("p_hcw", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PHCW]; });

    doublecolumn("c_hcw", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::CHCW]; });

    doublecolumn("d", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::D]; });

    doublecolumn("q", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::Q]; });

    doublecolumn("p_s", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PS]; });

    doublecolumn("rrd", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::RRD]; });

    doublecolumn("lambda", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::LAMBDA]; });

    doublecolumn("weight", [=] (const particle& particle) -> double {
        return particle.weight; });
}

void IOdatapipeline::WriteInferenceSimuTable(
    int Nparticle, const std::vector<particle>& particleList, Table *table)
{
    size_t vectorsize = 0;

    for (const auto& particle : particleList) {
        vectorsize += particle.simu_outs.size();
    }

    std::vector<int> int_values(vectorsize);

    auto intcolumn =
        [&] (const char *colname, std::function<int(const particle& particle, int var)> element) {
            size_t dst = 0;

            for (int kk = 0; kk < Nparticle; ++kk) {
                size_t nvar = particleList[kk].simu_outs.size();
                for (unsigned int var = 0; var < nvar; ++var) {
                    int_values[dst] = element(particleList[kk], var);
                    ++dst;
                }
            }

            table->add_column(colname, int_values);
        };

    intcolumn("iter", [=] (const particle& particle, int var) -> int {
        return particle.iter; });

    intcolumn("day", [=] (const particle& particle, int var) -> int {
        return var; });

    intcolumn("inc_case", [=] (const particle& particle, int var) -> int {
        return particle.simu_outs[var]; });

    intcolumn("inc_death_hospital", [=] (const particle& particle, int var) -> int {
        return particle.hospital_death_outs[var]; });

    intcolumn("inc_death", [=] (const particle& particle, int var) -> int {
        return particle.death_outs[var]; });
}

void IOdatapipeline::WriteInferenceEndsTable(
    int Nparticle, const std::vector<particle>& particleList, Table *table)
{
    size_t vectorsize = 0;

    for (const auto& particle : particleList) {
        vectorsize += particle.end_comps.size();
    }

    std::vector<int> int_values(vectorsize);

    auto intcolumn =
        [&] (const char *colname, std::function<int(const particle& particle, int age)> element) {
            size_t dst = 0;

            for (int kk = 0; kk < Nparticle; ++kk) {
                size_t nage = particleList[kk].end_comps.size();
                for (unsigned int age = 0; age < nage; ++age) {
                    int_values[dst] = element(particleList[kk], age);
                    ++dst;
                }
            }

            table->add_column(colname, int_values);
        };

    intcolumn("iter", [=] (const particle& particle, int age) -> int {
        return particle.iter; });

    intcolumn("age_group", [=] (const particle& particle, int age) -> int {
        return age; });

    intcolumn("S", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].S; });

    intcolumn("E", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].E; });

    intcolumn("E_t", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].E_t; });

    intcolumn("I_p", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_p; });

    intcolumn("I_t", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_t; });

    intcolumn("I1", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I1; });

    intcolumn("I2", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I2; });

    intcolumn("I3", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I3; });

    intcolumn("I4", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I4; });

    intcolumn("I_s1", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_s1; });

    intcolumn("I_s2", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_s2; });

    intcolumn("I_s3", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_s3; });

    intcolumn("I_s4", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].I_s4; });

    intcolumn("H", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].H; });

    intcolumn("R", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].R; });

    intcolumn("D", [=] (const particle& particle, int age) -> int {
        return particle.end_comps[age].D; });
}


void IOdatapipeline::WritePredictionsToFiles(std::vector<Status> statuses, const std::string &modelType)
{
    if (!datapipelineActive_)
    {
        IO::WritePredictionsToFiles(statuses, outDirPath_, log_);
    }
    else
    {
        // write_table("original/prediction/2020-08-29_10_11_00", "log_", log_table);
        // write_table("eera_outputs/original/prediction", "full", full);
        // write_table("eera_outputs/original/prediction", "simu", simu);
        
        std::string product = "outputs/" + modelType + "/prediction/" + timeStamp_;
        std::cout << "Writing product: " << product << "\n";

        // Simu

        {
            std::string component_simu = "simu";
            std::cout << "    " << component_simu << "\n";

            Table output_simu;
            WritePredictionSimuTable(statuses, &output_simu);
            dp_->write_table(product, component_simu, output_simu);
        }
    
        // Full

        {
            std::string component_full = "full";
            std::cout << "    " << component_full << "\n";

            Table output_full;
            WritePredictionFullTable(statuses, &output_full);
            dp_->write_table(product, component_full, output_full);
        }
    }
}

void IOdatapipeline::WritePredictionSimuTable(std::vector<Status> statuses, Table *table)
{
    size_t vectorsize = 0;

    for (const auto& status : statuses) {
        vectorsize += status.simulation.size();
    }

    std::vector<int> int_values(vectorsize);

    auto intcolumn =
        [&] (const char *colname, std::function<int(const Status& status, int iter, int day)> element) {
            size_t dst = 0, nstatus = statuses.size();

            for (int iter = 0; iter < nstatus; ++iter) {
                size_t nday = statuses[iter].simulation.size();
                for (unsigned int day = 0; day < nday; ++day) {
                    int_values[dst] = element(statuses[iter], iter, day);
                    ++dst;
                }
            }

            table->add_column(colname, int_values);
        };

    intcolumn("iter", [=] (const Status& status, int iter, int day) -> int {
        return iter; });

    intcolumn("day", [=] (const Status& status, int iter, int day) -> int {
        return day; });

    intcolumn("inc_case", [=] (const Status& status, int iter, int day) -> int {
        return status.simulation[day]; });

    intcolumn("inc_death_hospital", [=] (const Status& status, int iter, int day) -> int {
        return status.hospital_deaths[day]; });

    intcolumn("inc_death", [=] (const Status& status, int iter, int day) -> int {
        return status.deaths[day]; });
}

void IOdatapipeline::WritePredictionFullTable(std::vector<Status> statuses, Table *table)
{
    size_t vectorsize = 0;

    for (const auto& status : statuses) {
        for (const auto& age_group : status.pop_array) {
            vectorsize += age_group.size();
        }
    }

    std::vector<int> int_values(vectorsize);

    auto intcolumn =
        [&] (const char *colname, std::function<int(const Compartments& comp, int iter, int day, int age)> element) {
            size_t dst = 0, nstatus = statuses.size();

            for (int iter = 0; iter < nstatus; ++iter) {
                size_t nday = statuses[iter].pop_array.size();
                for (unsigned int day = 0; day < nday; ++day) {
                    size_t nage = statuses[iter].pop_array[day].size();
                    for (unsigned int age = 0; age < nage; ++age) {
                        int_values[dst] = element(statuses[iter].pop_array[day][age], iter, day, age);
                        ++dst;
                    }
                }
            }

            table->add_column(colname, int_values);
        };

    intcolumn("iter", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return iter; });

    intcolumn("day", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return day; });

    intcolumn("age_group", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return age; });

    intcolumn("S", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.S; });

    intcolumn("E", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.E; });

    intcolumn("E_t", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.E_t; });

    intcolumn("I_p", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_p; });

    intcolumn("I_t", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_t; });

    intcolumn("I1", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I1; });

    intcolumn("I2", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I2; });

    intcolumn("I3", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I3; });

    intcolumn("I4", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I4; });

    intcolumn("I_s1", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_s1; });

    intcolumn("I_s2", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_s2; });

    intcolumn("I_s3", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_s3; });

    intcolumn("I_s4", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.I_s4; });

    intcolumn("H", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.H; });

    intcolumn("R", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.R; });

    intcolumn("D", [=] (const Compartments& comp, int iter, int day, int age) -> int {
        return comp.D; });
}

// Some utilities ----------------------------------------

void IOdatapipeline::dpdistribution(
    const std::string& data_product, const std::string& component,
    std::string p1, double *a, std::string p2, double *b)
{
    (*log_) << "\t- (data pipeline) \"" << data_product << "\", \"" << component << "\"" << std::endl;

    Distribution input = dp_->read_distribution(data_product, component);
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