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

        std::string component_ends = component_prefix + "/ends";
        std::string component_simu = component_prefix + "/simu";

        std::cout << "Writing product: " << product << "\n";
        std::cout << "    " << component_ends << "\n";
        std::cout << "    " << component_simu << "\n";

        // write_table("original/inference/2020-08-29_10_11_00", "log_", log_table);
        // for step in 0:nsteps
        //   write_table("eera_outputs/original/inference/2020-08-29_10_11_00", "steps/step/ends", ends_step_n);
        //   write_table("eera_outputs/original/inference", "steps/step/simu", simu_step_n);
        //   write_table("eera_outputs/original/inference", "steps/step/particles", particles_step_n);
        // ends

        // So, these tables can only have columns added to them, but the other write routine is trying
        // to output rows, so some flippage needed, ho-hum...

        std::string component_particles = component_prefix + "/particles";
        std::cout << "    " << component_particles << "\n";

        Table output_part;
        WriteInferenceParticlesTable(Nparticle, particleList, &output_part);
        dp_->write_table(product, component_particles, output_part);


    //     std::stringstream namefile, namefile_simu, namefile_ends;
    //     namefile << (outDirPath_ + "/output_abc-smc_particles_step") << smc << "_shb"<< herd_id << "_" << log_->getLoggerTime() << ".txt";
    //     namefile_simu << (outDirPath_ + "/output_abc-smc_simu_step") << smc << "_shb"<< herd_id << "_" << log_->getLoggerTime() << ".txt";
    //     namefile_ends << (outDirPath_ + "/output_abc-smc_ends_step") << smc << "_shb"<< herd_id << "_" << log_->getLoggerTime() << ".txt";		

    //     std::ofstream output_step (namefile.str().c_str());
    //     std::ofstream output_simu (namefile_simu.str().c_str());
    //     std::ofstream output_ends (namefile_ends.str().c_str());
        
    //     //add the column names for each output list of particles
    //     WriteInferenceParticlesHeader(output_step);

    //     //add the column names for each output list of chosen simulations
    //     WriteSimuHeader(output_simu);
        
    //     //add the column names for each output list of the compartment values of the last day of the chosen simulations
    //     WriteInferenceEndsHeader(output_ends);	

    //     // outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
    //     for (int kk = 0; kk < Nparticle; ++kk) {
    //         WriteInferenceParticlesRow(output_step, kk, particleList[kk]);
            
    //         for (unsigned int var = 0; var < particleList[kk].simu_outs.size(); ++var) {
    //             WriteSimuRow(output_simu, particleList[kk].iter, var , particleList[kk].simu_outs[var],
    //                         particleList[kk].hospital_death_outs[var], particleList[kk].death_outs[var]);
    //         }
            
    //         for (unsigned int age = 0; age < particleList[kk].end_comps.size(); ++age) {
    //             WriteInferenceEndsRow(output_ends, particleList[kk].iter, age, particleList[kk].end_comps[age]);
    //         }
    //     }
        
    //     output_step.close();
    //     output_simu.close();
    //     output_ends.close();
    }
}

void IOdatapipeline::WriteInferenceParticlesTable(
    int Nparticle, const std::vector<particle>& particleList, Table *table)
{
    std::vector<int> int_values(Nparticle);
    std::vector<double> double_values(Nparticle);

    for (int kk = 0; kk < Nparticle; ++kk) {
        int_values[kk] = kk;
    }

    table->add_column("iter", int_values);

    auto buildcolumn = [&] (const char *colname, std::function<double(const particle& particle)> element) {
        for (int kk = 0; kk < Nparticle; ++kk) {
            const auto& particle = particleList[kk];
            double_values[kk] = element(particle);
        }

        table->add_column(colname, double_values);
    };

    buildcolumn("nsse_cases", [=] (const particle& particle) -> double {
        return particle.nsse_cases; });

    buildcolumn("nsse_deaths", [=] (const particle& particle) -> double {
        return particle.nsse_deaths; });

    buildcolumn("p_inf", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PINF]; });

    buildcolumn("p_hcw", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PHCW]; });

    buildcolumn("c_hcw", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::CHCW]; });

    buildcolumn("d", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::D]; });

    buildcolumn("q", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::Q]; });

    buildcolumn("p_s", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::PS]; });

    buildcolumn("rrd", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::RRD]; });

    buildcolumn("lambda", [=] (const particle& particle) -> double {
        return particle.parameter_set[Model::ModelParameters::LAMBDA]; });

    buildcolumn("weight", [=] (const particle& particle) -> double {
        return particle.weight; });
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

    // std::stringstream namefile_simu, namefile_full;
    // namefile_simu << (outDirPath_ + "/output_prediction_simu") << "_" << log_->getLoggerTime() << ".txt";
    // namefile_full << (outDirPath_ + "/output_prediction_full") << "_" << log_->getLoggerTime() << ".txt";

    // std::ofstream output_simu(namefile_simu.str());
    // std::ofstream output_full(namefile_full.str());

    // WriteSimuHeader(output_simu);
    // for (unsigned int iter = 0; iter < statuses.size(); ++iter) {
    //     const Status& status = statuses[iter];
        
    //     for (unsigned int day = 0; day < status.simulation.size(); ++day) {
    //         WriteSimuRow(output_simu, iter, day, status.simulation[day],
    //             status.hospital_deaths[day], status.deaths[day]);
    //     }
    // }

    // WritePredictionFullHeader(output_full);
    // for (unsigned int iter = 0; iter < statuses.size(); ++iter) {
    //     const Status& status = statuses[iter];
    //     const auto& pop_array = status.pop_array;
        
    //     for (unsigned int day = 0; day < pop_array.size(); ++day) {
    //         const auto& age_groups = pop_array[day];
            
    //         for (unsigned int age = 0; age < age_groups.size(); age++) {
    //             const auto& comp = age_groups[age];
    //             WritePredictionFullRow(output_full, iter, day, age, comp);
    //         }
    //     }
    // }
    }
}

// void WritePredictionFullHeader(std::ostream& os)
// {
//     os << "iter, day, age_group, S, E, E_t, I_p, I_t,"
//         " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D" << std::endl;
// }

// void WriteInferenceEndsHeader(std::ostream& os)
// {
//     os << "iter, age_group, S, E, E_t, I_p, I_t,"
//         " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D" << std::endl;
// }

// void WriteInferenceParticlesHeader(std::ostream& os)
// {
//     os << "iter, nsse_cases, nsse_deaths, p_inf, "
//         "p_hcw, c_hcw, d, q, p_s, rrd, lambda, weight" << std::endl;
// }

// void WritePredictionFullRow(std::ostream& os, int iter, int day, int age_group, const Compartments& comp)
// {
//     os << iter          << ", ";
//     os << day           << ", ";
//     os << age_group     << ", ";
//     os << comp.S        << ", ";
//     os << comp.E        << ", ";
//     os << comp.E_t      << ", ";
//     os << comp.I_p      << ", ";
//     os << comp.I_t      << ", ";
//     os << comp.I1       << ", ";
//     os << comp.I2       << ", ";
//     os << comp.I3       << ", ";
//     os << comp.I4       << ", ";
//     os << comp.I_s1     << ", ";
//     os << comp.I_s2     << ", ";
//     os << comp.I_s3     << ", ";
//     os << comp.I_s4     << ", ";
//     os << comp.H        << ", ";
//     os << comp.R        << ", ";
//     os << comp.D        << std::endl;
// }

// void WriteInferenceEndsRow(std::ostream& os, int iter, int age_group, const Compartments& comp)
// {
//     os << iter          << ", ";
//     os << age_group     << ", ";
//     os << comp.S        << ", ";
//     os << comp.E        << ", ";
//     os << comp.E_t      << ", ";
//     os << comp.I_p      << ", ";
//     os << comp.I_t      << ", ";
//     os << comp.I1       << ", ";
//     os << comp.I2       << ", ";
//     os << comp.I3       << ", ";
//     os << comp.I4       << ", ";
//     os << comp.I_s1     << ", ";
//     os << comp.I_s2     << ", ";
//     os << comp.I_s3     << ", ";
//     os << comp.I_s4     << ", ";
//     os << comp.H        << ", ";
//     os << comp.R        << ", ";
//     os << comp.D        << std::endl;
// }

// void WriteInferenceParticlesRow(std::ostream& os, int iter, const particle particle)
// {
//     os << iter                                                     << ", ";
//     os << particle.nsse_cases                                      << ", ";
//     os << particle.nsse_deaths                                     << ", ";
//     os << particle.parameter_set[Model::ModelParameters::PINF]     << ", ";
//     os << particle.parameter_set[Model::ModelParameters::PHCW]     << ", ";
//     os << particle.parameter_set[Model::ModelParameters::CHCW]     << ", ";
//     os << particle.parameter_set[Model::ModelParameters::D]        << ", ";
//     os << particle.parameter_set[Model::ModelParameters::Q]        << ", ";
//     os << particle.parameter_set[Model::ModelParameters::PS]       << ", ";
//     os << particle.parameter_set[Model::ModelParameters::RRD]      << ", ";
//     os << particle.parameter_set[Model::ModelParameters::LAMBDA]   << ", ";
//     os << particle.weight                                          << std::endl;
// }

// void WriteSimuHeader(std::ostream& os)
// {
//     os << "iter, day, " << "inc_case, " << "inc_death_hospital, " << "inc_death" << std::endl;
// }

// void WriteSimuRow(std::ostream& os, int iter, int day, int inc_case, int inc_death_hospital,
//     int inc_death) 
// {
//     os << iter  << ", ";
//     os << day   << ", ";
//     os << inc_case  << ", ";
//     os << inc_death_hospital << ", ";
//     os << inc_death << std::endl;
// }

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