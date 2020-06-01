#pragma once

#include "ModelTypes.h"
#include <vector>
#include <random>
#include "Utilities.h"
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace Model {

/**
 * @brief Run the model and inference framework
 * 
 * Runs the model based on the given input parameters, observations and seeded random number generator.
 * Places the outputs in the indicated directory.
 * 
 * @param modelInputParameters Input parameters to the model run
 * @param observations Input observations to the model run
 * @param r Seeded random number generator
 * @param gen Seeded random number generator for importance sampling
 * @param outDirPath Path to the directory in which the output files should be placed
 * @param log Pointer to the logger
 */
void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::InputObservations observations,
		 gsl_rng* r,
		 std::mt19937& gen,
		 const std::string& outDirPath,
		 EERAModel::Utilities::logging_stream::Sptr log);

// void model_select(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
// 	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
// 	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
// 	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
// 	int duration, seed seedlist, int day_shut, gsl_rng * r, const std::vector<int>& obsHosp,
// 	const std::vector<int>& obsDeaths);

void model_select(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
	::EERAModel::InputParametersObservations &parameters_in, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths);

} // namespace Model
} // namespace EERAModel
