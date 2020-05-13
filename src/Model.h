#pragma once

#include "ModelTypes.h"
#include <vector>
#include <random>
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
 */
void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::Observations observations,
		 gsl_rng* r,
		 std::mt19937& gen,
		 const std::string& outDirPath);

void model_select(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, int Npop, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths);

void select_obs(int& Npop, int& t_index, int& duration, int& day_intro, int& day_shut, 
	std::vector<int>& obsHosp_tmp, std::vector<int>& obsDeaths_tmp, 
	std::vector<std::vector<int> > data_tmp, std::vector<std::vector<int> > death_tmp, int herd_id, 
	int time_back);

} // namespace Model
} // namespace EERAModel