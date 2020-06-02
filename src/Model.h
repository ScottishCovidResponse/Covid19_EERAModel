#pragma once

#include "ModelTypes.h"
#include <vector>
#include <random>
#include "Utilities.h"
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace Model {

/**
 * @brief Structure to hold per-age group data vectors
 * 
 * Includes probabilities and mean number of contacts per age group information
 */
struct AgeGroupData
{
	
	std::vector<std::vector<double>> waifw_norm;	/*!< mean number of daily contacts per age group (overall). */	
	std::vector<std::vector<double>> waifw_home;	/*!< mean number of daily contacts per age group (home only). */
	std::vector<std::vector<double>> waifw_sdist;	/*!< mean number of daily contacts per age group (not school, not work). */
	std::vector<std::vector<double>> cfr_byage;		/*!< Case fatality ratio by age. */
	std::vector<double> pf_byage;					/*!< Frailty Probability by age. */
};

/**
 * @brief Structure to hold status objects
 * 
 * Contains vectors containing information on the status of the simulation itself, 
 * population and deaths for each simulation step.
 */
struct Status
{
	std::vector<int> simulation;		/*!< Status of the simulation. */
	std::vector<int> deaths;			/*!< Number of deaths. */
	std::vector<int> hospital_deaths;	/*!< Number of deaths in hospitals. */
	std::vector<std::vector<int>> ends;	/*!< Population per age group on last day. */
};

/**
 * @brief Calculate number of simulation steps
 * 
 * Calculates the time steps for the simulation based on the duration time
 * and the step duration.
 * 
 * @param duration Duration of the simulation
 * @param tau Duration of a single time step
 * 
 * @return Number of simulation steps
 */
int get_n_simulation_steps(const int& duration, const double& tau)
{
	return static_cast<int>(ceil(duration/tau));
}

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

void model_select(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths);

/**
 * @brief Construct the population seed
 * 
 * Builds a vector of integers based on the input.
 * 
 * @param age_nums Vector containing the number of people in each age group
 * 
 * @return Vector containing populations
 */
std::vector<double> build_population_seed(const std::vector<int>& age_nums);

} // namespace Model
} // namespace EERAModel
