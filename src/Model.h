#pragma once

#include "ModelTypes.h"
#include <vector>
#include <map>
#include <random>
#include <sstream>
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

struct Compartments
{
	int S = 0;
	int E = 0;
	int E_t = 0;
	int I_p = 0;
	int I_t = 0;
	int I1 = 0;
	int I2 = 0;
	int I3 = 0;
	int I4 = 0;
	int I_s1 = 0;
	int I_s2 = 0;
	int I_s3 = 0;
	int I_s4 = 0;
	int H = 0;
	int R = 0;
	int D = 0;
};

/**
 * @brief Accumulate total across all compartments
 * 
 * Totals all counters within every compartment in a Compartments struct
 * 
 * @param comp Compartments object
 * 
 * @return Total across all compartments
 */
int accumulate_compartments(const Compartments& comp)
{
	int _total = 0;
	_total += comp.S + comp.E + comp.E_t + comp.I_p;
	_total += comp.I_t + comp.I1 + comp.I2 + comp.I3;
	_total += comp.I4 + comp.I_s1 + comp.I_s2 + comp.I_s3;
	_total += comp.I_s4 + comp.H + comp.R + comp.D;

	return _total;
}

/**
 * @brief Convert Vector of Compartments struct to a vector of integers
 * 
 * NOTE: This is a temporary function to allow compatibility
 * converts the Compartments struct to a vector of integers
 * 
 * @param cmps_vec Vector of compartments struct containing population per category
 * 
 * @return Vector of population counters
 */
std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec)
{
	std::vector<std::vector<int>> _temp;

	for(auto cmps : cmps_vec)
	{
		_temp.push_back({cmps.S, cmps.E, cmps.E_t, cmps.I_p,
						cmps.I_t, cmps.I1, cmps.I2, cmps.I3,
						cmps.I4, cmps.I_s1, cmps.I_s2, cmps.I_s3,
						cmps.I_s4, cmps.H, cmps.R, cmps.D});
	}

	return _temp;
}

/**
 * @brief Structure to hold status objects
 * 
 * Contains vectors containing information on the status of the simulation itself, 
 * population and deaths for each simulation step.
 */
struct Status
{
	std::vector<int> simulation;				/*!< Status of the simulation. */
	std::vector<int> deaths;					/*!< Number of deaths. */
	std::vector<int> hospital_deaths;			/*!< Number of deaths in hospitals. */
	std::vector<Compartments> ends;					/*!< Population per Category per age group on last day. */
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

/**
 * @brief Construct the population array
 * 
 * Sets up an array for the population at each timestep in each age and disease category	
 * also set up the age distribution of old ages as target for disease introduction.
 * 
 * @param age_nums Vector containing the number of people in each age group
 * 
 * @return Vector of vectors containing compartment populations
 */ 
std::vector<Compartments> build_population_array(const std::vector<int>& age_nums);

/**
 * @brief Randomly assign movement of individuals between compartments
 * 
 * Uses a Poisson distribution for a given rate to deduce number of individuals
 * passing from one compartment to another. Also compares current value prior to flow with that after to ensure
 * flow is positive.
 * 
 * @param pops_from_val Number of individuals in starting compartment within prior compartments struct
 * @param pops_new_from_val Number of individuals in starting compartment within the compartments struct currently being modified
 * @param rate The rate of flow 
 * 
 * @return Number of individuals moving/flowing from one compartment to the other
 */
int flow(gsl_rng*r, const int& pops_from_val, const int& pops_new_from_val, const double& rate);

} // namespace Model
} // namespace EERAModel
