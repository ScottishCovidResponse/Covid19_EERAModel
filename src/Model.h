#pragma once

#include "ModelTypes.h"
#include <vector>
#include <random>
#include "Utilities.h"
#include "Random.h"

namespace EERAModel {
namespace Model {

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
 * @brief Run the model and inference framework
 * 
 * Runs the model based on the given input parameters, observations and seeded random number generator.
 * Places the outputs in the indicated directory.
 * 
 * @param modelInputParameters Input parameters to the model run
 * @param observations Input observations to the model run
 * @param rng Seeded random number generator
 * @param gen Seeded random number generator for importance sampling
 * @param outDirPath Path to the directory in which the output files should be placed
 * @param log Pointer to the logger
 */
void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::InputObservations observations,
		 Random::RNGInterface::Sptr rng,
		 const std::string& outDirPath,
		 EERAModel::Utilities::logging_stream::Sptr log);

void model_select(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, Random::RNGInterface::Sptr rng, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths, ModelStructureId structure);

/**
 * @brief Run the model with the given parameters and configurations
 * 
 * @param parameter_set: set of parameters that are being infered (i.e. particles)
 * @param fixed_parameters: set of fixed (known) parameters
 * @param per_age_data: age-structured known parameters (such as case fatality ratio (CFR) and probability of severe clinical outcomes )
 * @param seedlist: seeding method to initialise infection ("random": randomly allocate n infectious individuaals at the start, or "background": background transmission over a defined pre-lockdown period )
 * @param day_shut: day of the lock down
 * @param agenums: number of individuals in each age group in the study area
 * @param n_sim_steps: Number of steps to simulate
 * @param structure Model structure id
 * @param rng: Seeded random number generator
 * 
 * @return Status of model after run
 */
Status RunModel(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
				int n_sim_steps, ModelStructureId structure, Random::RNGInterface::Sptr rng);

/**
 * @brief Construct the population seed, based on a choice of model structure
 * 
 * The population seed is constructed according to the following criteria:
 *    - The final age category (HCW) are omitted from the seeding always
 *    - The first age category (< 20yo) are omiited in the original model structure; they are 
 *      retained in the Irish model structure
 * 
 * @param age_nums Vector containing the number of people in each age category
 * @param structure Model structure id
 * 
 * @return Seed population
 */
std::vector<double> BuildPopulationSeed(const std::vector<int>& age_nums, ModelStructureId structure);

/**
 * @brief Construct the population array, based on a choice of model structure
 * 
 * Sets up an array for the population at each timestep in each age and disease category	
 * also set up the age distribution of old ages as target for disease introduction.
 * 
 * @param rng Seeded random number generator
 * @param age_nums Vector containing the number of people in each age group
 * @param seedlist Seed object
 * @param structure Model structure id
 * 
 * @return Vector of vectors containing compartment populations
 */ 
std::vector<Compartments> BuildPopulationArray(Random::RNGInterface::Sptr rng,
    const std::vector<int>& age_nums, const seed& seedlist, ModelStructureId structure);

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
int Flow(Random::RNGInterface::Sptr rng, int pops_from_val, int pops_new_from_val, double rate);

/**
 * @brief Introduced diseased to the population, based on a choice of model structure
 * 
 * Compute the total number of susceptible and the number of susceptible per age class
 * 
 * @param rng Seeded random number generator
 * @param poparray Population array to be manipulated
 * @param seedarray Population seed array to be manipulated
 * @param bkg_lambda Lambda for generating number of diseased individuals
 * @param structure Model structure id
 */
void GenerateDiseasedPopulation(Random::RNGInterface::Sptr rng,
    std::vector<Compartments>& poparray, std::vector<double>& seedarray,
    const double& bkg_lambda, ModelStructureId structure);

/**
 * @brief Generate vector of lambda values for the age groups
 * 
 * Creates lambda values based on compartment occupancy for each age group
 * 
 * @param inf_hosp Number of hospitalised infected
 * @param parameter_set Set of model parameters
 * @param u_val
 * @param age_data Age group data set
 * @param pops Population array containing compartments for each age group
 * @param shut State of lockdown
 */
std::vector<double> GenerateForcesOfInfection(int& inf_hosp, const std::vector<double>& parameter_set, double u_val, 
			const AgeGroupData& age_data, const std::vector<Compartments>& pops, bool shut);

/**
 * @brief Generate an infection spread as per the structure of the original EERA model
 * 
 * Creates an infection spread state and counters number of people in different states
 * 
 * @param rng Seeded random number generator
 * @param pop Particular population age group
 * @param n_hospitalised Number of people in hospital prior to new spread
 * @param fixed_parameters Fixed model parameters
 * @param parameter_set Variable model parameters
 * @param cfr_tab Case Fatality Ratio table
 * @param pf_val Frailty Probability
 * @param lambda Rate of spread
 */
InfectionState GenerateInfectionSpreadOriginal(Random::RNGInterface::Sptr rng, Compartments& pop,
    const int& n_hospitalised, ::EERAModel::params fixed_parameters, 
    std::vector<double> parameter_set, std::vector<double> cfr_tab,
    double pf_val, double lambda);

/**
 * @brief Generate an infection spread as per the structure of the Irish model
 * 
 * Creates an infection spread state and counters number of people in different states
 * 
 * @param rng Seeded random number generator
 * @param pop Particular population age group
 * @param n_hospitalised Number of people in hospital prior to new spread
 * @param fixed_parameters Fixed model parameters
 * @param parameter_set Variable model parameters
 * @param cfr_tab Case Fatality Ratio table
 * @param pf_val Frailty Probability
 * @param lambda Rate of spread
 */
InfectionState GenerateInfectionSpreadIrish(Random::RNGInterface::Sptr rng, Compartments& pop,
    const int& n_hospitalised, ::EERAModel::params fixed_parameters, 
    std::vector<double> parameter_set, std::vector<double> cfr_tab,
    double pf_val, double lambda);

} // namespace Model
} // namespace EERAModel
