#pragma once

#include "ModelTypes.h"
#include <vector>
#include <random>
#include "Utilities.h"
#include "Random.h"

namespace EERAModel {
namespace Model {

/**
 * @brief Structure to hold per-age group data vectors
 * 
 * Includes probabilities and mean number of contacts per age group information
 */
struct AgeGroupData
{
	
	std::vector<std::vector<double>> waifw_norm;	/*!< mean number of daily contacts between age groups (overall). */	
	std::vector<std::vector<double>> waifw_home;	/*!< mean number of daily contacts between age groups (home only). */
	std::vector<std::vector<double>> waifw_sdist;	/*!< mean number of daily contacts between age groups (not school, not work). */
	std::vector<std::vector<double>> cfr_byage;		/*!< Case fatality ratio by age. */
	std::vector<double> pf_byage;					/*!< Frailty Probability by age. */
};

/**
 * @brief Structure containing counters for different population categories
 * 
 * Integers to count the number of people within the different compartments within the model 
 */
struct Compartments
{
	int S = 0;    	/*!< Number of susceptible individuals (not infected). */	
	int E = 0;		/*!< Number of infected individuals but not yet infectious (exposed). */	
	int E_t = 0;	/*!< Number of exposed individuals and tested positive. */	
	int I_p = 0;	/*!< Number of infected and infectious symptomatic individuals but at pre-clinical stage (show yet no symptoms). */	
	int I_t = 0;	/*!< Number of tested positive individuals that infectious. */	
	int I1 = 0;		/*!< Number of infected and infectious asymptomatic individuals: first stage. */	
	int I2 = 0;		/*!< Number of infected and infectious asymptomatic individuals: second stage.  */	
	int I3 = 0;		/*!< Number of infected and infectious asymptomatic individuals: third stage.  */	
	int I4 = 0;		/*!< Number of infected and infectious asymptomatic individuals: last stage.  */	
	int I_s1 = 0;	/*!< Number of infected and infectious symptomatic individuals: first stage. */	
	int I_s2 = 0;	/*!< Number of infected and infectious symptomatic individuals: second stage.  */	
	int I_s3 = 0;	/*!< Number of infected and infectious symptomatic individuals: thrid stage. */	
	int I_s4 = 0;	/*!< Number of infected and infectious symptomatic individuals: last stage.  */	
	int H = 0;		/*!< Number of infected individuals that are hospitalised.  */	
	int R = 0;		/*!< Number of infected individuals that are recovered from infection.   */	
	int D = 0;		/*!< Number of dead individuals due to disease. */	
};

/**
 * @brief Structure to hold status objects
 * 
 * Contains vectors containing information on the trajectories of the predicted disease status for individuals involved in the model 
 * (i.e. number of incident cases, number of incident deaths for each simulation step).
 */
struct Status
{
	std::vector<int> simulation;				/*!< Number of incident (newly detected/reported) cases of covid in each simulation step. */
	std::vector<int> deaths;					/*!< Overall number of incident deaths due to covid in each simulation step. */
	std::vector<int> hospital_deaths;			/*!< Number of incident deaths reported at hospital due to covid in each simulation step. */
	std::vector<Compartments> ends;				/*!< Population per epidemiological state and per age group on last day. */
};

/**
 * @brief Structure containing population counters after infection
 * 
 * Contains counters for the number of newly detected (due to testing or hospitalisation) cases, deaths and hospitalisations
 * at each time step.
 */
struct InfectionState
{
	int detected = 0;				/*!< Number of newly detected (due to testing or hospitalisation) cases at each time step. */
	int hospitalised = 0;			/*!< Number of newly detected cases due to their hospitalisation at each time step. */
	int deaths = 0;					/*!< Overall number of incident deaths due to covid at each time step.  */
	int hospital_deaths = 0;		/*!< Number of incident deaths reported at hospital due to covid at each time step. */
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
 * @brief Generate an indection spread and compute resulting populations
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
 * @param structure 
 */
InfectionState GenerateInfectionSpread(Random::RNGInterface::Sptr rng, Compartments& pop,
    const int& n_hospitalised, ::EERAModel::params fixed_parameters, 
    std::vector<double> parameter_set, std::vector<double> cfr_tab,
    double pf_val, double lambda,  ModelStructureId structure);

} // namespace Model
} // namespace EERAModel
