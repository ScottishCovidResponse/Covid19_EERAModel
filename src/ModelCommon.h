#pragma once

#include "ModelTypes.h"
#include "Random.h"
#include <vector>
#include <memory>

namespace EERAModel {
namespace Model {

/**
 * @class ModelInterface
 * @brief Abstract interface to a model
 */
class ModelInterface 
{
public:
    using Sptr = std::shared_ptr<ModelInterface>;

    virtual ~ModelInterface() = default;

    /**
     * @brief Run interface
     * 
     * Interface to model-specific implementation
     */
    virtual Status Run(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
				int n_sim_steps, ModelStructureId structure, Random::RNGInterface::Sptr rng) = 0;
};

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
 * @brief Accumulate total across all compartments
 * 
 * Totals all counters within every compartment in a Compartments struct
 * 
 * @param comp Compartments object
 * 
 * @return Total across all compartments
 */
int accumulate_compartments(const Compartments& comp);

/**
 * @brief Get the population of a region
 * 
 * Within the input cases data, the total population of the region is located in the first column
 * of the row corresponding to the region
 * 
 * @param obs Input observations
 * @param region_id Index of the region within the observation data set
 * 
 * @return The population of the region
 */
int GetPopulationOfRegion(const InputObservations& obs, int region_id);

/**
 * @brief Compute number of health care workers (HCW) in the region
 * 
 * Calculates the number of HCW in the region, assuming that the ratio of the region's population
 * to that of the country is same as that for HCWs only
 * 
 * @param regionalPopulation Population of the region
 * @param totalHCW Total number of HCWs in the country
 * @param observations Input observations
 */
int ComputeNumberOfHCWInRegion(int regionalPopulation, int totalHCW, const InputObservations& observations);

/**
 * @brief Compute the agenums 
 * 
 * @param Npop Population
 * @param Nhcw Number of health care workers
 * @param obs Model observations
 * 
 * @return agenums
 */
std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const InputObservations& obs);

/**
 * @brief Build the fixed parameters data structure
 * 
 * @param size Number of instances of the parameters to generate
 * @param parameters Source parameter list
 * 
 * @return Vector of @p size copies of @p parameters
 */
std::vector<params> BuildFixedParameters(unsigned int size, params parameters);

} // namespace Model
} // namespace EERAModel