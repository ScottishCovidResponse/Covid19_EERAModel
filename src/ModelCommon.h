#pragma once

#include "ModelTypes.h"
#include "Random.h"
#include <vector>
#include <memory>

/**
 * @brief Namespace for the EERA COVID-19 Model and member functions
 *
 * This namespace includes functions and variables used to access data,
 * initiate, run and evaluate the model.
 */
namespace EERAModel {
/**
 * @brief Namespace containing the model definitions and execution
 *
 * This namespace contains functions which run the model using the given
 * setup/parameters.
 */
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
    virtual Status Run(std::vector<double> parameter_set, seed& seedlist, int day_shut, int n_sim_steps) = 0;
};

/**
 * @class ModelParameters
 * @brief Describes the parameters used by the model
 * 
 * Provides constants associated with the use of model parameter sets. Model parameter sets are
 * stored as vectors, where each element of the vector corresponds to a particular parameter. In
 * order to avoid the use of lots of magic numbers in the code, this class provides index constants
 * with meaningful names.
 * 
 * @note Model parameters are distinct from "fixed parameters" to the model; the latter are also
 * used to configure the behaviour of the model, but are not varies and inferred as part of the
 * inference framework. The former are also referred to as "Inference Parameters" when speaking in 
 * the context of the inference framework.
 */
class ModelParameters 
{
public:
    /** @brief Enumeration for different parameter indices */
    using Enum = unsigned int;

    /** @brief Probability of infection (non-HCW) */
    static constexpr Enum PINF   = 0; 
    /** @brief Probability of infection (HCW) */
    static constexpr Enum PHCW   = 1;
    /** @brief Mean number of HCW contacts per day */
    static constexpr Enum CHCW   = 2;
    /** @brief Proportion of population observing social distancing */
    static constexpr Enum D      = 3;
    /** @brief Proportion of contacts made by self-isolaters */
    static constexpr Enum Q      = 4;
    /** @brief Age-dependent probability of developing symptoms */
    static constexpr Enum PS     = 5;
    /** @brief Risk of death if not hospitalised */
    static constexpr Enum RRD    = 6;
    /** @brief Background transmission rate */
    static constexpr Enum LAMBDA = 7;

    /** @brief Number of model parameters */
    static constexpr unsigned int NPARAMS = 8;
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
int GetPopulationOfRegion(const ObservationsForModels& obs, int region_id);

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
int ComputeNumberOfHCWInRegion(int regionalPopulation, int totalHCW, const ObservationsForModels& observations);

/**
 * @brief Compute the agenums 
 * 
 * @param Npop Population
 * @param Nhcw Number of health care workers
 * @param obs Model observations
 * 
 * @return agenums
 */
std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const ObservationsForModels& obs);

/**
 * @brief Build the fixed parameters data structure
 * 
 * @param size Number of instances of the parameters to generate
 * @param parameters Source parameter list
 * 
 * @return Vector of @p size copies of @p parameters
 */
std::vector<params> BuildFixedParameters(unsigned int size, params parameters);

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
std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec);

} // namespace Model
} // namespace EERAModel