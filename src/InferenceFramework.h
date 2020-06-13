#pragma once

#include "ModelTypes.h"
#include "Random.h"
#include "Utilities.h"

namespace EERAModel {
namespace Inference {

/**
 * @class InferenceFramework
 * @brief Parameter inference framework
 */
class InferenceFramework
{
public:
    /**
     * @brief Framework constructor
     *
     * @param modelInputParameters Model input parameters
     * @param observations Observations
     * @param rng Seeded random number generator
     * @param outDir Outputs directory path
     * @param log Logger
     */
    InferenceFramework(const ModelInputParameters& modelInputParameters,
        const InputObservations& observations,
        Random::RNGInterface::Sptr rng,
        const std::string& outDir,
        Utilities::logging_stream::Sptr log);

    /**
     * @brief Run the model within the inference framework
     */
     void Run();

private:
    /**
     * @private
     * @brief Run the model inside the inference framework
     */
    void ModelSelect(::EERAModel::particle &outvec, const std::vector<params>& fixed_parameters,
        const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
        const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
        const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
        int duration, seed seedlist, int day_shut, Random::RNGInterface::Sptr rng, const std::vector<int>& obsHosp,
        const std::vector<int>& obsDeaths, ModelStructureId structure);

    /**
     * @private
     * @brief Model input parameters
     */
    ModelInputParameters modelInputParameters_;

    /**
     * @private
     * @brief Model input observations
     */
    InputObservations observations_;

    /**
     * @private
     * @brief Random number generator
     */
    Random::RNGInterface::Sptr rng_;

    /**
     * @private
     * @brief Outputs directory path
     */
    std::string outDir_;

    /**
     * @private
     * @brief Logger
     */
    Utilities::logging_stream::Sptr log_;
};

/**
 * @brief Compute the Kernel Window
 * 
 * @param nPar Number of parameters to compute kernel for
 * @param particleList List of previously accepted particles
 * @param kernelFactor Common kernel multiplicative factor
 * @param vlimitKernel Storage for the computed kernel
 * @param vect_Max Storage for maximum values of ranges
 * @param vect_Min Storage for minimum values of ranges
 */
void ComputeKernelWindow(int nPar, const std::vector<particle>& particleList,
	double kernelFactor, std::vector<double>& vlimitKernel, std::vector<double>& vect_Max, 
	std::vector<double>& vect_Min);

/**
 * @brief Compute weight distribution
 * 
 * Create the discrete distribution of the weights for the "importance sampling" process
 * 
 * @param particleList List of previously accepted particles
 * 
 * @return Weight distribution
 */
std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList);


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

} // namespace Inference
} // namespace EERAModel
