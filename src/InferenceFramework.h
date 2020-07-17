#pragma once

#include "ModelCommon.h"
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
     * @param model Model interface
     * @param inferenceConfig Inference config
     * @param observations Observations
     * @param rng Seeded random number generator
     * @param outDir Outputs directory path
     * @param log Logger
     */
    InferenceFramework(Model::ModelInterface::Sptr model,
        const InferenceConfig& inferenceConfig,
        // const InputObservations& observations,
        Random::RNGInterface::Sptr rng,
        const std::string& outDir,
        Utilities::logging_stream::Sptr log);
    
    /**
     * @brief Calculate the time offset for the dataset
     * 
     * Determines the offset for the start of the dataset based on
     * the given parameters
     * 
     * @param modelInputParameters model parameters including seeding options
     * @param log Logging stream
     * 
     * @return integer offset value
     */
    static int GetTimeOffSet(const InferenceConfig& inferenceConfig);

    /**
     * @brief Run the model within the inference framework
     */
     void Run();

private:
    /**
     * @private
     * @brief Run the model inside the inference framework
     */
    void ModelSelect(EERAModel::particle& outvec, const int& n_sim_steps, seed seedlist,
        int day_shut, const std::vector<int>& obsHosp, const std::vector<int>& obsDeaths);
    /**
     * @private
     * @brief Model interface
     */
    Model::ModelInterface::Sptr model_;

    // /**
    //  * @private
    //  * @brief Model input parameters
    //  */
    // ModelInputParameters modelInputParameters_;

    InferenceConfig inferenceConfig_;

    // /**
    //  * @private
    //  * @brief Model input observations
    //  */
    // InputObservations observations_;

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


} // namespace Inference
} // namespace EERAModel
