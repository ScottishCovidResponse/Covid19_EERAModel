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
    int GetTimeOffSet(const ModelInputParameters& modelInputParameters);

    /**
     * @brief Run the model within the inference framework
     */
     void Run();

private:
    /**
     * @private
     * @brief Run the model inside the inference framework
     */
    void ModelSelect(EERAModel::particle& outvec, const std::vector<params>& fixed_parameters,
	    const AgeGroupData& per_age_data, std::vector <int> agenums, const int& n_sim_steps, 
	    seed seedlist, int day_shut, Random::RNGInterface::Sptr rng, const std::vector<int>& obsHosp, const std::vector<int>& obsDeaths, ModelStructureId structure);

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

} // namespace Inference
} // namespace EERAModel