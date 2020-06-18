#pragma once

#include "ModelCommon.h"
#include "ModelTypes.h"
#include "Random.h"
#include "Utilities.h"
#include <vector>
#include <functional>

namespace EERAModel {
namespace Prediction {

/**
 * @class Framework
 * @brief Model prediction framework
 */
class PredictionFramework
{
public:
    /**
     * @brief Framework constructor
     *
     * @param model Model interface
     * @param modelInputParameters Model input parameters
     * @param observations Observations
     * @param rng Seeded random number generator
     * @param log Logger
     */
    PredictionFramework(Model::ModelInterface::Sptr model,
        const ModelInputParameters& modelInputParameters,
        const InputObservations& observations,
        Random::RNGInterface::Sptr rng,
        Utilities::logging_stream::Sptr log);

    /**
     * @brief Run the model within the prediction framework
     * 
     * @param parameterSet Set of model parameters to use
     * @param nSimulationSteps Number of simulation steps to run
     * 
     * @param Status at the end of the model run
     */
     void Run(std::vector<double> parameterSet, int nSimulationSteps);

private:
    
    /**
     * @private
     * @brief Model interface
     */
    Model::ModelInterface::Sptr model_;

    /**
     * @private
     * @brief Fixed model parameters
     */
    std::vector<EERAModel::params> fixedParameters_;

    /**
     * @private 
     * @brief Distribution of population amongst different age groups
     */
    AgeGroupData ageGroupData_;

    /**
     * @private
     * @brief TBC
     */
    seed seedlist_;

    /**
     * @private
     * @brief Day at which the lockdown took effect
     */
    int dayShut_;

    std::vector<int> ageNums_;

    /**
     * @private
     * @brief Which model structure to use
     */
    ModelStructureId modelStructure_;

    /**
     * @private
     * @brief Random number generator
     */
    Random::RNGInterface::Sptr rng_;
};

} // namespace Prediction
} // namespace EERAModel