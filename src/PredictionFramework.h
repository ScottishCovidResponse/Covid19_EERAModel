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
     * @param config Prediction configuration
     * @param rng Seeded random number generator
     * @param log Logger
     */
    PredictionFramework(Model::ModelInterface::Sptr model,
        const PredictionConfig& config,
        Random::RNGInterface::Sptr rng,
        const std::string& outDir,
        Utilities::logging_stream::Sptr log);

    /**
     * @brief Run the model within the prediction framework
     *  
     * @param Status at the end of the model run
     */
     void Run();

private:
    
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

    /**
     * @private
     * @brief Prediction configuration
     */
    PredictionConfig config_;

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

} // namespace Prediction
} // namespace EERAModel