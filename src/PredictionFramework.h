#pragma once

#include "ModelCommon.h"
#include "ModelTypes.h"
#include "Random.h"
#include "Utilities.h"
#include "IO-datapipeline.h"
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
     * @param config Prediction configuration
     * @param rng Seeded random number generator
     * @param log Logger
     * @param datapipeline data pipeline handler
     */
    PredictionFramework(Model::ModelInterface::Sptr model,
        const PredictionConfig& config,
        Random::RNGInterface::Sptr rng,
        Utilities::logging_stream::Sptr log,
        IO::IOdatapipeline *datapipeline);

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
     * @brief Logger
     */
    Utilities::logging_stream::Sptr log_;

    /**
     * @private
     * @brief Data pipeline handler
     */
    IO::IOdatapipeline *datapipeline_;
};

} // namespace Prediction
} // namespace EERAModel