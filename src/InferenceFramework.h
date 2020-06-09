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
     * @param log Logger
     */
    InferenceFramework(const ModelInputParameters& modelInputParameters,
        const InputObservations& observations,
        Random::RNGInterface::Sptr rng,
        Utilities::logging_stream::Sptr log);

    /**
     * @brief Run the model within the inference framework
     */
     void Run();

private:
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
     * @brief Logger
     */
    Utilities::logging_stream::Sptr log_;
};

} // namespace Inference
} // namespace EERAModel