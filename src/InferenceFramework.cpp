#include "InferenceFramework.h"

namespace EERAModel {
namespace Inference {

InferenceFramework::InferenceFramework(const ModelInputParameters& modelInputParameters,
    const InputObservations& observations,
    Random::RNGInterface::Sptr rng,
    Utilities::logging_stream::Sptr log) {}

void InferenceFramework::Run() {}

} // namespace Inference
} // namespace EERAModel