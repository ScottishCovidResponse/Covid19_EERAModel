#include "PredictionFramework.h"

namespace EERAModel {
namespace Prediction {

Framework::Framework(ModelRunner modelRunner,
    ModelStructureId modelStructure,
    unsigned int size,
    AgeGroupData ageGroupData,
    const ModelInputParameters& modelInputParameters,
    int dayShut,
    std::vector<int> ageNums,
    Random::RNGInterface::Sptr rng)
        : modelRunner_(modelRunner),
        ageGroupData_(ageGroupData),
        seedlist_(modelInputParameters.seedlist),
        dayShut_(dayShut),
        ageNums_(ageNums),
        modelStructure_(modelStructure),
        rng_(rng)
{
    fixed_parameters_ = Model::BuildFixedParameters(size, modelInputParameters.paramlist);
}

} // namespace Prediction
} // namespace EERAModel