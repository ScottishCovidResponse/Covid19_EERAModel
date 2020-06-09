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
    fixedParameters_ = Model::BuildFixedParameters(size, modelInputParameters.paramlist);
}

Status Framework::RunModel(EERAModel::params parameterSet, int nSimulationSteps)
{
    return modelRunner_(parameterSet, fixedParameters_, ageGroupData_, seedlist_, dayShut_,
        ageNums_, nSimulationSteps, modelStructure_, rng_);
}

} // namespace Prediction
} // namespace EERAModel