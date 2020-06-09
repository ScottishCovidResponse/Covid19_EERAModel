#include "PredictionFramework.h"
#include "Model.h"

namespace EERAModel {
namespace Prediction {

Framework::Framework(
    const ModelInputParameters& modelInputParameters,
    const InputObservations& observations,
    Random::RNGInterface::Sptr rng,
    Utilities::logging_stream::Sptr log)
     : seedlist_(modelInputParameters.seedlist),
       dayShut_(modelInputParameters.day_shut),
       modelStructure_(modelInputParameters.model_structure),
       rng_(rng)
{
    fixedParameters_ = Model::BuildFixedParameters(
        observations.waifw_norm.size(), modelInputParameters.paramlist
    );
        
    ageGroupData_ = AgeGroupData{
            observations.waifw_norm,
            observations.waifw_home,
            observations.waifw_sdist,
            observations.cfr_byage,
            observations.pf_pop[modelInputParameters.herd_id - 1]
    };
    
    int regionalPopulation = Model::GetPopulationOfRegion(
        observations, modelInputParameters.herd_id
    );

    int healthCareWorkers = Model::ComputeNumberOfHCWInRegion(
        regionalPopulation, modelInputParameters.totN_hcw, observations
    );
    
    ageNums_ = Model::ComputeAgeNums(
        modelInputParameters.herd_id, regionalPopulation, healthCareWorkers, observations
    );
}

void Framework::Run(std::vector<double> parameterSet, int nSimulationSteps)
{
    Status status = Model::RunModel(
        parameterSet, fixedParameters_, ageGroupData_, seedlist_, dayShut_, ageNums_,
        nSimulationSteps, modelStructure_, rng_
    );

    // To do: Write model status outputs here
}

} // namespace Prediction
} // namespace EERAModel