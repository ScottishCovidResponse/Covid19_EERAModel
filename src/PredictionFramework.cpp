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

    // Duplicated from Model.cpp - needs refactoring
    int scotlandPopulation = 0;
	for (unsigned int region = 0; region < observations.cases.size() - 1; ++region) {
		scotlandPopulation += observations.cases[region][0];
	}
	double regionalProportion = static_cast<double>(regionalPopulation) / scotlandPopulation;
	int healthCareWorkers = round(modelInputParameters.totN_hcw * regionalProportion); 
    
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