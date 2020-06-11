#include "PredictionFramework.h"
#include "Model.h"
#include "IO.h"

namespace EERAModel {
namespace Prediction {

PredictionFramework::PredictionFramework(
    const ModelInputParameters& modelInputParameters,
    const InputObservations& observations,
    Random::RNGInterface::Sptr rng,
    const std::string& outDir,
    Utilities::logging_stream::Sptr log)
     : seedlist_(modelInputParameters.seedlist),
       dayShut_(modelInputParameters.day_shut),
       modelStructure_(modelInputParameters.model_structure),
       rng_(rng),
       outDir_(outDir),
       log_(log)
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

void PredictionFramework::Run(std::vector<double> parameterSet, int nSimulationSteps)
{
    clock_t startTime = clock();

    Status status = Model::RunModel(
        parameterSet, fixedParameters_, ageGroupData_, seedlist_, dayShut_, ageNums_,
        nSimulationSteps, modelStructure_, rng_
    );

    double time_taken;
    time_taken = double(clock() - startTime)/(double)CLOCKS_PER_SEC;

    (*log_) << " <computation time> " << time_taken << " seconds.\n";

    // To do: Write model prediction outputs here
    std::vector< std::vector<int> > end_comps;
	end_comps = Model::compartments_to_vector(status.ends);

    IO::WritePredictionsToFiles(status, end_comps, outDir_, log_);
}

} // namespace Prediction
} // namespace EERAModel