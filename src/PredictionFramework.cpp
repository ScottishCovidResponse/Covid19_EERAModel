#include "PredictionFramework.h"
#include "ModelCommon.h"
#include "Observations.h"
#include "IO.h"

namespace EERAModel {
namespace Prediction {

PredictionFramework::PredictionFramework(
    Model::ModelInterface::Sptr model,
    const ModelInputParameters& modelInputParameters,
    const InputObservations& observations,
    Random::RNGInterface::Sptr rng,
    const std::string& outDir,
    Utilities::logging_stream::Sptr log)
     : model_(model),
       modelInputParameters_(modelInputParameters),
       observations_(observations),
       rng_(rng),
       outDir_(outDir),
       log_(log) {}

void PredictionFramework::Run(std::vector<double> parameterSet, int nSimulationSteps)
{
    clock_t startTime = clock();

    (*log_) << "[Settings]:\n";
    (*log_) << "    number of parameters tested: "<< modelInputParameters_.nPar << std::endl;
    (*log_) << "    seeding method: "<< modelInputParameters_.seedlist.seedmethod<<  std::endl;
	if (modelInputParameters_.seedlist.seedmethod == "random"){
		(*log_) << "    number of seed: " << modelInputParameters_.seedlist.nseed << std::endl;
	} else if(modelInputParameters_.seedlist.seedmethod == "background"){
		(*log_) << "    duration of the high risk period (hrp): " << modelInputParameters_.seedlist.hrp << std::endl;
	}

    Status status = model_->Run(parameterSet, modelInputParameters_.seedlist,
        modelInputParameters_.day_shut, nSimulationSteps);

    double time_taken = double(clock() - startTime)/(double)CLOCKS_PER_SEC;

    (*log_) << "\n <computation time> " << time_taken << " seconds.\n";

    std::vector< std::vector<int> > end_comps;
	end_comps = Model::compartments_to_vector(status.ends);

    IO::WritePredictionsToFiles(status, end_comps, outDir_, log_);
}

} // namespace Prediction
} // namespace EERAModel