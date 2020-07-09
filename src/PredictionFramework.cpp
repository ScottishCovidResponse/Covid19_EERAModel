#include "PredictionFramework.h"
#include "ModelCommon.h"
#include "Observations.h"
#include "IO.h"

namespace EERAModel {
namespace Prediction {

PredictionFramework::PredictionFramework(
    Model::ModelInterface::Sptr model,
    const ModelInputParameters& modelInputParameters,
    const PredictionConfig& config,
    Random::RNGInterface::Sptr rng,
    const std::string& outDir,
    Utilities::logging_stream::Sptr log)
     : model_(model),
       modelInputParameters_(modelInputParameters),
       config_(config),
       rng_(rng),
       outDir_(outDir),
       log_(log) {}

void PredictionFramework::Run()
{
    clock_t startTime = clock();

    Status status = model_->Run(config_.posterior_parameters, modelInputParameters_.seedlist,
        modelInputParameters_.day_shut, config_.n_sim_steps);

    double time_taken = static_cast<double>(clock() - startTime)/static_cast<double>(CLOCKS_PER_SEC);

    (*log_) << "\n <computation time> " << time_taken << " seconds.\n";

    std::vector<std::vector<int>> end_comps = Model::compartments_to_vector(status.ends);

    IO::WritePredictionsToFiles(status, end_comps, outDir_, log_);
}

} // namespace Prediction
} // namespace EERAModel