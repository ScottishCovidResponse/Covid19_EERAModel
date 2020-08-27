#include "PredictionFramework.h"
#include "ModelCommon.h"
#include "Observations.h"
#include "IO.h"
#include "Timer.h"

namespace EERAModel {
namespace Prediction {

PredictionFramework::PredictionFramework(
    Model::ModelInterface::Sptr model,
    const PredictionConfig& config,
    Random::RNGInterface::Sptr rng,
    const std::string& outDir,
    Utilities::logging_stream::Sptr log)
     : model_(model),
       config_(config),
       rng_(rng),
       outDir_(outDir),
       log_(log) {}

void PredictionFramework::Run()
{
    SimpleTimer timer;

    Status status = model_->Run(config_.posterior_parameters, config_.seedlist,
        config_.day_shut, config_.n_sim_steps);

    (*log_) << "\n <computation time> " << timer.elapsedTime() << " seconds.\n";

    std::vector<std::vector<int>> end_comps = Model::compartments_to_vector(status.ends);

    IO::WritePredictionsToFiles(status, end_comps, outDir_, log_);
}

} // namespace Prediction
} // namespace EERAModel