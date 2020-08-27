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
    std::vector<Status> statuses(config_.n_iterations);
    SimpleTimer timer;

    (*log_) << "[Simulations]:" << std::endl;
    for (int iter = 0; iter < config_.n_iterations; ++iter) {
        statuses[iter] = model_->Run(config_.posterior_parameters, config_.seedlist,
            config_.day_shut, config_.n_sim_steps);
        
        // Output a marker for every 10 iterations run
        if (iter > 0 && iter % 10 == 0) (*log_) << "|" << std::flush;
    }
    
    (*log_) << "\n <computation time> " << timer.elapsedTime() << " seconds.\n";

    IO::WritePredictionsToFiles(statuses, outDir_, log_);
}

} // namespace Prediction
} // namespace EERAModel