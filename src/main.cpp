#include <iostream>
#include <time.h>
#include "ModelTypes.h"
#include "IO-datapipeline.h"
#include "IO.h"
#include "Random.h"
#include "ModelCommon.h"
#include "OriginalModel.h"
#include "IrishModel.h"
#include "TempModel.h"
#include "PredictionFramework.h"
#include "InferenceFramework.h"
#include "ArgumentParser.h"

using namespace EERAModel;

int main(int argc, char** argv) 
{
    ArgumentParser arg_parser(argc, argv);

    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";

    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);
    IO::LogGitVersionInfo(logger);
    IO::LogDependencyVersionInfo(logger);
    arg_parser.logArguments(logger);

    // Parameters ini file location
    const std::string params_addr = std::string(ROOT_DIR)+"/data/parameters.ini";

    // Model Observational data location
    std::string modelConfigDir(std::string(ROOT_DIR) + "/data");

    // Get the run time
    std::string timeStamp = logger->getLoggerTime();
    (*logger) << "[Run Time]:\n    " << timeStamp << "\n";

    // Start the data pipeline
    pybind11::scoped_interpreter guard{}; // start the interpreter and keep it alive
    IO::IOdatapipeline datapipeline{
        params_addr, modelConfigDir, out_dir, logger, arg_parser.getArgs().datapipeline_path, timeStamp};

    SupplementaryInputParameters supplementaryParameters = IO::ReadSupplementaryParameters(params_addr, logger);
    arg_parser.AppendOptions(supplementaryParameters);

    (*logger) << "[Parameters File]:\n    " << params_addr << std::endl;

    // Set up the random number generator, deciding what kind of seed to use
    unsigned long randomiser_seed;
    if (supplementaryParameters.seedlist.use_fixed_seed) {
        randomiser_seed = supplementaryParameters.seedlist.seed_value;
    } else {
        randomiser_seed = time(nullptr);
    }
    Random::RNG::Sptr rng = std::make_shared<Random::RNG>(randomiser_seed);
    IO::LogRandomiserSettings(supplementaryParameters, randomiser_seed, logger);

    (*logger) << "[Common Parameters]:" << std::endl;

    // Import common parameters for all models
    CommonModelInputParameters commonParameters = datapipeline.ReadCommonParameters();
    // CommonModelInputParameters commonParameters = IO::ReadCommonParameters(params_addr);

    (*logger) << "[Model Observations]:" << std::endl;

    // Import model observational data
    ObservationsForModels modelObservations = datapipeline.ReadModelObservations();
    // ObservationsForModels modelObservations = IO::ReadModelObservations(modelConfigDir, logger);

    // Log the disease seed settings
    IO::LogSeedSettings(supplementaryParameters.seedlist, logger);

    // Select the model structure to use
    Model::ModelInterface::Sptr model;
    if (ModelStructureId::ORIGINAL == supplementaryParameters.model_structure)
    {
        model = std::make_shared<Model::OriginalModel>(commonParameters, modelObservations, rng, logger);
    }
    else if (ModelStructureId::IRISH == supplementaryParameters.model_structure)
    {
        model = std::make_shared<Model::IrishModel>(commonParameters, modelObservations, rng, logger);
    }
    else
    {
        model = std::make_shared<Model::TempModel>(commonParameters, modelObservations, rng, logger);
    }

    // Select the mode to run in - prediction or inference    
    if (ModelModeId::PREDICTION == supplementaryParameters.run_type)
    {
        // std::string configDir(std::string(ROOT_DIR) + "/data");
        int index = arg_parser.parameterSetIndex();

        PredictionConfig predictionConfig = datapipeline.ReadPredictionConfig(index, commonParameters);
        IO::LogPredictionConfig(predictionConfig, logger);

        // Update the model with the fixed parameters from the prediction configuration
        model->SetFixedParameters(predictionConfig.fixedParameters);
        
        Prediction::PredictionFramework framework(model, predictionConfig, rng, logger, &datapipeline);

        framework.Run();
    }
    else
    {
        // std::string configDir(std::string(ROOT_DIR) + "/data");
        InferenceConfig inferenceConfig = datapipeline.ReadInferenceConfig(commonParameters, modelObservations);

        IO::LogFixedParameters(commonParameters.paramlist, logger);

        Inference::InferenceFramework framework(model, inferenceConfig, rng, logger, &datapipeline);
        
        framework.Run();
    }
}
