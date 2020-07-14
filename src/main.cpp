#include <iostream>
#include <time.h>
#include "ModelTypes.h"
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

	// Read in the model's input parameters
	arg_parser.logArguments(logger);

	const std::string params_addr = std::string(ROOT_DIR)+"/data/parameters.ini";

	// ModelInputParameters modelInputParameters = IO::ReadParametersFromFile(params_addr, logger);
    SupplementaryInputParameters supplementaryParameters = IO::ReadSupplementaryParameters(params_addr, logger);
    arg_parser.AppendOptions(supplementaryParameters);

	(*logger) << "[Parameters File]:\n    " << params_addr << std::endl;

	// Read in the observations
	InputObservations observations = IO::ReadObservationsFromFiles(logger);

	// Set up the random number generator, deciding what kind of seed to use
	unsigned long randomiser_seed;
	if (supplementaryParameters.seedlist.use_fixed_seed) {
		randomiser_seed = supplementaryParameters.seedlist.seed_value;
	} else {
        randomiser_seed = time(nullptr);
	}
    Random::RNG::Sptr rng = std::make_shared<Random::RNG>(randomiser_seed);
    IO::LogRandomiserSettings(supplementaryParameters, randomiser_seed, logger);

    CommonModelInputParameters commonParameters = IO::ReadCommonParameters(params_addr);

    // Log the fixed parameters
    IO::LogFixedParameters(commonParameters, logger);

    // Log the disease seed settings
    IO::LogSeedSettings(supplementaryParameters.seedlist, logger);
    
    // Select the model structure to use
    Model::ModelInterface::Sptr model;
    if (ModelStructureId::ORIGINAL == supplementaryParameters.model_structure)
    {
        model = std::make_shared<Model::OriginalModel>(commonParameters, observations, rng, logger);
    }
    else if (ModelStructureId::IRISH == supplementaryParameters.model_structure)
    {
        model = std::make_shared<Model::IrishModel>(commonParameters, observations, rng, logger);
    }
	else
	{
		model = std::make_shared<Model::TempModel>(commonParameters, observations, rng, logger);
	}

    // Select the mode to run in - prediction or inference    
    if (ModelModeId::PREDICTION == supplementaryParameters.run_type)
    {
        std::string configDir(std::string(ROOT_DIR) + "/data");
        PredictionConfig predictionConfig = IO::ReadPredictionConfig(configDir, logger);
        
        IO::LogPredictionConfig(predictionConfig, logger);

        Prediction::PredictionFramework framework(model, predictionConfig, 
            rng, out_dir, logger);

		framework.Run();
    }
    else
    {
        InferenceConfig inferenceConfig = IO::ReadInferenceConfig(params_addr, logger);

        Inference::InferenceFramework framework(model, inferenceConfig, observations, rng, out_dir, logger);
        
        framework.Run();
    }
}
