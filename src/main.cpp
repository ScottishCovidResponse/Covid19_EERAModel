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

	ModelInputParameters modelInputParameters = IO::ReadParametersFromFile(params_addr, logger);
    arg_parser.AppendOptions(modelInputParameters);

	(*logger) << "[Parameters File]:\n    " << params_addr << std::endl;
	
	// Read in the observations
	InputObservations observations = IO::ReadObservationsFromFiles(logger);

	// Set up the random number generator, deciding what kind of seed to use
	unsigned long randomiser_seed;
	if (modelInputParameters.seedlist.use_fixed_seed) {
		randomiser_seed = modelInputParameters.seedlist.seed_value;
	} else {
        randomiser_seed = time(nullptr);
	}
    Random::RNG::Sptr rng = std::make_shared<Random::RNG>(randomiser_seed);
    (*logger) << "[Seed]:\n    Type: ";
    (*logger) << ((modelInputParameters.seedlist.use_fixed_seed) ? "Fixed" : "Time based") << std::endl;
	(*logger) << "    Value: " << randomiser_seed << std::endl;

    // Select the model structure to use
    Model::ModelInterface::Sptr model;
    if (ModelStructureId::ORIGINAL == modelInputParameters.model_structure)
    {
        model = std::make_shared<Model::OriginalModel>(rng);
    }
    else if (ModelStructureId::IRISH == modelInputParameters.model_structure)
    {
        model = std::make_shared<Model::IrishModel>(rng);
    }
	else
	{
		model = std::make_shared<Model::TempModel>(rng);
	}

    // Select the mode to run in - prediction or inference    
    if (modelInputParameters.run_type == ModelModeId::PREDICTION)
    {
        const std::string posterior_params_addr = std::string(ROOT_DIR)+"/data/example_posterior_parameter_sets.txt";

        modelInputParameters.posterior_param_list = 
			IO::ReadPosteriorParametersFromFile(posterior_params_addr, 
				modelInputParameters.posterior_parameter_select);

        Prediction::PredictionFramework framework(model, modelInputParameters, observations, rng, out_dir, logger);

        int n_sim_steps = 100;

		framework.Run(modelInputParameters.posterior_param_list, n_sim_steps);
    }
    else
    {
        Inference::InferenceFramework framework(model, modelInputParameters, observations, rng, out_dir, logger);
        
        framework.Run();
    }
}
