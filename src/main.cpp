#include <iostream>
#include <time.h>
#include "ModelTypes.h"
#include "IO.h"
#include "Random.h"
#include "DataSourcing.h"
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

	Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(arg_parser.outputDir());

	// Read in the model's input parameters
	arg_parser.logArguments(logger);

    // Ordering of enum class SourceID important for this line as here we convert boolean isLocal to
    // the SourceID LOCAL/REMOTE
    DataSourcing::DataSource data_source = DataSourcing::getSource( DataSourcing::SourceID(arg_parser.runLocal()),
                                                                        logger,         
                                                                         arg_parser.localSourceDir()+"/data" );

    arg_parser.AppendOptions(data_source);

	// Set up the random number generator, deciding what kind of seed to use
	unsigned long randomiser_seed;
	if (data_source.getInputParameters().seedlist.use_fixed_seed) {
		randomiser_seed = data_source.getInputParameters().seedlist.seed_value;
	} else {
        randomiser_seed = time(nullptr);
	}
    Random::RNG::Sptr rng = std::make_shared<Random::RNG>(randomiser_seed);
    (*logger) << "[Seed]:\n    Type: ";
    (*logger) << ((data_source.getInputParameters().seedlist.use_fixed_seed) ? "Fixed" : "Time based") << std::endl;
	(*logger) << "    Value: " << randomiser_seed << std::endl;

    // Select the model structure to use
    Model::ModelInterface::Sptr model;
    if (ModelStructureId::ORIGINAL == data_source.getInputParameters().model_structure)
    {
        model = std::make_shared<Model::OriginalModel>(rng);
    }
    else if (ModelStructureId::IRISH == data_source.getInputParameters().model_structure)
    {
        model = std::make_shared<Model::IrishModel>(rng);
    }
	else
	{
		model = std::make_shared<Model::TempModel>(rng);
	}

    // Select the mode to run in - prediction or inference    
    if (data_source.getInputParameters().run_type == ModelModeId::PREDICTION)
    {

        Prediction::PredictionFramework framework(model, data_source.getInputParameters(), 
                    data_source.getInputObservations(), rng, arg_parser.outputDir(), logger);

        int n_sim_steps = 100;

		framework.Run(data_source.getInputParameters().posterior_param_list, n_sim_steps);
    }
    else
    {
        Inference::InferenceFramework framework(model, data_source.getInputParameters(), 
                                                data_source.getInputObservations(), rng, arg_parser.outputDir(), logger);
        
        framework.Run();
    }
	
}
