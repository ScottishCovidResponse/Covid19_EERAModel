/* Created on: 01 05 2020
 * Authors: Thibaud Porphyre
 *
 *
 *version 0.3.2.4
 *
 * version report: 
 *		  - change the whole model structure (changes in Lambda, infspread)
 *		  - fix issues in the selection processes
 * 		  - add/remove parameters to infer
 *        - add a void function _flow_ to simplify function _infspread_
 *  
 * fitting procedure: ABS-SMC
 * model to fit: spread, SEIIsHRD
 * number of herd: single
 * model type: stochastic, age-structured population, tau-leap
 *
 * time-step = 1 day
 *
 * selection measures: normalise sum squared error 
 *
 * fitted parameters: p_i, p_hcw, c_hcw, q and d, p_s, p_hf,rrdh,lambda
 *
 * main.cpp
 *
 *
 */


#include <iostream>
#include <time.h>
#include "ModelTypes.h"
#include "IO.h"
#include "Random.h"
#include "Model.h"
#include "PredictionFramework.h"
#include "InferenceFramework.h"

using namespace EERAModel;

int main() {

	const std::string out_dir = std::string(ROOT_DIR)+"/outputs";

	Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

	// Read in the model's input parameters
	std::cout << "PROJECT ROOT DIRECTORY:\t"+std::string(ROOT_DIR) << std::endl;
	const std::string params_addr = std::string(ROOT_DIR)+"/data/parameters.ini";

	ModelInputParameters modelInputParameters = IO::ReadParametersFromFile(params_addr, logger);

	(*logger) << "[Parameters File]:\n    " << params_addr << std::endl;
	
    // Read posterior particle parameters if run type is "Prediction"
	if (modelInputParameters.run_type == "Prediction")
	{
		const std::string posterior_params_addr = std::string(ROOT_DIR)+"/src/posterior_particle_params.csv";
		PosteriorParticleParameters PosteriorParticleParameters= IO::ReadPosteriorParametersFromFile(posterior_params_addr, logger);
		modelInputParameters.posterior_param_list = PosteriorParticleParameters.posterior_param_list;
	}

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

    if (modelInputParameters.run_type == "Prediction")
    {
        Prediction::PredictionFramework framework(modelInputParameters, observations, rng, out_dir, logger);

        int n_sim_steps = 100;
        std::vector<double> parameter_set(8, 0.0);

        // framework.Run(parameter_set, n_sim_steps);
		framework.Run(modelInputParameters.posterior_param_list, n_sim_steps);
    }
    else
    {
        Inference::InferenceFramework framework(modelInputParameters, observations, rng, out_dir, logger);
        
        framework.Run();
    }
	
}
