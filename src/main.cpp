/* Created on: 17 04 2020
 * Authors: Thibaud Porphyre
 *
 *
 *version 0.3.2
 *
 * version report: 
 *		  - fix the transmission flow between Is,H, D and R relative to %frail
 * 
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
 * fitted parameters: p_i, p_hcw, c_hcw, q and d
 *
 * main.cpp
 *
 *
 */

#include <random>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>

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

#include "Model.h"
#include "ModelTypes.h"
#include "IO.h"
#include "Utilities.h"

using namespace EERAModel;

int main() {

	const std::string out_dir = std::string(ROOT_DIR)+"/outputs";

	Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

	// Read in the model's input parameters
	std::cout << "PROJECT ROOT DIRECTORY:\t"+std::string(ROOT_DIR) << std::endl;
	const std::string params_addr = std::string(ROOT_DIR)+"/data/parameters.ini";

	ModelInputParameters modelInputParameters = IO::ReadParametersFromFile(params_addr, logger);

	(*logger) << "[Parameters File]:\n    " << params_addr << std::endl;

	// Read prior particle parameters if run type is "Prediction"
	if (modelInputParameters.run_type == "Prediction")
	{
		const std::string prior_params_addr = std::string(ROOT_DIR)+"data/prior_particle_params.csv";
		PriorParticleParameters priorParticleParameters= IO::ReadPriorParametersFromFile(prior_params_addr, logger);
		modelInputParameters.prior_param_list = priorParticleParameters.prior_param_list;
	}

	// Read in the observations
	InputObservations observations = IO::ReadObservationsFromFiles(logger);

	// Decide which kind of seed to use
	unsigned long randomiser_seed;
	if (modelInputParameters.seedlist.use_fixed_seed) {
		randomiser_seed = modelInputParameters.seedlist.seed_value;
	} else {
        randomiser_seed = time(NULL);
	}
	//initialise the gsl random number generator with a seed depending on the time of the run
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, randomiser_seed);

	//initialise the random number generator for importance sampling
	std::mt19937 gen(randomiser_seed);


	(*logger) << "[Seed]:\n    Type: ";
        (*logger) << ((modelInputParameters.seedlist.use_fixed_seed) ? "Fixed" : "Time based") << std::endl;
	(*logger) << "    Value: " << randomiser_seed << std::endl;

	Model::Run(modelInputParameters, observations, r, gen, out_dir, logger);
}
