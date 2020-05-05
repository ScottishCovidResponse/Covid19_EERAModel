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

#include "Model.h"
#include "ModelTypes.h"
#include "IO.h"

using namespace EERAModel;

int main(int argc, char **argv) {

	// Read in the model's input parameters
	ModelInputParameters modelInputParameters = IO::ReadParametersFromFile("./data/parameters.ini");

	// Read in the observations
	Observations observations = IO::ReadObservationsFromFiles();

	//initialise the gsl random number generator with a seed depending on the time of the run
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(NULL));

	//initialise the random number generator for importance sampling
	std::mt19937 gen(time(NULL));

	Model::Run(modelInputParameters, observations, r, gen, "./outputs");
}
