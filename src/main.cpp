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


#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include <omp.h>
#include <random>
#include <array>

#include "ModelTypes.h"
#include "IO.h"
#include "Utilities.h"
#include "DistanceComputation.h"
#include "Model.h"
#include "FittingProcess.h"

using namespace std;
using namespace EERAModel;
using namespace EERAModel::Utilities;
using namespace EERAModel::IO;
using namespace EERAModel::DistanceComputation;
using namespace EERAModel::Model;
using namespace EERAModel::FittingProcess;

int main(int argc, char **argv) {

	/*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;

	// declare vectors of inputs
	vector<int > obsHosp;
	vector<int > obsDeaths;

	// variable declaration
	
	double valMax, valmin;

	// Read in the model's input parameters
	ModelInputParameters modelInputParameters = ReadParametersFromFile("./data/parameters.ini");

	// Read in the observations
	Observations observations = ReadObservationsFromFiles();

    cout << "[Settings]:\n";
	cout<< "number of parameters tested: "<< modelInputParameters.nPar <<endl;
    cout<< "seeding method: "<< modelInputParameters.seedlist.seedmethod<<endl;
	if(modelInputParameters.seedlist.seedmethod == "random"){
		cout << "number of seed: " << modelInputParameters.seedlist.nseed <<endl;
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		cout<< "duration of the high risk period: " << modelInputParameters.seedlist.hrp <<endl;
	}
	
	//keep information for the health board if interest
	vector<double> pf_byage = observations.pf_pop[modelInputParameters.herd_id - 1];//define frailty structure of the shb of interest.

	//create vector of fixed parameters		
	vector<params> fixed_parameters(observations.waifw_norm.size());	
	for (unsigned int var = 0; var < fixed_parameters.size(); ++var) {
		fixed_parameters[var] = modelInputParameters.paramlist;
	}
	fixed_parameters[0].p_s = fixed_parameters[0].juvp_s;

	//Separate case information for each herd_id
	vector<int> obsHosp_tmp, obsDeaths_tmp;
	int Npop=0;
	int t_index=-1;
	modelInputParameters.seedlist.day_intro=0;
	
	int duration = 0;
	
	if(modelInputParameters.seedlist.seedmethod == "background"){
		select_obs(	Npop, t_index, duration, modelInputParameters.seedlist.day_intro, modelInputParameters.day_shut, 
			obsHosp_tmp, obsDeaths_tmp, observations.cases, observations.deaths, modelInputParameters.herd_id, 
			modelInputParameters.seedlist.hrp);
	} else{
		select_obs(	Npop, t_index, duration, modelInputParameters.seedlist.day_intro, modelInputParameters.day_shut, 
			obsHosp_tmp, obsDeaths_tmp, observations.cases, observations.deaths, modelInputParameters.herd_id, 
			modelInputParameters.paramlist.T_inf + modelInputParameters.paramlist.T_sym);
		if(modelInputParameters.seedlist.seedmethod!= "random"){
			cout << "Warning!! Unknown seeding method - applying _random_ seed method\n";
		}
	}

	obsHosp = obsHosp_tmp;
	obsDeaths = obsDeaths_tmp;
	obsHosp_tmp.clear();	
	obsDeaths_tmp.clear();	
	
	//define age structure and number of hcw of the population at risk
	//compute the number of hcw in the shb
	int N_scot = 0;
	for (unsigned int nn = 0; nn < observations.cases.size()-1; ++nn) {
		N_scot += observations.cases[nn][0];  //compute number of scots in scotland
	}
	double prop_scot = (double)Npop / (double)N_scot;  //proportion of Scots in each shb
	int N_hcw = round(modelInputParameters.totN_hcw * prop_scot); // modulate total number of hcw in Scotland to population in shb

	//adjust population structure to the current population size
	vector<double> agedist = observations.age_pop[modelInputParameters.herd_id - 1];//define age structure of the shb of interest. the -1 is to account for difference in number of rows between two datasets (age_pop does not have a row with column name)
	vector<int> agenums;
	for (unsigned int var = 0; var < agedist.size(); ++var) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(agedist[var] * (Npop - N_hcw))); 

	}		
	//push back the number of hcw in the shb of interest
	agenums.push_back(N_hcw);
	
    cout << "[Health Board settings]:\n";
	cout << "    SHB id: " << modelInputParameters.herd_id <<'\n';
	cout << "    Population size: " << Npop << '\n';
	cout << "    Number of HCW: " << N_hcw << '\n';
	cout << "    Simulation period: " << duration << "days\n";
	cout << "    time step: " << modelInputParameters.tau << "days\n";


	/*---------------------------------------
	 * Model settings
	 *---------------------------------------*/
	//initialise the gsl random number generator with a seed depending of time of the run
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); //�Mersenne Twister� random number generator
	gsl_rng_set(r, time(NULL));

	//initialise the random number generator for importance sampling
    //default_random_engine gen;
	mt19937 gen (time(NULL));

	//declare vectors for priors
	vector<double> flag1, flag2;
	
	flag1 = {
		modelInputParameters.prior_pinf_shape1,
		modelInputParameters.prior_phcw_shape1,
		modelInputParameters.prior_chcw_mean,
		modelInputParameters.prior_d_shape1,
		modelInputParameters.prior_q_shape1, 
		modelInputParameters.prior_rrdh_shape1,
		modelInputParameters.prior_rrdc_shape1,
		modelInputParameters.prior_rrh_shape1,
		modelInputParameters.prior_lambda_shape1
	};	
	flag2 = {
		modelInputParameters.prior_pinf_shape2,
		modelInputParameters.prior_phcw_shape2,
		modelInputParameters.prior_chcw_mean,
		modelInputParameters.prior_d_shape2,
		modelInputParameters.prior_q_shape2,
		modelInputParameters.prior_rrdh_shape2,
		modelInputParameters.prior_rrdc_shape2,
		modelInputParameters.prior_rrh_shape2,
		modelInputParameters.prior_lambda_shape2
	};			

	//declare vectors of outputs/inputs for ABC process
	vector<particle > particleList,particleList1;
	vector<double> weight_val;

	//declare intermediate vectors for ABC smc process
	int pick_val;
	double vlimitKernel[modelInputParameters.nPar];
	double vect_Max[modelInputParameters.nPar];
	double vect_min[modelInputParameters.nPar];

	//initialise the number of accepted particles
	int Nparticle = 0;
	int nsim_val = 0;

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	 cout << "[Simulations]:\n";
	for (int smc = 0; smc < modelInputParameters.nsteps; ++smc) {//todo:

		//the abort statement for keeping the number of particles less than 1000
		bool aborting = false;

		//initialise the counter for the number of simulation prior maximum particle is reached
		int nsim_count = 0;

		//initialise the number of accepted particles
		int counter = 0;

		//initialise the weight and particle lists
		if (smc>0) {
			//update the vectors
			particleList = particleList1;
			particleList1.clear();

			//compute the kernel window
			for (int i = 0; i < modelInputParameters.nPar; ++i) {
				particle valMax1 = *max_element(particleList.begin(),particleList.end(), [&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; } ); //find the element of a vector that is the biggest and return its value
				particle valmin1 = *min_element(particleList.begin(),particleList.end(), [&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; } ); //find the element of a vector  that is the smallest and return its value
				valMax = vect_Max[i] = valMax1.parameter_set[i];
				valmin = vect_min[i] = valmin1.parameter_set[i];
				vlimitKernel[i] = modelInputParameters.kernelFactor*fabs((valMax)-(valmin));
			}
		}


		//create the discrete distribution of the weights for the "importance sampling" process
		for (unsigned int i = 0; i < particleList.size(); ++i) {
			weight_val.push_back(particleList[i].weight);
		}
		//weight_val = weight_val1;
		discrete_distribution<int> weight_distr(weight_val.begin(), weight_val.end());
		weight_val.clear();

	/*---------------------------------------
	 * simulate the infection data set
	 *---------------------------------------*/
		//omp_set_num_threads(num_threads); // Use maximum num_threads threads for all consecutive parallel regions
		//#pragma omp parallel for default(shared), private(pick_val), firstprivate(weight_distr)
		for (int kk = 0; kk < modelInputParameters.nSim; ++kk) {//todo
			//abort statement if number of accepted particles reached nParticLimit particles
			//#pragma omp flush (aborting)
			if (!aborting) {

				//Update progress
				if (counter >= (modelInputParameters.nParticalLimit)) {
					aborting = true;
					//#pragma omp flush (aborting)
				}

				//declare and initialise output variables
				particle outs_vec;
				outs_vec.iter = kk;
				outs_vec.nsse_cases = 0.0;
				outs_vec.nsse_deaths = 0.0;
//				outs_vec.sum_sq = 1.e06;
				for (int i = 0; i < modelInputParameters.nPar; ++i) {
					outs_vec.parameter_set.push_back(0.0);
				}
				//pick the values of each particles
				if (smc==0) {
					//pick randomly and uniformly parameters' value from priors
					parameter_select_first_step(outs_vec.parameter_set,flag1, flag2, r, modelInputParameters.nPar);
				} else {
					//sample 1 particle from the previously accepted particles and given their weight (also named "importance sampling")
				    pick_val = weight_distr(gen);
				    parameter_select_nsteps(outs_vec.parameter_set, modelInputParameters.nPar, r, particleList, pick_val, vlimitKernel,vect_min,vect_Max);
				}
				//run the model and compute the different measures for each potential parameters value
				model_select(smc,outs_vec, fixed_parameters, observations.cfr_byage,pf_byage,
							observations.waifw_norm, observations.waifw_sdist, observations.waifw_home,
							agenums, modelInputParameters.tau, duration, modelInputParameters.seedlist,
							modelInputParameters.day_shut, Npop, r, obsHosp, obsDeaths);
				//count the number of simulations that were used to reach the maximum number of accepted particles
				//#pragma omp critical
					{
					if (counter < modelInputParameters.nParticalLimit) ++nsim_count;
					}
				//if the particle agrees with the different criteria defined for each ABC-smc step
				//if(counter < nParticLimit && outs_vec.nsse_cases <= toleranceLimit[smc]){
				//if(counter < nParticLimit && outs_vec.nsse_deaths <= toleranceLimit[smc]){
				
				if (counter < modelInputParameters.nParticalLimit && outs_vec.nsse_cases <= modelInputParameters.toleranceLimit[smc]) {	
					// Apply the death tolerance limit if applicable
					if (!modelInputParameters.apply_death_tolerance_limit || (modelInputParameters.apply_death_tolerance_limit && outs_vec.nsse_deaths <= modelInputParameters.toleranceLimit[smc]*1.5)) {
									
						//#pragma omp critical
						{
							weight_calc(smc,Nparticle, particleList, outs_vec, vlimitKernel, modelInputParameters.nPar);
							particleList1.push_back(outs_vec);
							++counter;		//count the number of accepted particles
							if(counter % 10 == 0) cout << "|" << flush;
							//cout << counter << " " ;
						}
					}
				}			
			}
		}


		Nparticle = counter;
		nsim_val = nsim_count;

		//time taken per step
		double time_taken;
		if(smc == 0){
			time_taken = double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		} else {
			time_taken = double( clock() - time_taken1 ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		}

		/*---------------------------------------
		 * Outputs
		 *---------------------------------------*/
		// Output on screen of the number of accepted particles, 
		// the number of simulations and the computation time at each step
		cout << "\nStep:" << smc
			<< ", <number of accepted particles> " << Nparticle
			<< "; <number of simulations> " << nsim_val
			<< "; <computation time> " <<  time_taken
			<< " seconds.\n";

		//break the ABC-smc at the step where no particles were accepted
		if (Nparticle > 0) {
			WriteOutputsToFiles(smc, modelInputParameters.herd_id, Nparticle, modelInputParameters.nPar, particleList1);
		}
	}

	//output on screen the overall computation time
	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
}
