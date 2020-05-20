#include "Model.h"
#include "IO.h"
#include "ModelTypes.h"
#include "FittingProcess.h"
#include "Observations.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <numeric>
#include <random>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <functional>

namespace EERAModel {
namespace Model {

static void flow(gsl_rng * r, unsigned int& pop_from, unsigned& pop_to, double rate, unsigned int& outs);

static void infspread(gsl_rng * r, std::vector<int>& pop, int& deaths, int& deathsH, int& detected, 
				::EERAModel::params fixed_parameters, std::vector<double> parameter_set, std::vector<double> cfr_tab,
				double pf_val, double lambda);

static std::vector<double> Lambda(std::vector<double> parameter_set,std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist,std::vector<std::vector<double>> waifw_home, 
				std::vector<std::vector<int>> pops, int shut);

static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage,
				std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist, std::vector<std::vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, std::vector<int> agenums, 
				double tau, gsl_rng * r, std::vector<int> &sim_status,std::vector<std::vector<int>> &ends,
				std::vector<int> &death_status,std::vector<int> &deathH_status);

/**
 * @brief Get the population of a region
 * 
 * Within the input cases data, the total population of the region is located in the first column
 * of the row corresponding to the region
 * 
 * @param obs Input observations
 * @param region_id Index of the region within the observation data set
 * 
 * @return The population of the region
 */
static inline int GetPopulationOfRegion(const InputObservations& obs, int region_id);

/**
 * @brief Compute the agenums 
 * 
 * @param Npop Population
 * @param Nhcw Number of health care workers
 * @param obs Model observations
 * 
 * @return agenums
 */
static std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const InputObservations& obs);

/**
 * @brief Compute the Kernel Window
 * 
 * @param nPar Number of parameters to compute kernel for
 * @param particleList List of previously accepted particles
 * @param kernelFactor Common kernel multiplicative factor
 * @param vlimitKernel Storage for the computed kernel
 * @param vect_Max Storage for maximum values of ranges
 * @param vect_Min Storage for minimum values of ranges
 */
static void ComputeKernelWindow(int nPar, const std::vector<particle>& particleList,
	double kernelFactor, std::vector<double>& vlimitKernel, std::vector<double>& vect_Max, 
	std::vector<double>& vect_Min);

/**
 * @brief Compute weight distribution
 * 
 * Create the discrete distribution of the weights for the "importance sampling" process
 * 
 * @param particleList List of previously accepted particles
 * 
 * @return Weight distribution
 */
static std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList);

void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::InputObservations observations,
		 gsl_rng* r,
		 std::mt19937& gen,
		 const std::string& outDirPath,
		 EERAModel::Utilities::logging_stream::Sptr log) {
	/*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;

    (*log) << "[Settings]:\n";
	(*log)<< "number of parameters tested: "<< modelInputParameters.nPar << std::endl;
    (*log)<< "seeding method: "<< modelInputParameters.seedlist.seedmethod<<  std::endl;
	if(modelInputParameters.seedlist.seedmethod == "random"){
		(*log) << "number of seed: " << modelInputParameters.seedlist.nseed << std::endl;
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		(*log)<< "duration of the high risk period: " << modelInputParameters.seedlist.hrp << std::endl;
	}

    (*log) << "[Parameter Settings]:\n";
    (*log)<< "parameter values: ";
	(*log)<< "  latent period (theta_l): " << modelInputParameters.paramlist.T_lat <<std::endl;
	(*log)<< "  pre-clinical period (theta_i): " << modelInputParameters.paramlist.T_inf <<std::endl;
	(*log)<< "  asymptomatic period (theta_r): " << modelInputParameters.paramlist.T_rec <<std::endl;
	(*log)<< "  symptomatic period (theta_s): " << modelInputParameters.paramlist.T_sym <<std::endl;
	(*log)<< "  hospitalisation stay (theta_h): " << modelInputParameters.paramlist.T_hos <<std::endl;
	(*log)<< "  pre-adult probability of symptoms devt (p_s[0]): " << modelInputParameters.paramlist.juvp_s <<std::endl;
	
	//keep information for the health board if interest
	std::vector<double> pf_byage = observations.pf_pop[modelInputParameters.herd_id - 1];//define frailty structure of the shb of interest.

	//create vector of fixed parameters		
	std::vector<params> fixed_parameters(observations.waifw_norm.size());	
	for (unsigned int var = 0; var < fixed_parameters.size(); ++var) {
		fixed_parameters[var] = modelInputParameters.paramlist;
	}

	//Separate case information for each herd_id
	std::vector<int> obsHosp_tmp, obsDeaths_tmp;
	modelInputParameters.seedlist.day_intro=0;
	
	int duration = 0;
	
	int time_back;
	if(modelInputParameters.seedlist.seedmethod == "background") {
		time_back = modelInputParameters.seedlist.hrp;
	} else {
		time_back = modelInputParameters.paramlist.T_inf + modelInputParameters.paramlist.T_sym;
		
		if (modelInputParameters.seedlist.seedmethod != "random") {
			(*log) << "Warning!! Unknown seeding method - applying _random_ seed method\n";
		}
	}

	int population = GetPopulationOfRegion(observations, modelInputParameters.herd_id);
	
	const std::vector<int>& regionalCases = observations.cases[modelInputParameters.herd_id];
	const std::vector<int>& regionalDeaths = observations.deaths[modelInputParameters.herd_id];
	const std::vector<int>& timeStamps = observations.cases[0];

	Observations::select_obs(duration, modelInputParameters.seedlist.day_intro, 
		modelInputParameters.day_shut, obsHosp_tmp, obsDeaths_tmp, timeStamps, regionalCases,
		regionalDeaths, time_back);

	(*log) << "Number of days of obs cases: " << obsHosp_tmp.size() << std::endl;
	(*log) << "Number of days of obs deaths: " << obsDeaths_tmp.size() << std::endl;

	std::vector<int> obsHosp(obsHosp_tmp);
	std::vector<int> obsDeaths(obsDeaths_tmp);
	
	//define age structure and number of hcw of the population at risk
	//compute the number of hcw in the shb
	int N_scot = 0;
	for (unsigned int nn = 0; nn < observations.cases.size()-1; ++nn) {
		N_scot += observations.cases[nn][0];  //compute number of scots in scotland
	}
	double prop_scot = (double)population / (double)N_scot;  //proportion of Scots in each shb
	int N_hcw = round(modelInputParameters.totN_hcw * prop_scot); // modulate total number of hcw in Scotland to population in shb

	std::vector<int> agenums = ComputeAgeNums(modelInputParameters.herd_id, population, N_hcw, observations);
	
    (*log) << "[Health Board settings]:\n";
	(*log) << "    SHB id: " << modelInputParameters.herd_id <<'\n';
	(*log) << "    Population size: " << population << '\n';
	(*log) << "    Number of HCW: " << N_hcw << '\n';
	(*log) << "    Simulation period: " << duration << "days\n";
	(*log) << "    time step: " << modelInputParameters.tau << "days\n";
	
	const std::vector<double> flag1 = {
		modelInputParameters.prior_pinf_shape1,
	 	modelInputParameters.prior_phcw_shape1,
		modelInputParameters.prior_chcw_mean, 
		modelInputParameters.prior_d_shape1, 
		modelInputParameters.prior_q_shape1,
		modelInputParameters.prior_ps_shape1,
		modelInputParameters.prior_phf_shape1, 
		modelInputParameters.prior_rrdh_shape1,
		modelInputParameters.prior_lambda_shape1
	};	
	const std::vector<double> flag2 = {
		modelInputParameters.prior_pinf_shape2, 
		modelInputParameters.prior_phcw_shape2, 
		modelInputParameters.prior_chcw_mean, 
		modelInputParameters.prior_d_shape2, 
		modelInputParameters.prior_q_shape2,
		modelInputParameters.prior_ps_shape2,
		modelInputParameters.prior_phf_shape2,
		modelInputParameters.prior_rrdh_shape2,
		modelInputParameters.prior_lambda_shape2
	};

	//declare vectors of outputs/inputs for ABC process
	std::vector<particle > particleList, particleList1;

	//initialise the number of accepted particles
	int prevAcceptedParticleCount = 0;

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	(*log) << "[Simulations]:\n";
	for (int smc = 0; smc < modelInputParameters.nsteps; ++smc) {//todo:

		//the abort statement for keeping the number of particles less than 1000
		bool aborting = false;

		//initialise the counter for the number of simulation prior maximum particle is reached
		int nsim_count = 0;

		//initialise the number of accepted particles
		int acceptedParticleCount = 0;

		//initialise the weight and particle lists
		std::vector<double> vect_Max(modelInputParameters.nPar, 0.0);
		std::vector<double> vect_min(modelInputParameters.nPar, 0.0);
		std::vector<double> vlimitKernel(modelInputParameters.nPar, 0.0);

		if (smc > 0) {
			//update the vectors
			particleList = particleList1;
			particleList1.clear();

			ComputeKernelWindow(modelInputParameters.nPar, particleList, 
				modelInputParameters.kernelFactor, vlimitKernel, vect_Max, vect_min);
		}

		std::discrete_distribution<int> weight_distr = ComputeWeightDistribution(particleList);

	/*---------------------------------------
	 * simulate the infection data set
	 *---------------------------------------*/
		//omp_set_num_threads(num_threads); // Use maximum num_threads threads for all consecutive parallel regions
		//#pragma omp parallel for default(shared), private(pick_val), firstprivate(weight_distr)
		for (int sim = 0; sim < modelInputParameters.nSim; ++sim) {//todo
			//abort statement if number of accepted particles reached nParticLimit particles
			//#pragma omp flush (aborting)
			if (!aborting) {

				//Update progress
				if (acceptedParticleCount >= modelInputParameters.nParticalLimit) {
					aborting = true;
					//#pragma omp flush (aborting)
				}

				//declare and initialise output variables
				particle outs_vec;
				outs_vec.iter = sim;
				outs_vec.nsse_cases = 0.0;
				outs_vec.nsse_deaths = 0.0;
//				outs_vec.sum_sq = 1.e06;
				for (int i = 0; i < modelInputParameters.nPar; ++i) {
					outs_vec.parameter_set.push_back(0.0);
				}
				//pick the values of each particles
				if (smc==0) {
					//pick randomly and uniformly parameters' value from priors
					outs_vec.parameter_set = FittingProcess::parameter_select_initial(flag1, flag2,
						r, modelInputParameters.nPar);
					
				} else {
					//sample 1 particle from the previously accepted particles and given their weight (also named "importance sampling")
				    int pick_val = weight_distr(gen);
				    outs_vec.parameter_set = FittingProcess::parameter_select(
						modelInputParameters.nPar, r, particleList, pick_val, vlimitKernel,vect_min,vect_Max);
				}

				//run the model and compute the different measures for each potential parameters value
				model_select(outs_vec, fixed_parameters, observations.cfr_byage, pf_byage,
							observations.waifw_norm, observations.waifw_sdist, observations.waifw_home,
							agenums, modelInputParameters.tau, duration, modelInputParameters.seedlist,
							modelInputParameters.day_shut, population, r, obsHosp, obsDeaths);

				//count the number of simulations that were used to reach the maximum number of accepted particles
				//#pragma omp critical
					{
						if (acceptedParticleCount < modelInputParameters.nParticalLimit) ++nsim_count;
					}
				//if the particle agrees with the different criteria defined for each ABC-smc step
				if (
					acceptedParticleCount < modelInputParameters.nParticalLimit &&
					outs_vec.nsse_cases <= modelInputParameters.toleranceLimit[smc] &&
					outs_vec.nsse_deaths <= modelInputParameters.toleranceLimit[smc]*1.5
					) {				
						//#pragma omp critical
						{
							FittingProcess::weight_calc(smc, prevAcceptedParticleCount, particleList, outs_vec, 
								vlimitKernel, modelInputParameters.nPar);
							particleList1.push_back(outs_vec);
							++acceptedParticleCount;
							if (acceptedParticleCount % 10 == 0) (*log) << "|" << std::flush;
							//(*log) << acceptedParticleCount << " " ;
						}
					
				}			
			}
		}

		prevAcceptedParticleCount = acceptedParticleCount;

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
		(*log) << "\nStep:" << smc
			<< ", <number of accepted particles> " << prevAcceptedParticleCount
			<< "; <number of simulations> " << nsim_count
			<< "; <computation time> " <<  time_taken
			<< " seconds.\n";

		//break the ABC-smc at the step where no particles were accepted
		if (prevAcceptedParticleCount > 0) {
			IO::WriteOutputsToFiles(smc, modelInputParameters.herd_id, prevAcceptedParticleCount,
				modelInputParameters.nPar, particleList1, outDirPath);
		}
	}

	//output on screen the overall computation time
	(*log) << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
}

void model_select(EERAModel::particle& outvec, const std::vector<params>& fixed_parameters,
	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, int population, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths) {

	//---------------------------------------
	// the root model
	//---------------------------------------
	std::vector<int> sim_status;
	/**@todo death_status is unused anywhere (although it is populated in my_model) */
	std::vector<int> death_status;
	std::vector<int> deathH_status;
	std::vector<std::vector<int>> ends;
	
	my_model(outvec.parameter_set, fixed_parameters, cfr_byage, pf_byage,waifw_norm, waifw_sdist,
		waifw_home,	duration, seedlist, day_shut, population,agenums, tau, r, sim_status, ends,
		death_status, deathH_status);

	//---------------------------------------
	// Compute the  sum of squared errors
	//---------------------------------------
	double sum_sq_cases = ::EERAModel::Utilities::sse_calc<int>(sim_status, obsHosp);
	double sum_sq_deaths= ::EERAModel::Utilities::sse_calc<int>(deathH_status, obsDeaths);
	
	//---------------------------------------
	// Compute the deviation of the fit expressed as % total number of cases
	//---------------------------------------
	//total numbers of cases and deaths (at hospital)
	int sum_obs_cases = std::accumulate(obsHosp.begin(), obsHosp.end(), 0);
	int sum_obs_deaths = std::accumulate(obsDeaths.begin(), obsDeaths.end(), 0);
	
	outvec.nsse_cases = sqrt(sum_sq_cases)/ (double) sum_obs_cases;
	outvec.nsse_deaths = sqrt(sum_sq_deaths)/ (double) sum_obs_deaths;

	//---------------------------------------
	// Return all selection measures
	//---------------------------------------
	outvec.simu_outs = sim_status;
	outvec.death_outs = deathH_status;
	outvec.end_comps = ends;
}
 
static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage,
				std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist, std::vector<std::vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, std::vector<int> agenums, 
				double tau, gsl_rng * r, std::vector<int> &sim_status,std::vector<std::vector<int>> &ends,
				std::vector<int> &death_status,std::vector<int> &deathH_status) {

	std::vector<std::vector<double>> parameter_fit(waifw_norm.size(), parameter_set);	
    parameter_fit[0][5] = fixed_parameters[0].juvp_s;

	//sets up an array for the population at each timestep in each age and disease category	
	//also set up the age distribution of old ages as target for disease introduction
	int n_agegroup = waifw_norm.size();
	int n_comparts = 18;
	std::vector<std::vector<int>> poparray(n_agegroup, std::vector<int>(n_comparts));	  
	for ( int age = 0; age < (n_agegroup); ++age) {
		for ( int st = 0; st < n_comparts; ++st) {
			if(st == 0) {
				poparray[age][st] = agenums[age]; //set up the starting population as fully susceptible
			} else {
				poparray[age][st] = 0;
			}
		}
	}

	//introduce disease at t=0. if seedmethod != "background"
	//WARNING!! MINOR BUG HERE: gsl multinomial sometimes freeze for unknown reasons. replaced by random flat, picking a value where to seed infection. will only works because seed=1.
	if(seedlist.seedmethod != "background"){
		std::vector<int> startdist = {0,0,0,0,0,0};
		startdist[(int) gsl_ran_flat(r, 0, 6)] = seedlist.nseed;

		for (int age = 1; age < (n_agegroup - 1); ++age) {
			poparray[age][0] -=  startdist[age - 1]; //take diseased out of S
			poparray[age][3] +=  startdist[age - 1]; //put diseased in I	
		}	
	}

  	//initialize saving of the detection for t=0
	sim_status.push_back(0); 
	death_status.push_back(0); 
	deathH_status.push_back(0); 

	//run the simulation
	for (int tt = 1; tt < (int)ceil(duration/tau); ++tt) {	
		//identify the lock down
		int inLockdown = tt > day_shut ? 1 : 0;
		
		//introduce disease from background infection until lockdown
	    if(!inLockdown) {
			if(seedlist.seedmethod == "background"){
				double bkg_lambda = parameter_set[parameter_set.size()-1];
				unsigned int startdz = gsl_ran_poisson(r, (double)Npop * bkg_lambda);	//how many diseased is introduced in each given day before lockdown
		
				std::vector<int> startdist = {0,0,0,0,0,0};
				int pickcomp = (int) gsl_ran_flat(r, 0, 6);
				while(poparray[pickcomp][0] < 1) {
					pickcomp = (int) gsl_ran_flat(r, 0, 6);
				} 
				startdist[pickcomp] = startdz;

				for (int age = 1; age < (n_agegroup - 1); ++age) {
					int nseed = std::min(poparray[age][0], startdist[age - 1]);
					poparray[age][0] -=  nseed; //take diseased out of S
					poparray[age][1] +=  nseed; //put diseased in E	
				}	
			}
	    }

		// Compute the forces of infection
		std::vector<double> lambda = Lambda(parameter_set, waifw_norm, waifw_sdist,	waifw_home,
			poparray, inLockdown);	
	
		// Step each agegroup through infections
		int deaths = 0;
		int deathsH = 0;
		int detected = 0;
		for (int age = 0; age < n_agegroup; ++age) {	
			infspread(r, poparray[age], deaths, deathsH, detected,fixed_parameters[age],parameter_fit[age],
				cfr_byage[age], pf_byage[age], lambda[age]);
		}

		sim_status.push_back(detected); 
		death_status.push_back(deaths); 
		deathH_status.push_back(deathsH);
    }	
	
	// save the population in each epi compt for the last day
	ends = poparray;
}

static void infspread(gsl_rng * r, std::vector<int>& pop, int& deaths, int& deathsH, int& detected, 
	params fixed_parameters, std::vector<double> parameter_set, std::vector<double> cfr_tab, double pf_val, double lambda){
		
    unsigned int S=pop[0], E=pop[1], E_t=pop[2], I_p=pop[3], I_t=pop[4],I1=pop[5],I2=pop[6],I3=pop[7],I4=pop[8];
	unsigned int I_s1=pop[9], I_s2=pop[10], I_s3=pop[11], I_s4=pop[12], I_fs=pop[13];
	unsigned int H=pop[14], H_f=pop[15], R=pop[16], D=pop[17] ;
	
	double T_lat= fixed_parameters.T_lat;
//	double p_s= fixed_parameters.p_s;
	double T_inf= fixed_parameters.T_inf;
	double T_rec= fixed_parameters.T_rec;
	double T_sym= fixed_parameters.T_sym;
	double T_hos= fixed_parameters.T_hos;
	
	double p_h= cfr_tab[0];
//	double cfr= cfr_tab[1];	
	double p_d= cfr_tab[2];	
	//double p_dc= cfr_tab[3];

	double p_s= parameter_set[5];
//	double T_hos= parameter_set[6];
	double p_hf = parameter_set[6];
	double rrdh = parameter_set[7];	

    // hospitalized  - non-frail
    unsigned int newdeathsH=0;
//	cout << "newdeathH: " << p_d * (1 / T_hos) << "\n";
	flow(r, H, D, p_d * (1 / T_hos), newdeathsH);
	
    unsigned int recoverH=0;
//	cout << "recoverH: " << (1 - p_d) * (1 / T_hos) << "\n";
	flow(r, H, R, (1 - p_d) * (1 / T_hos), recoverH);
	
	//hospitalize - frail
    unsigned int newdeathsH_f=0;
//	cout << "newdeathsH_f: " << rrdh * p_d * (1 / T_hos) << "\n";
	flow(r, H_f, D, rrdh * p_d * (1 / T_hos), newdeathsH_f);
	
    unsigned int recoverH_f=0;
//	cout << "recoverH_f: " << (1 - rrdh * p_d) * (1 / T_hos) << "\n";
	flow(r, H_f, R, (1 - rrdh * p_d) * (1 / T_hos), recoverH_f);
	
    // symptomatic - non-frail
    unsigned int hospitalize=0;
//	cout << "hospitalize: " << p_h * (4 / T_sym) << "\n";
	flow(r, I_s4, H, p_h * (4 / T_sym), hospitalize);

    unsigned int recoverI_s=0;
//	cout << "recoverI_s: " << (1 - p_h) * ( 4 / T_sym) << "\n";
	flow(r, I_s4, R, (1 - p_h) * ( 4 / T_sym), recoverI_s);
	
    unsigned int Is_from3_to_4=0;
//	cout << "Is_from3_to_4: " << (4 / T_sym) << "\n";
	flow(r, I_s3, I_s4, (4 / T_sym), Is_from3_to_4);
	
    unsigned int Is_from2_to_3=0;
//	cout << "Is_from2_to_3: " << (4 / T_sym) << "\n";
	flow(r, I_s2, I_s3, (4 / T_sym), Is_from2_to_3);
		
    unsigned int Is_from1_to_2=0;
//	cout << "Is_from1_to_2: " << (4 / T_sym) << "\n";
	flow(r, I_s1, I_s2, (4 / T_sym), Is_from1_to_2);

	// symptomatic -frail
    unsigned int hospitalize_f=0;
//	cout << "hospitalize_f: " << p_hf * (2 / T_sym) << "\n";
	flow(r, I_fs, H_f, p_hf * (1 / T_sym), hospitalize_f);	
	
    unsigned int newdeaths_f=0;
//	cout << "newdeaths_f: " << ( 1 - p_hf ) * ( 2 / T_sym) << "\n";
	flow(r, I_fs, D, ( 1 - p_hf ) * ( 1 / T_sym), newdeaths_f);
	
	// asymptomatic
    unsigned int recoverI=0;
//	cout << "recoverI: " << ( 4/T_rec ) << "\n";
	flow(r, I4, R, ( 4/T_rec ), recoverI);
	
    unsigned int I_from3_to_4=0;
//	cout << "I_from3_to_4: " << ( 4/T_rec ) << "\n";
	flow(r, I3, I4, ( 4/T_rec ), I_from3_to_4);
	
    unsigned int I_from2_to_3=0;
//	cout << "I_from2_to_3: " << ( 4/T_rec ) << "\n";
	flow(r, I2, I3, ( 4/T_rec ), I_from2_to_3);
 
    unsigned int I_from1_to_2=0;
//	cout << "I_from1_to_2: " << ( 4/T_rec ) << "\n";
	flow(r, I1, I2, ( 4/T_rec ), I_from1_to_2);
	
	// infectious - pre-clinical
    unsigned int newasymptomatic=0;
//	cout << "newasymptomatic: " << p_d * (1 / T_hos) << "\n";
	flow(r, I_p, I1, p_d * (1 / T_hos), newasymptomatic);	

    unsigned int newsymptomatic=0;
//	cout << "newsymptomatic: " << (1 - pf_val) * p_s * ( 1/T_inf ) << "\n";
	flow(r, I_p, I_s1, (1 - pf_val) * p_s * ( 1/T_inf ), newsymptomatic);	
	
    unsigned int newinffrail=0;
//	cout << "newinffrail: " << pf_val * ( 1/T_inf ) << "\n";
	flow(r, I_p, I_fs, pf_val * ( 1/T_inf ), newinffrail);		
	
    // latent
    unsigned int infectious=0;
//	cout << "infectious: " << 1/T_lat << "\n";
	flow(r, E, I_p, ( 1/T_lat ), infectious);
	
    unsigned int infectious_t=0;
//	cout << "infectious_t: " << (1 / T_lat) << "\n";
	flow(r, E_t, I_t, ( 1/T_lat ), infectious_t);
	
    // susceptible
    unsigned int newinfection=0;
//	cout << "newinfection: " << lambda << "\n";
	flow(r, S, E, lambda, newinfection);
	
    // recover population
	std::vector<int> newpop= {(int)S,(int)E,(int)E_t,(int)I_p,(int)I_t, (int)I1,(int)I2,(int)I3,(int)I4,(int)I_s1,(int)I_s2,(int)I_s3,(int)I_s4,(int)I_fs,(int)H,(int)H_f,(int)R,(int)D};
	
	pop = newpop;
	deaths += (newdeathsH + newdeathsH_f + newdeaths_f) ;
	deathsH += newdeathsH + newdeathsH_f;
	detected += (hospitalize + hospitalize_f + infectious_t);	
}

static std::vector<double> Lambda(std::vector<double> parameter_set,std::vector<std::vector<double>> waifw_norm,
			std::vector<std::vector<double>> waifw_sdist,std::vector<std::vector<double>> waifw_home, 
			std::vector<std::vector<int>> pops, int shut) {
  
	int n_agegroup = waifw_norm.size();

	double p_i = parameter_set[0];	
	double p_hcw = parameter_set[1];	
	double c_hcw = parameter_set[2];	
	double d_val = parameter_set[3];					
	double q_val = parameter_set[4];					

	double quarantined = 0.0; //mean proportion reduction in contacts due to quarantine
	
	//contact matrix
	std::vector<std::vector<double>> waifw(n_agegroup, std::vector<double> (n_agegroup)); 

	//age-dependent transmission rate
	std::vector<std::vector<double>> beta(n_agegroup, std::vector<double> (n_agegroup)); 

	for (int from = 0; from < n_agegroup; ++from) {
		for (int to = 0; to < n_agegroup; ++to) {
			
			if(shut == 0){
				// contact network during normal period
				waifw[from][to]  = waifw_norm[from][to] ;
			} else {
				//contact network during shutdown period, assuming a proportion d will properly do it
				waifw[from][to] = (1.0 - d_val) * waifw_sdist[from][to] + (d_val) * waifw_home[from][to];
			}
			
			beta[from][to] = waifw[from][to] * p_i;
			if(waifw[from][to] > 0) {
				quarantined += waifw_home[from][to] / waifw[from][to];
			}
		}
	} 
  
	// mean proportion reduction in contacts due to quarantine accounting for compliance
	quarantined = 1 - (quarantined / (double)(n_agegroup * n_agegroup) ) * (1 - q_val);
	
	// compute the pressure from infectious groups, normalizes by group size
	// compute the pressure from hospitalised
	int inf_hosp = 0;
	std::vector<double> I_mat(n_agegroup);
	
	for (int from = 0; from < n_agegroup; ++from) {
		if (from < (n_agegroup - 1)) { //n-agegroup-1 to not account for hcw (as different contact)
			// are considered those that are infectious pre-clinical, asymp, tested or symptomatic only
			// asym are infectious pre-clinical and asymptomatic
			// sym are infectious frail and not frail plus those that are tested positive
			int tot_asym = pops[from][3] + pops[from][5] + pops[from][6] + pops[from][7] + pops[from][8];
			int tot_sym = pops[from][4] + pops[from][9] + pops[from][10] + pops[from][11] + pops[from][12] + pops[from][13];
			
			I_mat[from] = (double)tot_asym + (1-quarantined) * (double)tot_sym;
			
			//normalisation shouldnt account for dead individuals	
			int tot_pop = accumulate(pops[from].begin(), pops[from].end(), 0);
			I_mat[from] = (double)I_mat[from] / (double)(tot_pop); 		
		}

		//sum up of number of hospitalised cases (frail and non-frail)
		inf_hosp += pops[from][14] + pops[from][15];
	}

	std::vector<double> lambda(n_agegroup);
	
	// Sum up infectious pressure from each age group
	for (int to = 0; to < n_agegroup; ++to) {
		// Assume perfect self isolation and regular testing of infected HCW if infected
		for (int from = 0; from < (n_agegroup-1); ++from) {
			lambda[to] += (beta[to][from] * I_mat[from]);
		}
	}

	int tot_pop = accumulate(pops[n_agegroup - 1].begin(), pops[n_agegroup - 1].end(), 0);
	lambda[n_agegroup-1] = p_hcw * c_hcw * ((double)inf_hosp/(double)(tot_pop) );

	return lambda;
}

void flow(gsl_rng * r, unsigned int& pop_from, unsigned int& pop_to, double rate, unsigned int& outs){
    outs = gsl_ran_poisson(r, rate * (double)pop_from); //symptomatic become hospitalized
    outs = std::min(pop_from, outs);
	pop_to += outs;
	pop_from -= outs;
}

static std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const InputObservations& obs) {
	std::vector<int> agenums;
	
	// define age structure of the shb of interest. the -1 is to account for difference in number of
	// rows between two datasets (age_pop does not have a row with column name)
	const auto& agedist = obs.age_pop[shb_id - 1];
	
	for (const auto& var : agedist) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(var * (Npop - N_hcw))); 

	}		
	agenums.push_back(N_hcw);

	return agenums;
}

static void ComputeKernelWindow(int nPar, const std::vector<particle>& particleList, double kernelFactor,
	std::vector<double>& vlimitKernel, std::vector<double>& vect_Max, std::vector<double>& vect_Min) {

	//compute the kernel window
	for (int i = 0; i < nPar; ++i) {
		
		std::function<bool(particle, particle)> compare = 
			[&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; };

		particle valMax1 = *std::max_element(particleList.begin(), particleList.end(), compare);
		
		particle valmin1 = *std::min_element(particleList.begin(), particleList.end(), compare);

		vect_Max[i] = valMax1.parameter_set[i];
		vect_Min[i] = valmin1.parameter_set[i];
		
		vlimitKernel[i] = kernelFactor * fabs(vect_Max[i] - vect_Min[i]);
	}	
}

static std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList) {
	
	std::vector<double> weight_val;
	for (auto p : particleList) {
		weight_val.push_back(p.weight);
	}
	
	return std::discrete_distribution<int>(weight_val.begin(), weight_val.end());
}

static inline int GetPopulationOfRegion(const InputObservations& obs, int region_id)
{
	return obs.cases[region_id][0];
}

} // namespace Model
} // namespace EERAModel
