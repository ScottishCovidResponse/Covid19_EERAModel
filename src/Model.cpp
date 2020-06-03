#include "Model.h"
#include "IO.h"
#include "ModelTypes.h"
#include "FittingProcess.h"
#include "Observations.h"
#include "InferenceParameters.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_gamma.h>

#include <functional>

namespace EERAModel {
namespace Model {

static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
				int n_sim_steps, Random::RNGInterface::Sptr rng, Status& status);

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

static void weekly(std::vector<int>& reduced, const std::vector<int>& original);

void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::InputObservations observations,
		 Random::RNGInterface::Sptr rng,
		 const std::string& outDirPath,
		 EERAModel::Utilities::logging_stream::Sptr log) {
	/*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;

    (*log) << "[Settings]:\n";
	(*log) << "    number of parameters tested: "<< modelInputParameters.nPar << std::endl;
    (*log) << "    seeding method: "<< modelInputParameters.seedlist.seedmethod<<  std::endl;
	if (modelInputParameters.seedlist.seedmethod == "random"){
		(*log) << "    number of seed: " << modelInputParameters.seedlist.nseed << std::endl;
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		(*log) << "    duration of the high risk period (hrp): " << modelInputParameters.seedlist.hrp << std::endl;
	}

    (*log) << "[Fixed parameter values]:\n";
	(*log) << "    latent period (theta_l): " << modelInputParameters.paramlist.T_lat <<std::endl;
	(*log) << "    pre-clinical period (theta_i): " << modelInputParameters.paramlist.T_inf <<std::endl;
	(*log) << "    asymptomatic period (theta_r): " << modelInputParameters.paramlist.T_rec <<std::endl;
	(*log) << "    symptomatic period (theta_s): " << modelInputParameters.paramlist.T_sym <<std::endl;
	(*log) << "    hospitalisation stay (theta_h): " << modelInputParameters.paramlist.T_hos <<std::endl;
	(*log) << "    pre-adult probability of symptoms devt (p_s[0]): " << modelInputParameters.paramlist.juvp_s <<std::endl;
	(*log) << "    bed capacity at hospital (K): " << modelInputParameters.paramlist.K <<std::endl;
	(*log) << "    relative infectiousness of asymptomatic (u): " << modelInputParameters.paramlist.inf_asym <<std::endl;
	
	//keep information for the health board if interest
	std::vector<double> pf_byage = observations.pf_pop[modelInputParameters.herd_id - 1];//define frailty structure of the shb of interest.

	//create vector of fixed parameters		
	std::vector<params> fixed_parameters(observations.waifw_norm.size());	
	for (unsigned int var = 0; var < fixed_parameters.size(); ++var) {
		fixed_parameters[var] = modelInputParameters.paramlist;
	}

	//Separate case information for each herd_id
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

	std::vector<int> obsHosp, obsDeaths;
	Observations::SelectObservations(duration, modelInputParameters.seedlist.day_intro, 
		modelInputParameters.day_shut, obsHosp, obsDeaths, timeStamps, regionalCases,
		regionalDeaths, time_back, log);

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
		modelInputParameters.prior_rrd_shape1,		
		modelInputParameters.prior_lambda_shape1
	};	
	const std::vector<double> flag2 = {
		modelInputParameters.prior_pinf_shape2, 
		modelInputParameters.prior_phcw_shape2, 
		modelInputParameters.prior_chcw_mean, 
		modelInputParameters.prior_d_shape2, 
		modelInputParameters.prior_q_shape2,
		modelInputParameters.prior_ps_shape2,
		modelInputParameters.prior_rrd_shape2,		
		modelInputParameters.prior_lambda_shape2
	};

    Inference::InferenceParameterGenerator inferenceParameterGenerator(
        modelInputParameters.nPar, rng, flag1, flag2);

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
				for (int i{0}; i < modelInputParameters.nPar; ++i) {
					outs_vec.parameter_set.push_back(0.0);
				}
				//pick the values of each particles
				if (smc==0) {
					//pick randomly and uniformly parameters' value from priors

                    outs_vec.parameter_set = inferenceParameterGenerator.GenerateInitial();
					
				} else {
					// sample 1 particle from the previously accepted particles
                    // and given their weight (also named "importance sampling")
				    int pick_val = weight_distr(rng->MT19937());
				    outs_vec.parameter_set = inferenceParameterGenerator.GenerateWeighted(
                        particleList[pick_val].parameter_set, vlimitKernel, vect_Max, vect_min);
				}

				//run the model and compute the different measures for each potential parameters value
				model_select(outs_vec, fixed_parameters, observations.cfr_byage, pf_byage,
							observations.waifw_norm, observations.waifw_sdist, observations.waifw_home,
							agenums, modelInputParameters.tau, duration, modelInputParameters.seedlist,
							modelInputParameters.day_shut, rng, obsHosp, obsDeaths);

				//count the number of simulations that were used to reach the maximum number of accepted particles
				//#pragma omp critical
					{
						if (acceptedParticleCount < modelInputParameters.nParticalLimit) ++nsim_count;
					}
				//if the particle agrees with the different criteria defined for each ABC-smc step
				if (
					acceptedParticleCount < modelInputParameters.nParticalLimit &&
					outs_vec.nsse_cases <= modelInputParameters.toleranceLimit[smc] &&
					outs_vec.nsse_deaths <= modelInputParameters.toleranceLimit[smc]//*1.5
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
				modelInputParameters.nPar, particleList1, outDirPath, log);
		}
	}

	//output on screen the overall computation time
	(*log) << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
}

void model_select(EERAModel::particle& outvec, const std::vector<params>& fixed_parameters,
	const std::vector<std::vector<double>>& cfr_byage, const std::vector<double>& pf_byage, 
	const std::vector<std::vector<double>>& waifw_norm, const std::vector<std::vector<double>>& waifw_sdist,
	const std::vector<std::vector<double>>& waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, Random::RNGInterface::Sptr rng, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths) {

	//---------------------------------------
	// the root model
	//---------------------------------------

	//initialize saving of the detection for t=0 in simulations, deaths, hospital deaths
	Status status = {{0}, {0}, {0}, {}};
	/**@todo status.deaths is unused anywhere (although it is populated in my_model) */

	const AgeGroupData per_age_data = {waifw_norm, waifw_home, waifw_sdist, cfr_byage, pf_byage};
	const int n_sim_steps = static_cast<int>(ceil(duration/tau));
	
	my_model(outvec.parameter_set, fixed_parameters, per_age_data, seedlist, day_shut,
			agenums, n_sim_steps, rng, status);

	//---------------------------------------
	// compute the  sum of squared errors for daily observations
	//---------------------------------------
	double sum_sq_cases = Utilities::sse_calc<int>(status.simulation, obsHosp);

	//---------------------------------------
	// compute the  sum of squared errors for weekly observations
	//---------------------------------------	
	// aggregate daily values to weekly values 
	// due to reporting process creating uncertaincies within weeks
	std::vector<int> obs_death_red, sim_death_red;
	weekly(obs_death_red, obsDeaths);
	weekly(sim_death_red, status.hospital_deaths);
	
	double sum_sq_deaths = Utilities::sse_calc<int>(sim_death_red, obs_death_red);

	//---------------------------------------
	// Compute the deviation of the fit expressed as % total number of observations
	//---------------------------------------
	//total numbers of cases and deaths (at hospital)
	int sum_obs_cases = std::accumulate(obsHosp.begin(), obsHosp.end(), 0);
    //int sum_obs_deaths = std::accumulate(obsDeaths.begin(), obsDeaths.end(), 0);	
	int sum_obs_deaths = std::accumulate(obs_death_red.begin(), obs_death_red.end(), 0);
		
	// std::cout <<"daily obs: " << sum_obs_deaths<< " , red: " << sum_obs_deaths_red << 
	// " , length red: " <<	obs_death_red.size() << std::endl;
	outvec.nsse_cases = sqrt(sum_sq_cases) / static_cast<double>(sum_obs_cases);
	outvec.nsse_deaths = sqrt(sum_sq_deaths) / static_cast<double>(sum_obs_deaths);

	//---------------------------------------
	// Return all selection measures
	//---------------------------------------
	outvec.simu_outs = status.simulation;
	outvec.death_outs = status.hospital_deaths;
	outvec.end_comps = compartments_to_vector(status.ends);
}
 
std::vector<double> build_population_seed(const std::vector<int>& age_nums)
{
	std::vector<double> _temp = {};
	for ( int age{0}; age < age_nums.size(); ++age)
	{
		// seed only occurs in >20yo and not in HCW
		if( age > 0 && age < age_nums.size()-1 ) 
		{
			_temp.push_back(static_cast<double>(age_nums[age]));
		}
	}

	return _temp;
}

std::vector<Compartments> build_population_array(Random::RNGInterface::Sptr rng, const std::vector<int>& age_nums, const seed& seedlist)
{
	std::vector<Compartments> _temp(age_nums.size(), Compartments());
	for ( int age{0}; age < age_nums.size(); ++age) {
		_temp[age].S = age_nums[age]; //set up the starting population as fully susceptible
	}

	//introduce disease at t=0. if seedmethod != "background"
	//WARNING!! MINOR BUG HERE: gsl multinomial sometimes freeze for unknown reasons. replaced by random flat, picking a value where to seed infection. will only works because seed=1.
	if(seedlist.seedmethod != "background"){
		//	size_t k = sizeof(seed_pop);
		//	unsigned int startdist[k];
		std::vector<int> startdist(6,0);
		startdist[static_cast<int>(rng->Flat(0, 6))] = seedlist.nseed;
		//	gsl_ran_multinomial(r, k, startdz, seed_pop,startdist); //distribute the diseased across the older age categories

		for ( int age{1}; age < age_nums.size()-1; ++age) {
			_temp[age].S -=  startdist[age-1];	// take diseased out of S
			_temp[age].I_p +=  startdist[age-1];	// put diseased in I	
		}	
	}

	return _temp;
}

void generate_diseased(Random::RNGInterface::Sptr rng, std::vector<Compartments>& poparray, std::vector<double>& seedarray, const double& bkg_lambda)
{
	int n_susc = 0;

	for ( int age{1}; age < poparray.size()-1; ++age) 
	{
		seedarray[age-1] = static_cast<double>(poparray[age].S);
		n_susc += poparray[age].S;	
	}

	//how many diseased is introduced in each given day before lockdown
	//as a proportion of number of Susceptible available for background infection (not total population, only 20-70 individuals)
	int startdz = rng->Poisson(static_cast<double>(n_susc) * bkg_lambda);

	//size_t k = sizeof(seed_pop) / sizeof(seed_pop[0]);

	unsigned int startdist[seedarray.size()];

/*	std::vector<int> startdist = {0,0,0,0,0,0};
	int pickcomp = (int)gsl_ran_flat(r, 0, 6);
	while(poparray[pickcomp][0]<1){
		pickcomp = (int)gsl_ran_flat(r, 0, 6);
	} 
	startdist[pickcomp] = startdz;
*/				
	rng->Multinomial(seedarray.size(), startdz, &seedarray[0], startdist); //distribute the diseased across the older age categories

	for ( int age{1}; age < poparray.size()-1; ++age) {
		int nseed = startdist[age-1];
		//std::min(poparray[age][0], startdist[age-1]);
		poparray[age].S -=  nseed;	// take diseased out of S
		poparray[age].E +=  nseed;	// put diseased in E
	}
}

static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
				int n_sim_steps, Random::RNGInterface::Sptr rng, Status& status) {


///	std::cout<< "top0..\n";

	const int n_agegroup = per_age_data.waifw_norm.size();

	// Start without lockdown
	bool inLockdown = false;

	// Assumes that the number of age groups matches the size of the 'agenums' vector
	std::vector<double> seed_pop = build_population_seed(agenums);

	std::vector<std::vector<double>> parameter_fit(per_age_data.waifw_norm.size(), parameter_set);	
	parameter_fit[0][5] = fixed_parameters[0].juvp_s;

//	std::cout<< "top1..\n";
	std::vector<Compartments> poparray = build_population_array(rng, agenums, seedlist);
	

//	std::cout<< "top2..\n";

//	std::cout<< "top3..\n";
	
	//run the simulation
	for (int tt{1}; tt < n_sim_steps; ++tt) {
		//initialize return value
		InfectionState infection_state;	
		
//		std::cout<< "top4..\n";
		//identify the lock down
		if(tt > day_shut){ inLockdown = true; }
		
	
		//introduce disease from background infection until lockdown
		if(!inLockdown && seedlist.seedmethod == "background")
		{
			const double bkg_lambda = parameter_set[parameter_set.size()-1];

			//compute the total number of susceptible and the number of susceptible per age class			
			generate_diseased(rng, poparray, seed_pop, bkg_lambda);
		}

		//compute the forces of infection
		std::vector<double> lambda = generate_lambda_vector(infection_state.hospitalised, parameter_set, fixed_parameters[0].inf_asym, per_age_data,
				poparray, inLockdown);	

		//step each agegroup through infections
		for ( int age{0}; age < (n_agegroup); ++age) {	
			const InfectionState new_spread = infection_spread(rng, poparray[age], infection_state.hospitalised, fixed_parameters[age],parameter_fit[age],
															per_age_data.cfr_byage[age], per_age_data.pf_byage[age],lambda[age]);

			infection_state.deaths += new_spread.deaths;
			infection_state.hospital_deaths += new_spread.hospital_deaths;
			infection_state.detected += new_spread.detected;
		}
//std::cout<< "top9..\n";
//		std::cout << tt << " , " <<  detected << " , " << deaths << " , " << deathsH << '\n';

		status.simulation.push_back(infection_state.detected); 
		status.deaths.push_back(infection_state.deaths); 
		status.hospital_deaths.push_back(infection_state.hospital_deaths);
    }	
	//save the population in each epi compt for the last day
	for ( int age{0}; age < n_agegroup; ++age) {
		status.ends.push_back(poparray[age]);
	}
}

InfectionState infection_spread(Random::RNGInterface::Sptr rng, Compartments& pop, const int& n_hospitalised,
	params fixed_parameters, std::vector<double> parameter_set, std::vector<double> cfr_tab, double pf_val, double lambda)
{
		
/*    unsigned int S=pop[0], E=pop[1], E_t=pop[2], I_p=pop[3], I_t=pop[4],I1=pop[5],I2=pop[6],I3=pop[7],I4=pop[8];
	unsigned int I_s1=pop[9], I_s2=pop[10], I_s3=pop[11], I_s4=pop[12];//, I_fs=pop[13]
	unsigned int H=pop[13], R=pop[14], D=pop[15] ;//, H_f=pop[15]
*/	
	//define the vector of return values
	Compartments newpop(pop);

	//Define infection state for this spread
	InfectionState infection_state;
	
	//define the position of each compartment in pop
	// I_fs = pop.H;
	// H_f  = pop.D;
	
	//fixed parameters
	const double T_lat= fixed_parameters.T_lat;
//	double p_s= fixed_parameters.p_s;
	const double T_inf= fixed_parameters.T_inf;
	const double T_rec= fixed_parameters.T_rec;
	const double T_sym= fixed_parameters.T_sym;
	const double T_hos= fixed_parameters.T_hos;
	const double K=fixed_parameters.K;
	
	double capacity = n_hospitalised/K;
	capacity = std::min(1.0,capacity);
	
	const double p_h= cfr_tab[0];
//	double cfr= cfr_tab[1];	
	const double p_d= cfr_tab[2];	
	//double p_dc= cfr_tab[3];

	const double p_s= parameter_set[5];
	const double rrd= parameter_set[6];
//	double T_hos= parameter_set[6];
//	double p_hf = parameter_set[6];
//	double rrdh = parameter_set[7];	

    // hospitalized  - non-frail
    const int newdeathsH= flow(rng, pop.H, newpop.H, p_d * (1.0 / T_hos));
	newpop.H -= newdeathsH;
	newpop.D += newdeathsH;
	
    const int recoverH = flow(rng, pop.H, newpop.H, (1.0 - p_d) * (1.0 / T_hos));
	newpop.H -= recoverH;
	newpop.R += recoverH;
	
	//hospitalize - frail
//    unsigned int newdeathsH_f=0;
//	std::cout << "newdeathsH_f: " << rrdh * p_d * (1 / T_hos) << "\n";
//	flow(rng, H_f, D, rrdh * p_d * (1 / T_hos), newdeathsH_f);
	
//    unsigned int recoverH_f=0;
//	std::cout << "recoverH_f: " << (1 - rrdh * p_d) * (1 / T_hos) << "\n";
//	flow(rng, H_f, R, (1 - rrdh * p_d) * (1 / T_hos), recoverH_f);
	
    // symptomatic - non-frail
    const int hospitalize = flow(rng, pop.I_s4, newpop.I_s4, p_h  * (1.0 - capacity) * (4.0 / T_sym));
	newpop.I_s4 -= hospitalize;
	newpop.H += hospitalize;
	
    const int newdeathsI_s = flow(rng, pop.I_s4, newpop.I_s4, p_h  * p_d * rrd * capacity * ( 4.0 / T_sym));
	newpop.I_s4 -= newdeathsI_s;
	newpop.D += newdeathsI_s;
	
    const int recoverI_s = flow(rng, pop.I_s4, newpop.I_s4, ( (1.0 - p_h) + p_h  * (1 - p_d * rrd) * capacity) * ( 4.0 / T_sym));
	newpop.I_s4 -= recoverI_s;
	newpop.R += recoverI_s;
	

    const int Is_from3_to_4 = flow(rng, pop.I_s3, newpop.I_s3, (4.0 / T_sym));
	newpop.I_s3 -= Is_from3_to_4;
	newpop.I_s4 += Is_from3_to_4;
	
    const int Is_from2_to_3 = flow(rng, pop.I_s2, newpop.I_s2, (4.0 / T_sym));
	newpop.I_s2 -= Is_from2_to_3;
	newpop.I_s3 += Is_from2_to_3;
		
    const int Is_from1_to_2 = flow(rng, pop.I_s1, newpop.I_s1, (4.0 / T_sym));
	newpop.I_s1 -= Is_from1_to_2;
	newpop.I_s2 += Is_from1_to_2;

	// symptomatic -frail
//   unsigned int hospitalize_f=0;
//	std::cout << "hospitalize_f: " << p_hf * (2 / T_sym) << "\n";
//	flow(rng, I_fs, H_f, p_hf * (1 / T_sym), hospitalize_f);	
	
//   unsigned int newdeaths_f=0;
//	std::cout << "newdeaths_f: " << ( 1 - p_hf ) * ( 2 / T_sym) << "\n";
//	flow(rng, I_fs, D, ( 1 - p_hf ) * ( 1 / T_sym), newdeaths_f);
	
	// asymptomatic
    const int recoverI = flow(rng, pop.I4, newpop.I4, ( 4.0 / T_rec ));
	newpop.I4 -= recoverI;
	newpop.R += recoverI;
	
    const int I_from3_to_4 = flow(rng, pop.I3, newpop.I3, ( 4.0 / T_rec ));
	newpop.I3 -= I_from3_to_4;
	newpop.I4 += I_from3_to_4;
	
    const int I_from2_to_3=flow(rng, pop.I2, newpop.I2, ( 4.0 / T_rec ));
	newpop.I2 -= I_from2_to_3;
	newpop.I3 += I_from2_to_3;
 
    const int I_from1_to_2 = flow(rng, pop.I1, newpop.I1, ( 4.0 / T_rec ));
	newpop.I1 -= I_from1_to_2;
	newpop.I2 += I_from1_to_2;
	
	// infectious - pre-clinical
    const int newasymptomatic = flow(rng, pop.I_p, newpop.I_p, (1.0 - p_s) * ( 1.0 / T_inf ));
	newpop.I_p -= newasymptomatic;
	newpop.I1 += newasymptomatic;
		
//flow(rng, I_p, I1, p_d * ( 1/T_sym ), newasymptomatic);

    const int newsymptomatic = flow(rng, pop.I_p, newpop.I_p, p_s * ( 1.0 / T_inf ));
	newpop.I_p -= newsymptomatic;
	newpop.I_s1 += newsymptomatic;
	
//    unsigned int newinffrail=0;
//	std::cout << "newinffrail: " << pf_val * ( 1/T_inf ) << "\n";
//	flow(rng, I_p, I_fs, pf_val * ( 1/T_inf ), newinffrail);		
	
    // latent
    const int infectious = flow(rng, pop.E, newpop.E, ( 1.0 / T_lat ));
	newpop.E -= infectious;
	newpop.I_p += infectious;
	
    const int infectious_t = flow(rng, pop.E_t, newpop.E_t, ( 1.0 / T_lat ));
	newpop.E_t -= infectious_t;
	newpop.I_t += infectious_t;
	
    // susceptible
    const int newinfection = flow(rng, pop.S, newpop.S, lambda);
	newpop.S -= newinfection;
	newpop.E += newinfection;
	
    // update population
//	std::vector<int> newpop= {(int)S,(int)E,(int)E_t,(int)I_p,(int)I_t, (int)I1,(int)I2,(int)I3,(int)I4,(int)I_s1,(int)I_s2,(int)I_s3,(int)I_s4,(int)I_fs,(int)H,(int)H_f,(int)R,(int)D};
//	std::vector<int> newpop= {(int)S,(int)E,(int)E_t,(int)I_p,(int)I_t, (int)I1,(int)I2,(int)I3,(int)I4,(int)I_s1,(int)I_s2,(int)I_s3,(int)I_s4,(int)H,(int)R,(int)D};
	
//	std::cout << S<<","<<E<<","<<E_t<<","<<I_p<<","<<I_t<<","<<I1<<","<<I2<<","<<I3<<","<<I4<<","<<I_s1<<","<<I_s2<<","<<I_s3<<","<<I_s4<<","<<I_fs<<","<<H<<","<<H_f<<","<<R<<","<<D<<std::endl;
	pop = newpop;
	infection_state.deaths += (newdeathsH + newdeathsI_s) ;
	infection_state.hospital_deaths += newdeathsH ;
	infection_state.detected += (hospitalize + infectious_t);	

	return infection_state;
}

std::vector<double> generate_lambda_vector(int& inf_hosp, const std::vector<double>& parameter_set, const double& u_val, 
			const AgeGroupData& age_data, const std::vector<Compartments>& pops, const bool& shut) 
{

	double p_i = parameter_set[0];	
	double p_hcw = parameter_set[1];	
	double c_hcw = parameter_set[2];	
	double d_val = parameter_set[3];					
	double q_val = parameter_set[4];
	
	//lambda=rep(NA,nrow(pops))
	int n_agegroup = age_data.waifw_norm.size();

	std::vector<double> lambda(n_agegroup);

	double quarantined = 0.0;	// mean proportion reduction in contacts due to quarantine

	std::vector<double> I_mat(n_agegroup); 

	//contact matrix
	std::vector<std::vector<double>> waifw = (shut) ? std::vector<std::vector<double>>(n_agegroup, std::vector<double> (n_agegroup)) : age_data.waifw_norm; 

	//age-dependent transmission rate
	std::vector<std::vector<double>> beta(n_agegroup, std::vector<double> (n_agegroup));

  if(shut)
  {
    for ( int from{0}; from < n_agegroup; ++from) 
    {
      for ( int to{0}; to < n_agegroup; ++to) 
      {
        //contact network during shutdown period, assuming a proportion d will properly do it
        waifw[from][to] = (1.0-d_val) * age_data.waifw_sdist[from][to] + (d_val) * age_data.waifw_home[from][to];
      }
    }
  }
  
  for ( int from{0}; from < n_agegroup; ++from) {
    for ( int to{0}; to < n_agegroup; ++to) 
    {
      if(waifw[from][to]>0)
      {
        quarantined += age_data.waifw_home[from][to]/waifw[from][to];
      } 
			
      else
      {
         quarantined += 0;
      }
      beta[from][to] = waifw[from][to] * p_i ;
	}
  } 
  
  //mean proportion reduction in contacts due to quarantine accounting for compliance
  quarantined = 1 - (quarantined / (double)(n_agegroup * n_agegroup) ) * (1 - q_val);
  
  //compute the pressure from infectious groups, normalizes by group size
  //compute the pressure from hospitalised
  for ( int from{0}; from < n_agegroup; ++from) {
	if(from < (n_agegroup-1)) { //n-agegroup-1 to not account for hcw (as different contact)
	 	// are considered those that are infectious pre-clinical, asymp, tested or symptomatic only
		// asym are infectious pre-clinical and asymptomatic
		// sym are infectious plus those that are tested positive
		int tot_asym = pops[from].I_p + u_val * (pops[from].I1 + pops[from].I2 + pops[from].I3 + pops[from].I4 );
		int tot_sym = pops[from].I_t + pops[from].I_s1 + pops[from].I_s2 + pops[from].I_s3 + pops[from].I_s4;
		
		I_mat[from] = (double)tot_asym + (1-quarantined) * (double)tot_sym;
		//normalisation shouldnt account for dead individuals	
		const int tot_pop = accumulate_compartments(pops[from]);
		I_mat[from] = (double)I_mat[from] / (double)(tot_pop); 
				
	}
	//sum up of number of hospitalised cases (frail and non-frail)
	inf_hosp+= pops[from].H ;
  }

  //sum up infectious pressure from each age group
  for ( int to{0}; to < (n_agegroup); ++to) {
	  for(int from{0}; from < (n_agegroup-1); ++from){ //assume perfect self isolation and regular testing of infected HCW if infected
	  	lambda[to] += ( beta[to][from]*I_mat[from] );
	  }
  }

  const int tot_pop = accumulate_compartments(pops[n_agegroup-1]);
  //int rem_pop = pops[n_agegroup-1][4]+pops[n_agegroup-1][5]+pops[n_agegroup-1][6]+pops[n_agegroup-1][8];
  lambda[n_agegroup-1] = p_hcw * c_hcw * ((double)inf_hosp/(double)(tot_pop) );

  return lambda;
}

/*void flow(Random::RNGInterface::Sptr rng, unsigned int& pop_from, unsigned int& pop_to, double rate, unsigned int& outs){
    outs = rng->Poisson(rate * (double)pop_from); //symptomatic become hospitalized
    outs = std::min(pop_from, outs);
	pop_to += outs;
	pop_from -= outs;
}*/

int flow(Random::RNGInterface::Sptr rng, const int& pops_from_val, const int& pops_new_from_val, const double& rate)
{
	int outs = rng->Poisson(rate * static_cast<double>(pops_from_val)); //symptomatic become hospitalized
	outs = std::min(pops_new_from_val, outs);
	return outs;
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
	for (int i{0}; i < nPar; ++i) {
		
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

//reduce vector of daily  value into vector of weekly value
static void weekly(std::vector<int>& reduced, const std::vector<int>& original){

	int weekly_value=0;
	for(unsigned int ttime =0; ttime< original.size(); ++ttime){
		weekly_value+=original[ttime];
				
		if( (ttime % 7) == 6){
			reduced.push_back(weekly_value);
			weekly_value=0;
		}
	}
}

} // namespace Model
} // namespace EERAModel
