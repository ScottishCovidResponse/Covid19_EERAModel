#include "Model.h"
#include "IO.h"
#include "ModelTypes.h"
#include "DistanceComputation.h"
#include "Utilities.h"
#include "FittingProcess.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <numeric>
#include <random>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

namespace EERAModel {
namespace Model {

static void infspread(gsl_rng * r, std::vector<int>& pop, int& deaths, int& deathsH, int& detected, 
				::EERAModel::params fixed_parameters, std::vector<double> parameter_set, std::vector<double> cfr_tab,
				double pf_val, double lambda);

static void Lambda(std::vector<double> &lambda, std::vector<double> parameter_set,std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist,std::vector<std::vector<double>> waifw_home, 
				std::vector<std::vector<int>> pops, int shut);

static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage,
				std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist, std::vector<std::vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, std::vector<int> agenums, 
				double tau, gsl_rng * r, std::vector<int> &sim_status,std::vector<std::vector<int>> &ends,
				std::vector<int> &death_status,std::vector<int> &deathH_status);

static void compute_incidence(std::vector<int> v, std::vector<int>& r_val);

static void correct_incidence(std::vector<int>& v, std::vector<int> cumv);

void Run(EERAModel::ModelInputParameters& modelInputParameters,
         EERAModel::Observations observations,
		 gsl_rng* r,
		 std::mt19937& gen,
		 const std::string& outDirPath) {
	/*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;

	// declare vectors of inputs
	std::vector<int > obsHosp;
	std::vector<int > obsDeaths;

	// variable declaration
	double valMax, valmin;

    std::cout << "[Settings]:\n";
	std::cout<< "number of parameters tested: "<< modelInputParameters.nPar << std::endl;
    std::cout<< "seeding method: "<< modelInputParameters.seedlist.seedmethod<<  std::endl;
	if(modelInputParameters.seedlist.seedmethod == "random"){
		std::cout << "number of seed: " << modelInputParameters.seedlist.nseed << std::endl;
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		std::cout<< "duration of the high risk period: " << modelInputParameters.seedlist.hrp << std::endl;
	}
	
	//keep information for the health board if interest
	std::vector<double> pf_byage = observations.pf_pop[modelInputParameters.herd_id - 1];//define frailty structure of the shb of interest.

	//create vector of fixed parameters		
	std::vector<params> fixed_parameters(observations.waifw_norm.size());	
	for (unsigned int var = 0; var < fixed_parameters.size(); ++var) {
		fixed_parameters[var] = modelInputParameters.paramlist;
	}
	fixed_parameters[0].p_s = fixed_parameters[0].juvp_s;

	//Separate case information for each herd_id
	std::vector<int> obsHosp_tmp, obsDeaths_tmp;
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
			std::cout << "Warning!! Unknown seeding method - applying _random_ seed method\n";
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
	std::vector<double> agedist = observations.age_pop[modelInputParameters.herd_id - 1];//define age structure of the shb of interest. the -1 is to account for difference in number of rows between two datasets (age_pop does not have a row with column name)
	std::vector<int> agenums;
	for (unsigned int var = 0; var < agedist.size(); ++var) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(agedist[var] * (Npop - N_hcw))); 

	}		
	//push back the number of hcw in the shb of interest
	agenums.push_back(N_hcw);
	
    std::cout << "[Health Board settings]:\n";
	std::cout << "    SHB id: " << modelInputParameters.herd_id <<'\n';
	std::cout << "    Population size: " << Npop << '\n';
	std::cout << "    Number of HCW: " << N_hcw << '\n';
	std::cout << "    Simulation period: " << duration << "days\n";
	std::cout << "    time step: " << modelInputParameters.tau << "days\n";
	//declare vectors for priors
	std::vector<double> flag1, flag2;
	
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
	std::vector<particle > particleList,particleList1;
	std::vector<double> weight_val;

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
	 std::cout << "[Simulations]:\n";
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
		std::discrete_distribution<int> weight_distr(weight_val.begin(), weight_val.end());
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
					FittingProcess::parameter_select_first_step(outs_vec.parameter_set,flag1, flag2, r, modelInputParameters.nPar);
				} else {
					//sample 1 particle from the previously accepted particles and given their weight (also named "importance sampling")
				    pick_val = weight_distr(gen);
				    FittingProcess::parameter_select_nsteps(outs_vec.parameter_set, modelInputParameters.nPar, r, particleList, pick_val, vlimitKernel,vect_min,vect_Max);
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
							FittingProcess::weight_calc(smc,Nparticle, particleList, outs_vec, vlimitKernel, modelInputParameters.nPar);
							particleList1.push_back(outs_vec);
							++counter;		//count the number of accepted particles
							if(counter % 10 == 0) std::cout << "|" << std::flush;
							//std::cout << counter << " " ;
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
		std::cout << "\nStep:" << smc
			<< ", <number of accepted particles> " << Nparticle
			<< "; <number of simulations> " << nsim_val
			<< "; <computation time> " <<  time_taken
			<< " seconds.\n";

		//break the ABC-smc at the step where no particles were accepted
		if (Nparticle > 0) {
			IO::WriteOutputsToFiles(smc, modelInputParameters.herd_id, Nparticle,
				modelInputParameters.nPar, particleList1, outDirPath);
		}
	}

	//output on screen the overall computation time
	std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
}

void model_select(int smc, EERAModel::particle &outvec, std::vector<params> fixed_parameters,
	std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage, 
	std::vector<std::vector<double>> waifw_norm, std::vector<std::vector<double>> waifw_sdist,
	std::vector<std::vector<double>> waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, int Npop, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths) {

	//---------------------------------------
	// the root model
	//---------------------------------------
	//declare the vector of outputs		int t, t_int, i, j, yr,yr_int,increment;
	std::vector<int> sim_status;
	std::vector<int> death_status;
	std::vector<int> deathH_status;
	std::vector<std::vector<int>> ends;
	
	//main loop
	my_model(outvec.parameter_set, fixed_parameters, cfr_byage, pf_byage,waifw_norm, waifw_sdist, waifw_home,
					duration, seedlist, day_shut, Npop,agenums, tau, r, sim_status, ends,
					death_status,deathH_status);

	//---------------------------------------
	// compute the  sum of squared errors
	//---------------------------------------
	double sum_sq_cases = 1000000.0,sum_sq_deaths = 1000000.0,nsse_cases=1000000.0,nsse_deaths=1000000.0;
	std::vector<int> obs_case = obsHosp;
	std::vector<int> obs_death = obsDeaths;
	sum_sq_cases = ::EERAModel::DistanceComputation::sse_calc_int(sim_status,obs_case);
		
	sum_sq_deaths= ::EERAModel::DistanceComputation::sse_calc_int(deathH_status,obs_death);
	//---------------------------------------
	// compute the deviation of the fit expressed as % total number of cases
	//---------------------------------------
	//total numbers of cases and deaths (at hospital)
	int sum_obs_cases = std::accumulate(obsHosp.begin(), obsHosp.end(), 0);
	int sum_obs_deaths = std::accumulate(obsDeaths.begin(), obsDeaths.end(), 0);
	
	nsse_cases = sqrt(sum_sq_cases)/ (double)sum_obs_cases;
	nsse_deaths = sqrt(sum_sq_deaths)/ (double)sum_obs_deaths;
 	
	//cout << "nsse_cases: " << nsse_cases << " , nsse_deaths: " << nsse_deaths <<'\n';

	//---------------------------------------
	// Return a vector with all selection measures
	//---------------------------------------
	//	outvec.sum_sq = sum_sq;
	outvec.nsse_cases = nsse_cases;
	outvec.nsse_deaths = nsse_deaths;
	for (unsigned int ii = 0; ii < sim_status.size(); ++ii) {
		outvec.simu_outs.push_back(sim_status[ii]);
		outvec.death_outs.push_back(deathH_status[ii]);
	}
	for (unsigned int ii = 0; ii < ends.size(); ++ii) {
		outvec.end_comps.push_back(ends[ii]);
	}		

}

 
void select_obs(int& Npop, int& t_index, int& duration, int& day_intro, int& day_shut, 
	std::vector<int>& obsHosp_tmp, std::vector<int>& obsDeaths_tmp, 
	std::vector<std::vector<int> > data_tmp, std::vector<std::vector<int> > death_tmp, int herd_id, 
	int time_back) {
			
	int maxTime=0;
	std::vector<int> seqTime;
		
	//define population size
	Npop = data_tmp[herd_id][0];
	
	//create the vector of cases (cummulative) and define time of first detection (index)
	for (unsigned int iv = 1; iv < data_tmp[0].size(); ++iv) {
		if(data_tmp[herd_id][iv]>=0){
			seqTime.push_back(data_tmp[0][iv]);
			maxTime = ::EERAModel::Utilities::max_calc_int(maxTime,(int)data_tmp[0][iv]);
			obsHosp_tmp.push_back(data_tmp[herd_id][iv]);
			obsDeaths_tmp.push_back(death_tmp[herd_id][iv]);
			if(data_tmp[herd_id][iv]>0 && t_index<0){
				t_index = (int)data_tmp[0][iv];//identify when the first case is detected/hospitalised 
			}		
		}
	}

	//identify the first day of infectiousness for the index case, which will be the start of our simulation, and add days in the observation, and define duration of the disease process
	int intro = (t_index+1)-(time_back);

	//define the duration of the study period
	if(intro < 0) duration = maxTime - intro + 1;
	else duration = maxTime + 1;

	//define the day of the incursion
	if(intro < 0 ) day_intro = 0;
	else day_intro=intro;
		
	//add the extra information on the observations
	if(duration > obsHosp_tmp.size()){
		int extra_time = duration - obsHosp_tmp.size();
		std::vector<int> extra_cases, extra_deaths;
		for (int extra = 0; extra < extra_time; ++extra) {
			extra_cases.push_back(0);
			extra_deaths.push_back(0);
		}		
		for (unsigned int iv = 0; iv < obsHosp_tmp.size(); ++iv) {
			extra_cases.push_back(obsHosp_tmp[iv]);
			extra_deaths.push_back(obsDeaths_tmp[iv]);
		}		
		day_shut = day_shut+ extra_time;
		obsHosp_tmp = extra_cases;
		obsDeaths_tmp = extra_deaths;
		extra_cases.clear();	
		extra_deaths.clear();	
	}
	
	std::cout << "Number of days of obs cases: " << obsHosp_tmp.size() << std::endl;
	std::cout << "Number of days of obs deaths: " << obsDeaths_tmp.size() << std::endl;
	//transform  cumulative numbers into incident cases
	std::vector<int> obsHosp_tmp2;
	compute_incidence(obsHosp_tmp,obsHosp_tmp2);
	correct_incidence(obsHosp_tmp2,obsHosp_tmp);	
	obsHosp_tmp = obsHosp_tmp2;
	obsHosp_tmp2.clear();	
	
	//transform  cumulative numbers into incident deaths
	std::vector<int> obsDeaths_tmp2;
	compute_incidence(obsDeaths_tmp,obsDeaths_tmp2);
	correct_incidence(obsDeaths_tmp2,obsDeaths_tmp);	
	obsDeaths_tmp = obsDeaths_tmp2;
	obsDeaths_tmp2.clear();	
}

static void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage,
				std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist, std::vector<std::vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, std::vector<int> agenums, 
				double tau, gsl_rng * r, std::vector<int> &sim_status,std::vector<std::vector<int>> &ends,
				std::vector<int> &death_status,std::vector<int> &deathH_status) {

	int n_agegroup = waifw_norm.size();
	int inLockdown = 0;
	double seed_pop[6];
	

	//sets up an array for the population at each timestep in each age and disease category	
	//also set up the age distribution of old ages as target for disease introduction
	std::vector<std::vector<int>> poparray(n_agegroup, std::vector<int>(9));	  
	for ( int age = 0; age < (n_agegroup); ++age) {
		for ( int st = 0; st < 9; ++st) {
			if(st ==0) {
				poparray[age][st]=agenums[age]; //set up the starting population as fully susceptible
			} else {
				poparray[age][st] = 0;
			}
		}
		if(age >0 & age < (n_agegroup-1) ){ //seed only occurs in >20yo and not in HCW
			seed_pop[age-1] = (double)agenums[age];
		}
	}


	//introduce disease at t=0. if seedmethod != "background"
	//WARNING!! MINOR BUG HERE: gsl multinomial sometimes freeze for unknown reasons. replaced by random flat, picking a value where to seed infection. will only works because seed=1.
	if(seedlist.seedmethod != "background"){
	//	size_t k = sizeof(seed_pop);
	//	unsigned int startdist[k];
		std::vector<int> startdist = {0,0,0,0,0,0};
		startdist[(int) gsl_ran_flat(r, 0, 6)] = seedlist.nseed;
	//	gsl_ran_multinomial(r, k, startdz, seed_pop,startdist); //distribute the diseased across the older age categories

		for ( int age =1; age < (n_agegroup-1); ++age) {
			poparray[age][0] -=  startdist[age-1];//take diseased out of S
			poparray[age][3] +=  startdist[age-1];//put diseased in I	
		}	
	}


  	//initialize saving of the detection for t=0
	sim_status.push_back(0); 
	death_status.push_back(0); 
	deathH_status.push_back(0); 
	//run the simulation
	for (int tt = 1; tt < (int)ceil(duration/tau); ++tt) {
		//initialize return value
		int deaths=0;
		int deathsH=0;
		int detected=0;
		std::vector<int> popcomp;
		
		//identify the lock down
		if(tt >= day_shut)  inLockdown = 1;
		
		
		//introduce disease from background infection until lockdown
	    if(inLockdown<1){
			if(seedlist.seedmethod == "background"){
				double bkg_lambda = parameter_set[parameter_set.size()-1];
				unsigned int startdz = gsl_ran_poisson(r, (double)Npop * bkg_lambda);	//how many diseased is introduced in each given day before lockdown
		
			//	size_t k = sizeof(seed_pop);
			//	unsigned int startdist[k];
				std::vector<int> startdist = {0,0,0,0,0,0};
				int pickcomp = (int)gsl_ran_flat(r, 0, 6);
				while(poparray[pickcomp][0]<1){
					pickcomp = (int)gsl_ran_flat(r, 0, 6);
				} 
				startdist[pickcomp] = startdz;
			//	gsl_ran_multinomial(r, k, startdz, seed_pop,startdist); //distribute the diseased across the older age categories

				for ( int age =1; age < (n_agegroup-1); ++age) {
					int nseed = std::min(poparray[age][0], startdist[age-1]);
					poparray[age][0] -=  nseed;//take diseased out of S
					poparray[age][1] +=  nseed;//put diseased in E	
				}	
			}
	    }


		//compute the forces of infection
		std::vector<double> lambda(n_agegroup);
		Lambda(lambda, parameter_set,waifw_norm, waifw_sdist,waifw_home, poparray, inLockdown);	
	
		//step each agegroup through infections
		for ( int age = 0; age < (n_agegroup); ++age) {	
			infspread(r, poparray[age], deaths, deathsH, detected,fixed_parameters[age],parameter_set,
				cfr_byage[age],pf_byage[age],lambda[age]);
		}

		//cout << tt << " , " <<  detected << " , " << obsHosp[tt]<< " , "<< deaths << '\n';

		sim_status.push_back(detected); 
		death_status.push_back(deaths); 
		deathH_status.push_back(deathsH);
    }	
	//save the population in each epi compt for the last day
	for ( int age = 0; age < n_agegroup; ++age) {
		ends.push_back(poparray[age]);
	}
}

static void infspread(gsl_rng * r, std::vector<int>& pop, int& deaths, int& deathsH, int& detected, 
	params fixed_parameters, std::vector<double> parameter_set, std::vector<double> cfr_tab, double pf_val, double lambda){
		
    unsigned int S=pop[0], E=pop[1], E_t=pop[2], I=pop[3], I_t=pop[4], I_s=pop[5], H=pop[6], R=pop[7], D=pop[8] ;
	
	double T_lat= fixed_parameters.T_lat;
	double p_s= fixed_parameters.p_s;
	double T_inf= fixed_parameters.T_inf;
	double T_rec= fixed_parameters.T_rec;
	double T_sym= fixed_parameters.T_sym;
	double T_hos= fixed_parameters.T_hos;
	
	double p_h= cfr_tab[0];
	double cfr= cfr_tab[1];	
	double p_d= cfr_tab[2];	
	//double p_dc= cfr_tab[3];
		
	double rrdh = parameter_set[5];
	double rrdc = parameter_set[6];
	double rrh = parameter_set[7];	
	
	
	double C_val = pf_val  * rrh * p_h ;
	double D_val = (1 - pf_val) * p_h;
	double A_val = ( ( C_val * rrdh + D_val ) / ( C_val + D_val ) ) * p_d;
	double B_val = ( (C_val * ( 1 -  rrdh * p_d ) + D_val * ( 1 - p_d ) ) / (C_val+D_val) );
	double E_val = (( pf_val - C_val ) * rrdc + ((1 - pf_val ) - D_val ) ) * cfr;
	double F_val = ( pf_val - C_val ) * (1 - rrdc * cfr ) + ((1 - pf_val ) - D_val ) * ( 1 - cfr ) ;
	
    // hospitalized
    unsigned int newdeathsH = gsl_ran_poisson(r, A_val * (1 / T_hos) * (double)H);
	newdeathsH =  std::min(H,newdeathsH);
    D += newdeathsH;
	H -= newdeathsH;
    unsigned int recoverH = gsl_ran_poisson(r,  B_val * (1 / T_hos) * (double)H);//hospitalized recover
	recoverH= std::min(H,recoverH) ;
    R += recoverH;
	H -= recoverH;
	
    // symptomatic
    unsigned int hospitalize = gsl_ran_poisson(r, (C_val + D_val) * (1 / T_sym) * (double)I_s); //symptomatic become hospitalized
    hospitalize = std::min(I_s, hospitalize);
	H += hospitalize;
	I_s -= hospitalize;
	
    unsigned int newdeathsC = gsl_ran_poisson(r, E_val * ( 2 / T_hos) * (double)I_s); //symptomatic die at home / in communities	
    newdeathsC = std::min(I_s, newdeathsC);
    D += newdeathsC;
	I_s -= newdeathsC;
	
    unsigned int recoverI_s = gsl_ran_poisson(r, F_val * ( 1 / T_rec) * (double)I_s); //symptomatic recover at home / in communities	
    recoverI_s = std::min(I_s, recoverI_s);
    R += recoverI_s;
	I_s -= recoverI_s;
	
    // infectious
    unsigned int symptomsI = gsl_ran_poisson(r, p_s * ( 1/T_inf ) * (double)I ) ; //infectious become symptomatic
    symptomsI = std::min(I, symptomsI);
	unsigned int symptomsI_t=gsl_ran_poisson(r, p_s * ( 1/T_inf ) * (double)I_t); //tested infectious become symptomatic
    symptomsI_t=std::min(I_t,symptomsI_t);
	I_s += (symptomsI+symptomsI_t);
	I -= symptomsI;
	I_t -= symptomsI_t;
	
    unsigned int recoverI = gsl_ran_poisson(r, (1 - p_s) * ( 1/T_inf ) * (double)I); //infectious recover
    recoverI = std::min(I, recoverI);
	unsigned int recoverI_t = gsl_ran_poisson(r, (1 - p_s) * ( 1/T_inf ) * (double)I_t); //tested infectious recover
	recoverI_t = std::min(I_t, recoverI_t);
    R += (recoverI+recoverI_t) ;
	I -= recoverI;
	I_t -= recoverI_t;
	
    // latent
    unsigned int infectious = gsl_ran_poisson(r, (double)E/T_lat); //latent become infectious
	infectious = std::min(E, infectious);
    I += infectious;
	E -= infectious;
    unsigned int infectious_t = gsl_ran_poisson(r, (double)E_t/T_lat); //tested latent become infectious
    infectious_t = std::min(E_t, infectious_t);
	I_t += infectious_t;
	E_t -= infectious_t;
	
    // susceptible
    unsigned int newinfection= gsl_ran_poisson(r, (double)S*lambda); //susceptible become infected
    newinfection= std::min(S, newinfection);
	S -= newinfection;
	E += newinfection;
	
    // recover population
	std::vector<int> newpop= {(int)S,(int)E,(int)E_t,(int)I,(int)I_t,(int)I_s,(int)H,(int)R,(int)D};
	pop = newpop;
	deaths += (newdeathsC + newdeathsH);
	deathsH += newdeathsH;
	detected += (H + I_t);	
	
}

static void Lambda(std::vector<double> &lambda, std::vector<double> parameter_set,std::vector<std::vector<double>> waifw_norm,
			std::vector<std::vector<double>> waifw_sdist,std::vector<std::vector<double>> waifw_home, 
			std::vector<std::vector<int>> pops, int shut) {

  double p_i = parameter_set[0];	
  double p_hcw = parameter_set[1];	
  double c_hcw = parameter_set[2];	
  double d_val = parameter_set[3];					
  double q_val = parameter_set[4];					
	
  //lambda=rep(NA,nrow(pops))
  int n_agegroup = waifw_norm.size();
  double quarantined = 0.0;//mean proportion reduction in contacts due to quarantine
  int inf_hosp=0;
  std::vector<double> I_mat(n_agegroup); 
  //contact matrix
  std::vector<std::vector<double>> waifw(n_agegroup, std::vector<double> (n_agegroup)); 
  //age-dependent transmission rate
  std::vector<std::vector<double>> beta(n_agegroup, std::vector<double> (n_agegroup)); 
  
  for ( int from = 0; from < n_agegroup; ++from) {
	  for ( int to = 0; to < n_agegroup; ++to) {
		  if(shut==0){
		  	//contact network during normal period
		    waifw[from][to]  = waifw_norm[from][to] ;
		  } else {
			//contact network during shutdown period, assuming a proportion d will properly do it
		  	waifw[from][to] = (1.0-d_val) * waifw_sdist[from][to] + (d_val) * waifw_home[from][to];
		  }
		  beta[from][to] = waifw[from][to] * p_i ;
		  if(waifw[from][to]>0) {
			  quarantined += waifw_home[from][to]/waifw[from][to];
		  } else {
			  quarantined += 0;
		  }
	  }
  } 
  
  //mean proportion reduction in contacts due to quarantine accounting for compliance
  quarantined = 1 - (quarantined / (double)(n_agegroup * n_agegroup) ) * (1 - q_val);
  
  //compute the pressure from infectious groups, normalizes by group size
  //compute the pressure from hospitalised, normalized by group size
  for ( int from = 0; from < n_agegroup; ++from) {
	if(from < (n_agegroup-1)) { //n-agegroup-1 to not account for hcw (as different contact)
	 	// are considered those that are infectious asymp, tested or symptomatic only
		I_mat[from] = (double)pops[from][3] + (1-quarantined) * ( (double)pops[from][4] + (double)pops[from][5] );
		//normalisation shouldnt account for dead individuals	
		int tot_pop = std::accumulate(pops[from].begin(), pops[from].end(), 0);
		I_mat[from] = (double)I_mat[from] / (double)(tot_pop- pops[from][pops[from].size()-1] ); 
				
	}
	//sum up of number of hospitalised cases 
	inf_hosp+= pops[from][6];
  }

  //sum up infectious pressure from each age group
  for ( int from = 0; from < (n_agegroup); ++from) {
	  for(int ii = 0; ii < (n_agegroup-1); ++ii){
	  	lambda[from] += ( beta[from][ii]*I_mat[ii] );
	  }
  }

  int tot_pop = std::accumulate(pops[n_agegroup-1].begin(), pops[n_agegroup-1].end(), 0);
  int rem_pop = pops[n_agegroup-1][4]+pops[n_agegroup-1][5]+pops[n_agegroup-1][6]+pops[n_agegroup-1][8];
  lambda[n_agegroup-1] = p_hcw * c_hcw * ((double)inf_hosp/(double)(tot_pop - rem_pop) );
 
}


//transform the timeseries of cummulative cases into incidence
static void compute_incidence(std::vector<int> v, std::vector<int>& r_val){
	for (unsigned int xx = 0; xx < v.size(); ++xx) {
		int tmp_inc = 0;
		if(xx>0) tmp_inc = v[xx] - v[xx-1];
		r_val.push_back(tmp_inc);
	}	
}

//correct the timeseries to avoid negative incidence records
static void correct_incidence(std::vector<int>& v, std::vector<int> cumv){
	if ( any_of(v.begin(), v.end(), [](int i){return i<0;}) ){
		for (unsigned int ii = 1; ii < v.size(); ++ii) {
			if(v[ii]<0){
				int case_bef = cumv[ii-1];//look next day number of cases
				int case_aft = cumv[ii+1];//look at the previous number of cases	
				
      		  	//check if the next tot cases is bigger than tot cases before
    			//if(bigger) compute difference and remove it from day before
		        if(case_bef<case_aft){
		          if(v[ii-1] >= abs(v[ii])){
		            cumv[ii-1] += v[ii];
		          } else {
		          	cumv[ii]=0;
		          }
		        } else{
		          //if smaller
		          //check when the day prior today was smaller than the next day and remove extra cases to this day
		          int counter_check=1;
		          while( cumv[ii-counter_check]>case_aft){
		            ++counter_check ;
		          }
		          if(v[ii-counter_check+1] >= abs(v[ii])){
				  	for ( int jj = (ii-counter_check+1); jj < (ii); ++jj) {
						cumv[jj] += v[ii];
					}			
		          } else {
		            cumv[ii]=0;
		          }
		        }
			}

		}
		//recompute the incident cases
		v.clear();
		compute_incidence(cumv,v);
	}
}

} // namespace Model
} // namespace EERAModel