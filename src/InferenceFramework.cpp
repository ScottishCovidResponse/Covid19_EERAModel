#include "InferenceFramework.h"
#include "Observations.h"
#include "IO.h"
#include "ModelCommon.h"
#include <functional>
#include <algorithm>

namespace EERAModel {
namespace Inference {

InferenceFramework::InferenceFramework(Model::ModelInterface::Sptr model,
    const ModelInputParameters& modelInputParameters,
    const InputObservations& observations,
    Random::RNGInterface::Sptr rng,
    const std::string& outDir,
    Utilities::logging_stream::Sptr log)
    : model_(model),
      modelInputParameters_(modelInputParameters),
      observations_(observations),
      rng_(rng),
      outDir_(outDir),
      log_(log),
      toleranceLimits_(modelInputParameters.toleranceLimit) {
    
    std::vector<double> flag1 = {
		modelInputParameters.prior_pinf_shape1,
	 	modelInputParameters.prior_phcw_shape1,
		modelInputParameters.prior_chcw_mean, 
		modelInputParameters.prior_d_shape1, 
		modelInputParameters.prior_q_shape1,
		modelInputParameters.prior_ps_shape1,
		modelInputParameters.prior_rrd_shape1,		
		modelInputParameters.prior_lambda_shape1
	};	
	std::vector<double> flag2 = {
		modelInputParameters.prior_pinf_shape2, 
		modelInputParameters.prior_phcw_shape2, 
		modelInputParameters.prior_chcw_mean, 
		modelInputParameters.prior_d_shape2, 
		modelInputParameters.prior_q_shape2,
		modelInputParameters.prior_ps_shape2,
		modelInputParameters.prior_rrd_shape2,		
		modelInputParameters.prior_lambda_shape2
	};

    inferenceParameterGenerator_ = std::make_shared<InferenceParameterGenerator>(rng, flag1, flag2);
}

int InferenceFramework::GetTimeOffSet(const ModelInputParameters& modelInputParameters)
{
	int time_back = 0;
	if(modelInputParameters.seedlist.seedmethod == "background") {
		time_back = modelInputParameters.seedlist.hrp;
	} else {
		time_back = modelInputParameters.paramlist.T_inf + modelInputParameters.paramlist.T_sym;
		
	}

	return time_back;
}

void InferenceFramework::Run()
{
    const unsigned int nInferenceParams = Model::ModelParameters::NPARAMS;

    /*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;

    (*log_) << "[Settings]:\n";
	(*log_) << "    number of parameters tested: "<< nInferenceParams << std::endl;
    (*log_) << "    seeding method: "<< modelInputParameters_.seedlist.seedmethod<<  std::endl;
	if (modelInputParameters_.seedlist.seedmethod == "random"){
		(*log_) << "    number of seed: " << modelInputParameters_.seedlist.nseed << std::endl;
	} else if(modelInputParameters_.seedlist.seedmethod == "background"){
		(*log_) << "    duration of the high risk period (hrp): " << modelInputParameters_.seedlist.hrp << std::endl;
	}
    (*log_) << "    model structure: " << 
        ((modelInputParameters_.model_structure == ModelStructureId::ORIGINAL) ? "Original" : "Irish") << std::endl;

    (*log_) << "[Fixed parameter values]:\n";
	(*log_) << "    latent period (theta_l): " << modelInputParameters_.paramlist.T_lat <<std::endl;
	(*log_) << "    pre-clinical period (theta_i): " << modelInputParameters_.paramlist.T_inf <<std::endl;
	(*log_) << "    asymptomatic period (theta_r): " << modelInputParameters_.paramlist.T_rec <<std::endl;
	(*log_) << "    symptomatic period (theta_s): " << modelInputParameters_.paramlist.T_sym <<std::endl;
	(*log_) << "    hospitalisation stay (theta_h): " << modelInputParameters_.paramlist.T_hos <<std::endl;
	(*log_) << "    pre-adult probability of symptoms devt (p_s[0]): " << modelInputParameters_.paramlist.juvp_s <<std::endl;
	(*log_) << "    bed capacity at hospital (K): " << modelInputParameters_.paramlist.K <<std::endl;
	(*log_) << "    relative infectiousness of asymptomatic (u): " << modelInputParameters_.paramlist.inf_asym <<std::endl;
	
	//keep information for the health board if interest
	const std::vector<double> pf_byage = observations_.pf_pop[modelInputParameters_.herd_id - 1];//define frailty structure of the shb of interest.

	const AgeGroupData per_age_data = {observations_.waifw_norm, observations_.waifw_home, observations_.waifw_sdist, 
										observations_.cfr_byage, pf_byage};

	//create vector of fixed parameters		
	std::vector<params> fixed_parameters = Model::BuildFixedParameters(observations_.waifw_norm.size(),
        modelInputParameters_.paramlist);
    
	const int time_back = GetTimeOffSet(modelInputParameters_);

	int population = Model::GetPopulationOfRegion(observations_, modelInputParameters_.herd_id);
	
	const std::vector<int>& regionalCases = observations_.cases[modelInputParameters_.herd_id];
	const std::vector<int>& regionalDeaths = observations_.deaths[modelInputParameters_.herd_id];
	const std::vector<int>& timeStamps = observations_.cases[0];

	const Observations::ObsSelect obs_selections = Observations::SelectObservations(
		modelInputParameters_.day_shut, timeStamps, regionalCases,
		regionalDeaths, time_back, log_);
	
	modelInputParameters_.seedlist.day_intro = obs_selections.sim_time.day_intro;

	int N_hcw = Model::ComputeNumberOfHCWInRegion(population, modelInputParameters_.totN_hcw, observations_);

	std::vector<int> agenums = Model::ComputeAgeNums(modelInputParameters_.herd_id, population, N_hcw, observations_);
	
    (*log_) << "[Health Board settings]:\n";
	(*log_) << "    SHB id: " << modelInputParameters_.herd_id <<'\n';
	(*log_) << "    Population size: " << population << '\n';
	(*log_) << "    Number of HCW: " << N_hcw << '\n';
	(*log_) << "    Simulation period: " << obs_selections.sim_time.duration << "days\n";
	(*log_) << "    time step: " << modelInputParameters_.tau << "days\n";

	const int n_sim_steps = static_cast<int>(ceil(obs_selections.sim_time.duration/modelInputParameters_.tau));

	// currentParticles is the list of particles accepted on the current ABC-smc loop
    // previousParticles is the list accepted on the previous loop
	std::vector<particle> currentParticles, previousParticles;

	//initialise the number of accepted particles
	int prevAcceptedParticleCount = 0;

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	(*log_) << "[Simulations]:\n";
	for (int smc = 0; smc < modelInputParameters_.nsteps; ++smc) {//todo:

		//the abort statement for keeping the number of particles less than 1000
		bool aborting = false;

		//initialise the counter for the number of simulation prior maximum particle is reached
		int nsim_count = 0;

		//initialise the number of accepted particles
		int acceptedParticleCount = 0;

		//initialise the weight and particle lists
		std::vector<double> vect_Max(nInferenceParams, 0.0);
		std::vector<double> vect_min(nInferenceParams, 0.0);
		std::vector<double> vlimitKernel(nInferenceParams, 0.0);

        
        if (smc > 0) {
			//update the vectors
			previousParticles = currentParticles;
			currentParticles.clear();

			ComputeKernelWindow(nInferenceParams, previousParticles, 
				modelInputParameters_.kernelFactor, vlimitKernel, vect_Max, vect_min);
		}

		std::discrete_distribution<int> weight_distr = ComputeWeightDistribution(previousParticles);

	/*---------------------------------------
	 * simulate the infection data set
	 *---------------------------------------*/
		//omp_set_num_threads(num_threads); // Use maximum num_threads threads for all consecutive parallel regions
		//#pragma omp parallel for default(shared), private(pick_val), firstprivate(weight_distr)
		for (int sim = 0; sim < modelInputParameters_.nSim; ++sim) {//todo
			//abort statement if number of accepted particles reached nParticLimit particles
			//#pragma omp flush (aborting)
			if (!aborting) {
				//Update progress
				if (acceptedParticleCount >= modelInputParameters_.nParticalLimit) {
					aborting = true;
					//#pragma omp flush (aborting)
				}

				//declare and initialise output variables
				particle outs_vec;
				outs_vec.iter = sim;				
//				outs_vec.sum_sq = 1.e06;
				for (int i{0}; i < nInferenceParams; ++i) {
					outs_vec.parameter_set.push_back(0.0);
				}
				//pick the values of each particles
				if (smc==0) {
					//pick randomly and uniformly parameters' value from priors

                    outs_vec.parameter_set = inferenceParameterGenerator_->GenerateInitial();
					
				} else {
					// sample 1 particle from the previously accepted particles
                    // and given their weight (also named "importance sampling")
				    int pick_val = weight_distr(rng_->MT19937());
				    outs_vec.parameter_set = inferenceParameterGenerator_->GenerateWeighted(
                        previousParticles[pick_val].parameter_set, vlimitKernel, vect_Max, vect_min);
				}

				//run the model and compute the different measures for each potential parameters value
				ModelSelect(outs_vec, fixed_parameters, per_age_data,
							agenums, n_sim_steps, modelInputParameters_.seedlist,
							modelInputParameters_.day_shut, obs_selections.hospitalised, obs_selections.deaths);

                if (acceptedParticleCount < modelInputParameters_.nParticalLimit) {
                    //count the number of simulations that were used to reach the maximum number of accepted particles
                    ++nsim_count;

                    // Record the particle if it passes the tolerance criteria
                    if (ParticlePassesTolerances(outs_vec, smc)) {				
                        //#pragma omp critical
                        {
                            if (smc == 0)
                            {
                                outs_vec.weight = 1;
                            }
                            else
                            {
                                ComputeParticleWeight(previousParticles, outs_vec, vlimitKernel);
                            }
                            currentParticles.push_back(outs_vec);
                            ++acceptedParticleCount;
                            
                            if (acceptedParticleCount % 10 == 0) (*log_) << "|" << std::flush;
                        }
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
		(*log_) << "\nStep:" << smc
			<< ", <number of accepted particles> " << prevAcceptedParticleCount
			<< "; <number of simulations> " << nsim_count
			<< "; <computation time> " <<  time_taken
			<< " seconds.\n";

		//break the ABC-smc at the step where no particles were accepted
		if (prevAcceptedParticleCount > 0) {
			IO::WriteOutputsToFiles(smc, modelInputParameters_.herd_id, prevAcceptedParticleCount,
				nInferenceParams, currentParticles, outDir_, log_);
		}
	}

	//output on screen the overall computation time
	(*log_) << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
}

bool InferenceFramework::ParticlePassesTolerances(const particle& p, int smc) {
    return p.nsse_cases <= toleranceLimits_[smc] && p.nsse_deaths <= toleranceLimits_[smc];
}

void InferenceFramework::ModelSelect(EERAModel::particle& outvec, const std::vector<params>& fixed_parameters,
	const AgeGroupData& per_age_data, std::vector <int> agenums, const int& n_sim_steps, 
	seed seedlist, int day_shut, const std::vector<int>& obsHosp, 
	const std::vector<int>& obsDeaths) {

	//---------------------------------------
	// the root model
	//---------------------------------------
	
	Status status = model_->Run(outvec.parameter_set, fixed_parameters, per_age_data, seedlist, day_shut,
							agenums, n_sim_steps);

	//---------------------------------------
	// compute the  sum of squared errors for daily observations
	//---------------------------------------
	double sum_sq_cases = Utilities::sse_calc<int>(status.simulation, obsHosp);

	//---------------------------------------
	// compute the  sum of squared errors for weekly observations
	//---------------------------------------	
	// aggregate daily values to weekly values 
	// due to reporting process creating uncertaincies within weeks

	const int week_length = 7;

	const std::vector<int> obs_death_red = Utilities::AccumulateEveryN(obsDeaths, week_length);
	const std::vector<int> sim_hospital_death_red = Utilities::AccumulateEveryN(status.hospital_deaths, week_length);
	
	double sum_sq_deaths = Utilities::sse_calc<int>(sim_hospital_death_red, obs_death_red);

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
	outvec.hospital_death_outs = status.hospital_deaths;
	outvec.death_outs = status.deaths;
	outvec.end_comps = Model::compartments_to_vector(status.ends);
}

void InferenceFramework::ComputeParticleWeight(std::vector<EERAModel::particle> pastPart,
	EERAModel::particle &currentPart, const std::vector<double>& vlimitKernel) {

    int pastNpart = pastPart.size();
    unsigned int nPar = vlimitKernel.size();
    double denom = 0;
    for (int jj = 0; jj < pastNpart; ++jj) {
        int rval_dist=0;
        double m = 0.0;
        //for each parameter of each particle
        for (int ii = 0; ii < nPar; ++ii) {
            //compute the distance between the jjth previous particle and the chosen particle for the iith parameters
            m= std::fabs(currentPart.parameter_set[ii] - pastPart[jj].parameter_set[ii]);
            //add 1 to rval_dist if iith parameter is within the range of the corresponding parameters of the jjth particle
            if(m <= vlimitKernel[ii]){
                ++rval_dist;
            }
        }
        //if all values of parameter of the particle are within the range of perturbation of all parameters of a given source particle
        if(rval_dist==nPar){
            //add the weight of this particle to denom
            denom += pastPart[jj].weight;
        }
    }

    if(denom == 0) denom = 1.0 ;
    //the weight of each particle is 1 over the sum of the weights from the particles at s-1 that may generate this particle
    currentPart.weight = 1.0 / denom;
}

void ComputeKernelWindow(int nPar, const std::vector<particle>& particleList, double kernelFactor,
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

std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList) {
	
	std::vector<double> weight_val;
	for (auto p : particleList) {
		weight_val.push_back(p.weight);
	}
	
	return std::discrete_distribution<int>(weight_val.begin(), weight_val.end());
}

} // namespace Inference
} // namespace EERAModel
