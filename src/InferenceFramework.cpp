#include "InferenceFramework.h"
#include "Observations.h"
#include "IO.h"
#include "ModelCommon.h"
#include <functional>
#include <algorithm>
#include <cmath>

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

    inferenceParticleGenerator_ = std::make_shared<InferenceParticleGenerator>(
        static_cast<unsigned int>(modelInputParameters_.nPar), modelInputParameters.kernelFactor, rng);
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

    logSettings(modelInputParameters_, log_);
    logFixedParameters(modelInputParameters_, log_);

	//keep information for the health board if interest
	const std::vector<double> pf_byage = observations_.pf_pop[modelInputParameters_.herd_id - 1];//define frailty structure of the shb of interest.

	const AgeGroupData per_age_data = {observations_.waifw_norm, observations_.waifw_home, observations_.waifw_sdist, 
										observations_.cfr_byage, pf_byage};

	std::vector<params> fixed_parameters = Model::BuildFixedParameters(
        observations_.waifw_norm.size(), modelInputParameters_.paramlist);
    
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

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	(*log_) << "[Simulations]:\n";
	for (int smc = 0; smc < modelInputParameters_.nsteps; ++smc) {//todo:

		//the abort statement for keeping the number of particles less than 1000
		bool aborting = false;

		//initialise the counter for the number of simulation prior maximum particle is reached
		int nsim_count = 0;

        if (smc > 0) {
			previousParticles = currentParticles;
			currentParticles.clear();

            inferenceParticleGenerator_->Update(previousParticles);
		}

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
				if (currentParticles.size() >= modelInputParameters_.nParticalLimit) {
					aborting = true;
					//#pragma omp flush (aborting)
				}

				// Generate a new inference particle
				particle par = inferenceParticleGenerator_->GenerateNew(smc, sim,
                    inferenceParameterGenerator_, previousParticles);

				//run the model and compute the different measures for each potential parameters value
				ModelSelect(par, fixed_parameters, per_age_data,
							agenums, n_sim_steps, modelInputParameters_.seedlist,
							modelInputParameters_.day_shut, obs_selections.hospitalised, obs_selections.deaths);

                if (currentParticles.size() < modelInputParameters_.nParticalLimit) {
                    //count the number of simulations that were used to reach the maximum number of accepted particles
                    ++nsim_count;

                    // Record the particle if it passes the tolerance criteria
                    if (ParticlePassesTolerances(par, smc)) {				
                        //#pragma omp critical
                        {
                            par.weight = ComputeParticleWeight(smc, previousParticles,
                                par, inferenceParticleGenerator_->KernelWindows());

                            currentParticles.push_back(par);
                            
                            if (currentParticles.size() % 10 == 0) (*log_) << "|" << std::flush;
                        }
                    }
                }			
			}
		}

		//time taken per step
		double time_taken;
		if(smc == 0){
			time_taken = double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		} else {
			time_taken = double( clock() - time_taken1 ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		}

		(*log_) << "\nStep:" << smc
			<< ", <number of accepted particles> " << currentParticles.size()
			<< "; <number of simulations> " << nsim_count
			<< "; <computation time> " <<  time_taken
			<< " seconds.\n";

		if (currentParticles.size() > 0) {
			IO::WriteOutputsToFiles(smc, modelInputParameters_.herd_id, currentParticles.size(),
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

double InferenceFramework::ComputeParticleWeight(int smc, const std::vector<EERAModel::particle>& previousParticles,
	const particle& currentParticle, const std::vector<KernelWindow> kernelWindows) {

    double weight = 1.0;
       
    if (smc > 0) 
    {
        double denom = 0;
        
        for (const auto& previousParticle : previousParticles) 
        {
            if (ParticlesAreClose(currentParticle, previousParticle, kernelWindows))
                denom += previousParticle.weight;
        }

        if (denom == 0) 
            denom = 1.0;
        
        weight = 1.0 / denom;
    }

    return weight;
}

bool InferenceFramework::ParticlesAreClose(const particle& first, const particle& second,
    const std::vector<KernelWindow> kernelWindows) {
       
    bool close = true;
    for (unsigned int i = 0; i < first.parameter_set.size(); ++i) 
    {
        double distance = std::fabs(first.parameter_set[i] - second.parameter_set[i]);
        
        if (distance > kernelWindows[i].kernel) 
            close = false;
    }

    return close;
}

std::vector<KernelWindow> ComputeKernelWindow(int nPar, const std::vector<particle>& particles, double kernelFactor) {

    std::vector<KernelWindow> kernelWindows(nPar);

	for (int i = 0; i < nPar; ++i) 
    {	
		std::function<bool(particle, particle)> compare = 
			[&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; };

		particle max = *std::max_element(particles.begin(), particles.end(), compare);
		
		particle min = *std::min_element(particles.begin(), particles.end(), compare);

		kernelWindows[i].max = max.parameter_set[i];
		kernelWindows[i].min = min.parameter_set[i];
		
		kernelWindows[i].kernel = kernelFactor * std::fabs(kernelWindows[i].max - kernelWindows[i].min);
	}

    return kernelWindows;
}

std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList) {
	
	std::vector<double> weight_val;
	for (auto p : particleList) {
		weight_val.push_back(p.weight);
	}
	
	return std::discrete_distribution<int>(weight_val.begin(), weight_val.end());
}

void logSettings(const ModelInputParameters& params, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Settings]:\n";
	(*log) << "    number of parameters tested: "<< params.nPar << std::endl;
    (*log) << "    seeding method: "<< params.seedlist.seedmethod<<  std::endl;
	if (params.seedlist.seedmethod == "random"){
		(*log) << "    number of seed: " << params.seedlist.nseed << std::endl;
	} else if(params.seedlist.seedmethod == "background"){
		(*log) << "    duration of the high risk period (hrp): " << params.seedlist.hrp << std::endl;
	}
    (*log) << "    model structure: " << 
        ((params.model_structure == ModelStructureId::ORIGINAL) ? "Original" : "Irish") << std::endl;
}

void logFixedParameters(const ModelInputParameters& params, Utilities::logging_stream::Sptr log)
{
    (*log) << "[Fixed parameter values]:\n";
	(*log) << "    latent period (theta_l): " << params.paramlist.T_lat <<std::endl;
	(*log) << "    pre-clinical period (theta_i): " << params.paramlist.T_inf <<std::endl;
	(*log) << "    asymptomatic period (theta_r): " << params.paramlist.T_rec <<std::endl;
	(*log) << "    symptomatic period (theta_s): " << params.paramlist.T_sym <<std::endl;
	(*log) << "    hospitalisation stay (theta_h): " << params.paramlist.T_hos <<std::endl;
	(*log) << "    pre-adult probability of symptoms devt (p_s[0]): " << params.paramlist.juvp_s <<std::endl;
	(*log) << "    bed capacity at hospital (K): " << params.paramlist.K <<std::endl;
	(*log) << "    relative infectiousness of asymptomatic (u): " << params.paramlist.inf_asym <<std::endl;
}

InferenceParticleGenerator::InferenceParticleGenerator(unsigned int nInferenceParams, double kernelFactor,
        Random::RNGInterface::Sptr rng)
         : nInferenceParams_(nInferenceParams),
           kernelFactor_(kernelFactor),
           rng_(rng) {}

void InferenceParticleGenerator::Update(std::vector<particle> particles)
{
    kernelWindows_ = ComputeKernelWindow(nInferenceParams_, particles, kernelFactor_);
    weightDistribution_ = ComputeWeightDistribution(particles);
}

particle InferenceParticleGenerator::GenerateNew(int smc, int sim,
    InferenceParameterGenerator::Sptr parameterGenerator, const std::vector<particle>& previousParticles) 
{
    particle par;
    par.iter = sim;	

    if (smc > 0) {
        int pick_val = weightDistribution_(rng_->MT19937());
        par.parameter_set = parameterGenerator->GenerateWeighted(
            previousParticles[pick_val].parameter_set, kernelWindows_);
    } else {
        par.parameter_set = parameterGenerator->GenerateInitial();
    }

    return par;
}

} // namespace Inference
} // namespace EERAModel
