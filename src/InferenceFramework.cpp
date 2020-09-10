#include "InferenceFramework.h"
#include "Observations.h"
#include "IO.h"
#include "ModelCommon.h"
#include "Timer.h"
#include <functional>
#include <algorithm>
#include <cmath>

namespace EERAModel {
namespace Inference {

InferenceFramework::InferenceFramework(Model::ModelInterface::Sptr model,
    const InferenceConfig& inferenceConfig,
    Random::RNGInterface::Sptr rng,
    Utilities::logging_stream::Sptr log,
	IO::IOdatapipeline *datapipeline)
    : model_(model),
      inferenceConfig_(inferenceConfig),
      rng_(rng),
      log_(log),
      datapipeline_(datapipeline) {

    const unsigned int nInferenceParams = Model::ModelParameters::NPARAMS;

    std::vector<double> flag1 = {
		inferenceConfig.prior_pinf_shape1,
	 	inferenceConfig.prior_phcw_shape1,
		inferenceConfig.prior_chcw_mean, 
		inferenceConfig.prior_d_shape1, 
		inferenceConfig.prior_q_shape1,
		inferenceConfig.prior_ps_shape1,
		inferenceConfig.prior_rrd_shape1,		
		inferenceConfig.prior_lambda_shape1
	};	
	std::vector<double> flag2 = {
		inferenceConfig.prior_pinf_shape2, 
		inferenceConfig.prior_phcw_shape2, 
		inferenceConfig.prior_chcw_mean, 
		inferenceConfig.prior_d_shape2, 
		inferenceConfig.prior_q_shape2,
		inferenceConfig.prior_ps_shape2,
		inferenceConfig.prior_rrd_shape2,		
		inferenceConfig.prior_lambda_shape2
	};

    inferenceParameterGenerator_ = std::make_shared<InferenceParameterGenerator>(rng, flag1, flag2);

    inferenceParticleGenerator_ = std::make_shared<InferenceParticleGenerator>(
        nInferenceParams, inferenceConfig_.kernelFactor, rng);
}

int InferenceFramework::GetTimeOffSet(const InferenceConfig& inferenceConfig)
{
	int time_back = 0;
	if(inferenceConfig.seedlist.seedmethod == "background") {
		time_back = inferenceConfig.seedlist.hrp;
	} else {
		time_back = static_cast<int>(inferenceConfig.paramlist.T_inf + inferenceConfig.paramlist.T_sym);
		
	}

	return time_back;
}

void InferenceFramework::Run()
{
    const unsigned int nInferenceParams = Model::ModelParameters::NPARAMS;
    
    SimpleTimer runTimer;

	const int time_back = GetTimeOffSet(inferenceConfig_);
	const std::vector<int>& regionalCases = inferenceConfig_.observations.cases[inferenceConfig_.herd_id];
	const std::vector<int>& regionalDeaths = inferenceConfig_.observations.deaths[inferenceConfig_.herd_id];
	const std::vector<int>& timeStamps = inferenceConfig_.observations.cases[0];

	const Observations::ObsSelect obs_selections = Observations::SelectObservations(
		inferenceConfig_.day_shut, timeStamps, regionalCases,
		regionalDeaths, time_back, log_);
	
	inferenceConfig_.seedlist.day_intro = obs_selections.sim_time.day_intro;
	
    (*log_) << "[Health Board settings]:\n";
	(*log_) << "    SHB id: " << inferenceConfig_.herd_id <<'\n';
	(*log_) << "    Simulation period: " << obs_selections.sim_time.duration << "days\n";
	(*log_) << "    time step: " << inferenceConfig_.tau << "days\n";

	const int n_sim_steps = static_cast<int>(ceil(obs_selections.sim_time.duration/inferenceConfig_.tau));

	// currentParticles is the list of particles accepted on the current abc-smc loop
    // previousParticles is the list accepted on the previous loop
	std::vector<particle> currentParticles, previousParticles;

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	(*log_) << "[Simulations]:\n";
	for (int smc = 0; smc < inferenceConfig_.nsteps; ++smc) {

		SimpleTimer stepTimer;
       
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
		for (int sim = 0; sim < inferenceConfig_.nSim; ++sim) {//todo
			if (!aborting) {
				// Abort if we've reached the limit on the number of accepted particles
				if (currentParticles.size() >= inferenceConfig_.nParticleLimit) {
					aborting = true;
				}

				// Generate a new inference particle
				particle par = inferenceParticleGenerator_->GenerateNew(smc, sim,
                    inferenceParameterGenerator_, previousParticles);

				// Run the model and compute the different measures for each potential parameters value
				ModelSelect(par, n_sim_steps, inferenceConfig_.seedlist,
							inferenceConfig_.day_shut, obs_selections.hospitalised, obs_selections.deaths);

                // Only record the particle if we haven't reached our limit for this run
                if (currentParticles.size() < inferenceConfig_.nParticleLimit) {
                    //count the number of simulations that were used to reach the maximum number of accepted particles
                    ++nsim_count;

                    // Record the particle if it passes the tolerance criteria
                    if (ParticlePassesTolerances(par, smc)) {				
                        par.weight = ComputeParticleWeight(smc, previousParticles,
                            par, inferenceParticleGenerator_->KernelWindows());

                        currentParticles.push_back(par);
                        
                        // Output a marker for every 10 particles accepted
                        if (currentParticles.size() % 10 == 0) (*log_) << "|" << std::flush;
                    }
                }
			}
		}

		(*log_) << "\nStep:" << smc
			<< ", <number of accepted particles> " << currentParticles.size()
			<< "; <number of simulations> " << nsim_count
			<< "; <computation time> " <<  stepTimer.elapsedTime()
			<< " seconds.\n";

		// Record the list of accepted particles in the output files.
		datapipeline_->WriteOutputsToFiles(smc, inferenceConfig_.herd_id, currentParticles.size(),
            nInferenceParams, currentParticles, model_->ModelName());
	}

	// Output on screen the overall computation time
	(*log_) << runTimer.elapsedTime() << " seconds." << std::endl;

	datapipeline_->WriteLog("inference", model_->ModelName());
}

bool InferenceFramework::ParticlePassesTolerances(const particle& p, int smc) 
{
    const std::vector<double>& toleranceLimit = inferenceConfig_.toleranceLimit;
    
    return p.nsse_cases <= toleranceLimit[smc] && p.nsse_deaths <= toleranceLimit[smc];
}

void InferenceFramework::ModelSelect(EERAModel::particle& outvec, int n_sim_steps, 
	seed seedlist, int day_shut, const std::vector<int>& obsHosp, const std::vector<int>& obsDeaths) {

	//---------------------------------------
	// the root model
	//---------------------------------------
	
	Status status = model_->Run(outvec.parameter_set, seedlist, day_shut, n_sim_steps);

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

	const std::vector<int> obs_death_red = Utilities::AccumulateEveryN<int>(obsDeaths, week_length);
	const std::vector<int> sim_hospital_death_red = Utilities::AccumulateEveryN<int>(status.hospital_deaths, week_length);
	
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
	outvec.end_comps = status.ends;
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

	for (int i = 0; i < nPar; ++i) {	
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
	
	std::vector<double> weight_val(particleList.size(), 0.0);

	for (unsigned int i = 0; i < weight_val.size(); i++) {
		weight_val[i] = particleList[i].weight;
	}
	
	return std::discrete_distribution<int>(weight_val.begin(), weight_val.end());
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
