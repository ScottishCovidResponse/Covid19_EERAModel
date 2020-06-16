#include "PredictionFramework.h"
#include "ModelCommon.h"
#include "Observations.h"
#include "IO.h"

namespace EERAModel {
namespace Prediction {

PredictionFramework::PredictionFramework(
    Model::ModelInterface::Sptr model,
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
       log_(log)
{
    fixedParameters_ = Model::BuildFixedParameters(
        observations.waifw_norm.size(), modelInputParameters.paramlist
    );
        
    ageGroupData_ = AgeGroupData{
            observations.waifw_norm,
            observations.waifw_home,
            observations.waifw_sdist,
            observations.cfr_byage,
            observations.pf_pop[modelInputParameters.herd_id - 1]
    };
    
    regionalPopulation_ = Model::GetPopulationOfRegion(
        observations, modelInputParameters.herd_id
    );

    healthCareWorkers_ = Model::ComputeNumberOfHCWInRegion(
        regionalPopulation_, modelInputParameters.totN_hcw, observations
    );
    
    ageNums_ = Model::ComputeAgeNums(
        modelInputParameters.herd_id, regionalPopulation_, healthCareWorkers_, observations
    );
}

void PredictionFramework::Run(std::vector<double> parameterSet, int nSimulationSteps)
{
    clock_t startTime = clock();

    (*log_) << "[Settings]:\n";
    (*log_) << "    number of parameters tested: "<< modelInputParameters_.nPar << std::endl;
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

    Status status = model_->Run(parameterSet, fixedParameters_, ageGroupData_, modelInputParameters_.seedlist,
        modelInputParameters_.day_shut, ageNums_, nSimulationSteps);

    double time_taken = double(clock() - startTime)/(double)CLOCKS_PER_SEC;

    (*log_) << "\n <computation time> " << time_taken << " seconds.\n";

    std::vector< std::vector<int> > end_comps;
	end_comps = compartments_to_vector(status.ends);

    IO::WritePredictionsToFiles(status, end_comps, outDir_, log_);
}

std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec)
{
	std::vector<std::vector<int>> _temp;

	for(auto cmps : cmps_vec)
	{
		_temp.push_back({cmps.S, cmps.E, cmps.E_t, cmps.I_p,
						cmps.I_t, cmps.I1, cmps.I2, cmps.I3,
						cmps.I4, cmps.I_s1, cmps.I_s2, cmps.I_s3,
						cmps.I_s4, cmps.H, cmps.R, cmps.D});
	}

	return _temp;
}

} // namespace Prediction
} // namespace EERAModel