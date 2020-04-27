#pragma once

#include <vector>

#include "ModelTypes.h"
#include "IniFile.h"

namespace EERAModel {

/**
 * @brief Structure containing the inputs to a model run
 */
struct ModelInputParameters
{
	int herd_id;
	double tau;
	int num_threads;
	bool apply_death_tolerance_limit;
	int nsteps;
	int nParticalLimit;
	int nSim;
	double kernelFactor;
	std::vector<double> toleranceLimit;
	params paramlist;
	seed seedlist;
	int day_shut;
	int totN_hcw;
	int nPar;
	double prior_pinf_shape1;
	double prior_pinf_shape2;
	double prior_phcw_shape1;
	double prior_phcw_shape2;
	double prior_chcw_mean;
	double prior_d_shape1;
	double prior_d_shape2;
	double prior_q_shape1;
	double prior_q_shape2;
	double prior_rrdh_shape1;
	double prior_rrdh_shape2;
	double prior_rrdc_shape1;
	double prior_rrdc_shape2;
	double prior_rrh_shape1;
	double prior_rrh_shape2;
	double prior_lambda_shape1;
	double prior_lambda_shape2;
};

/**
 * @brief Read model input parameters from an INI file
 * 
 * @param filePath Path to the INI file
 * 
 * @return Model parameters
 */
ModelInputParameters ReadParametersFromFile(const std::string& filePath)
{
	ModelInputParameters modelInputParameters;

	CIniFile parameters;

	//Settings
	modelInputParameters.herd_id = atoi(parameters.GetValue("shb_id", "Settings", filePath).c_str());

	//time step for the modelling process
	modelInputParameters.tau = atof(parameters.GetValue("tau", "Settings", filePath).c_str());

	//number of threads used for computation
	modelInputParameters.num_threads = atoi(parameters.GetValue("num_threads", "Settings", filePath).c_str());

	// Whether or not to apply the death tolerance limit
	modelInputParameters.apply_death_tolerance_limit = static_cast<bool>(
		atoi(parameters.GetValue("apply_death_tolerance_limit", "Settings", filePath).c_str())
	);

	//Seed settings
	modelInputParameters.seedlist.seedmethod = parameters.GetValue("seedmethod", "Seed settings", filePath).c_str();
	if(modelInputParameters.seedlist.seedmethod == "random"){
		modelInputParameters.seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", filePath).c_str());
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		modelInputParameters.seedlist.hrp = atoi(parameters.GetValue("hrp", "Seed settings", filePath).c_str());
	} else {
		cout<< "Warning!!! Unknown method - using random seed method instead." << endl;
		modelInputParameters.seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", filePath).c_str());
	}
	
	//Fit settings
	modelInputParameters.nsteps = atoi(parameters.GetValue("nsteps", "Fit settings", filePath).c_str());
	modelInputParameters.nParticalLimit = atoi(parameters.GetValue("nParticLimit", "Fit settings", filePath).c_str());
	modelInputParameters.nSim = atoi(parameters.GetValue("nSim", "Fit settings", filePath).c_str());
	modelInputParameters.kernelFactor = atof(parameters.GetValue("kernelFactor", "Fit settings", filePath).c_str());

	//Tolerance settings
	for (int ii = 1; ii <= modelInputParameters.nsteps; ++ii) {
		modelInputParameters.toleranceLimit.push_back(0.0);
	}

	for (int ii = 1; ii <= modelInputParameters.nsteps; ++ii) {
		if(ii==1) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key1", "Tolerance settings", filePath).c_str());
		if(ii==2) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key2", "Tolerance settings", filePath).c_str());
		if(ii==3) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key3", "Tolerance settings", filePath).c_str());
		if(ii==4) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key4", "Tolerance settings", filePath).c_str());
		if(ii==5) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key5", "Tolerance settings", filePath).c_str());
		if(ii==6) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key6", "Tolerance settings", filePath).c_str());
		if(ii==7) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key7", "Tolerance settings", filePath).c_str());
		if(ii==8) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key8", "Tolerance settings", filePath).c_str());
		if(ii==9) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key9", "Tolerance settings", filePath).c_str());
		if(ii==10) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key10", "Tolerance settings", filePath).c_str());
	}

	//Fixed parameters
	modelInputParameters.paramlist.T_lat = atof(parameters.GetValue("T_lat", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.p_s = atof(parameters.GetValue("p_s", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.juvp_s = atof(parameters.GetValue("juvp_s", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_inf = atof(parameters.GetValue("T_inf", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_rec = atof(parameters.GetValue("T_rec", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_sym = atof(parameters.GetValue("T_sym", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_hos = atof(parameters.GetValue("T_hos", "Fixed parameters", filePath).c_str());
	modelInputParameters.day_shut = atoi(parameters.GetValue("day_shut", "Fixed parameters", filePath).c_str());
	modelInputParameters.totN_hcw = atoi(parameters.GetValue("totN_hcw", "Fixed parameters", filePath).c_str());

	//priors settings
	modelInputParameters.nPar = atoi(parameters.GetValue("nPar", "Priors settings", filePath).c_str());
	modelInputParameters.prior_pinf_shape1 = atof(parameters.GetValue("prior_pinf_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_pinf_shape2 = atof(parameters.GetValue("prior_pinf_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_phcw_shape1 = atof(parameters.GetValue("prior_phcw_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_phcw_shape2 = atof(parameters.GetValue("prior_phcw_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_chcw_mean = atof(parameters.GetValue("prior_chcw_mean", "Priors settings", filePath).c_str());
	modelInputParameters.prior_d_shape1 = atof(parameters.GetValue("prior_d_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_d_shape2 = atof(parameters.GetValue("prior_d_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_q_shape1 = atof(parameters.GetValue("prior_q_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_q_shape2 = atof(parameters.GetValue("prior_q_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrdh_shape1 = atof(parameters.GetValue("prior_rrdh_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrdh_shape2=atof(parameters.GetValue("prior_rrdh_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrdc_shape1=atof(parameters.GetValue("prior_rrdc_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrdc_shape2=atof(parameters.GetValue("prior_rrdc_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrh_shape1=atof(parameters.GetValue("prior_rrh_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrh_shape2=atof(parameters.GetValue("prior_rrh_shape2", "Priors settings", filePath).c_str());	
	
	modelInputParameters.prior_lambda_shape1 = atof(parameters.GetValue("prior_lambda_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_lambda_shape2 = atof(parameters.GetValue("prior_lambda_shape2", "Priors settings", filePath).c_str());

	return modelInputParameters;
}

} // namespace EERAModel