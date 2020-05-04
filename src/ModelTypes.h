#pragma once

#include <vector>
#include <string>

namespace EERAModel {

/**
 * @brief Structure for particles generated during inference
 */
struct particle {
	double nsse_cases;
	double nsse_deaths;
	std::vector<double> parameter_set;
	int iter;
	double weight;
	std::vector<int> simu_outs;
	std::vector<int> death_outs;
	std::vector< std::vector<int> > end_comps;
};

/**
 * @brief Structure for fixed parameters
 */
struct params {
	double T_lat;
	double p_s;
	double juvp_s;
	double T_inf;
	double T_rec;
	double T_sym;
	double T_hos;
};

/**
 * @brief Structure for seed settings
 */
struct seed {
	std::string seedmethod;
	int nseed;
	int hrp;
	int day_intro;
	double lambda;
};

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
 * @brief Model input observations
 */
struct Observations {
	std::vector<std::vector<int>> cases;
	std::vector<std::vector<int>> deaths;
	std::vector<std::vector<double>> age_pop;
	std::vector<std::vector<double>> waifw_norm;
	std::vector<std::vector<double>> waifw_home;
	std::vector<std::vector<double>> waifw_sdist;
	std::vector<std::vector<double>> cfr_byage;
	std::vector<std::vector<double> > pf_pop;
};

} // namespace EERAModel
