#pragma once

#include <vector>
#include <string>

namespace EERAModel {

/**
 * @brief Structure containing counters for different population categories
 * 
 * Integers to count the number of people within the different compartments within the model 
 */
struct Compartments
{
	int S = 0;    	/*!< Number of susceptible individuals (not infected). */	
	int E = 0;		/*!< Number of infected individuals but not yet infectious (exposed). */	
	int E_t = 0;	/*!< Number of exposed individuals and tested positive. */	
	int I_p = 0;	/*!< Number of infected and infectious symptomatic individuals but at pre-clinical stage (show yet no symptoms). */	
	int I_t = 0;	/*!< Number of tested positive individuals that infectious. */	
	int I1 = 0;		/*!< Number of infected and infectious asymptomatic individuals: first stage. */	
	int I2 = 0;		/*!< Number of infected and infectious asymptomatic individuals: second stage.  */	
	int I3 = 0;		/*!< Number of infected and infectious asymptomatic individuals: third stage.  */	
	int I4 = 0;		/*!< Number of infected and infectious asymptomatic individuals: last stage.  */	
	int I_s1 = 0;	/*!< Number of infected and infectious symptomatic individuals: first stage. */	
	int I_s2 = 0;	/*!< Number of infected and infectious symptomatic individuals: second stage.  */	
	int I_s3 = 0;	/*!< Number of infected and infectious symptomatic individuals: thrid stage. */	
	int I_s4 = 0;	/*!< Number of infected and infectious symptomatic individuals: last stage.  */	
	int H = 0;		/*!< Number of infected individuals that are hospitalised.  */	
	int R = 0;		/*!< Number of infected individuals that are recovered from infection.   */	
	int D = 0;		/*!< Number of dead individuals due to disease. */	
};

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
	std::vector<int> hospital_death_outs;
	std::vector<int> death_outs;
	std::vector< Compartments > end_comps;
};

/**
 * @brief Structure for fixed parameters
 */
struct params {
	double T_lat;
	double juvp_s;
	double T_inf;
	double T_rec;
	double T_sym;
	double T_hos;
	int K;
	double inf_asym;
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
	bool use_fixed_seed;
	unsigned long int seed_value;
};

/**
 * @brief Enumeration identifying which model structure the code should use
 */
enum class ModelStructureId
{
    UNKNOWN,
    ORIGINAL,
    IRISH,
	TEMP
};

/**
 * @brief Enumeration identifying which mode to run the model
 */
enum class ModelModeId
{
	UNKNOWN,
    INFERENCE,
	PREDICTION
};

/**
 * @brief Parameters used for validation of observation data
 */
struct ValidationParameters {
	int nHealthBoards;
	int nAgeGroups;
	int nCfrCategories;
	int nCasesDays;
};

/**
 * @brief Observations for Inference framework
 */
struct ObservationsForInference {
	std::vector<std::vector<int>> cases;
	std::vector<std::vector<int>> deaths;
};

/**
 * @brief Observations for Models
 */
struct ObservationsForModels {
	std::vector<std::vector<int>> cases;
	std::vector<std::vector<double>> age_pop;
	std::vector<std::vector<double>> waifw_norm;
	std::vector<std::vector<double>> waifw_home;
	std::vector<std::vector<double>> waifw_sdist;
	std::vector<std::vector<double>> cfr_byage;
	std::vector<std::vector<double> > pf_pop;
};

/**
 * @brief Structure to hold per-age group data vectors
 * 
 * Includes probabilities and mean number of contacts per age group information
 */
struct AgeGroupData
{
	
	std::vector<std::vector<double>> waifw_norm;	/*!< mean number of daily contacts between age groups (overall). */	
	std::vector<std::vector<double>> waifw_home;	/*!< mean number of daily contacts between age groups (home only). */
	std::vector<std::vector<double>> waifw_sdist;	/*!< mean number of daily contacts between age groups (not school, not work). */
	std::vector<std::vector<double>> cfr_byage;		/*!< Case fatality ratio by age. */
	std::vector<double> pf_byage;					/*!< Frailty Probability by age. */
};

/**
 * @brief Structure to hold status objects
 * 
 * Contains vectors containing information on the trajectories of the predicted disease status for individuals involved in the model 
 * (i.e. number of incident cases, number of incident deaths for each simulation step).
 */
struct Status
{
	std::vector<int> simulation;				/*!< Number of incident (newly detected/reported) cases of covid in each simulation step. */
	std::vector<int> deaths;					/*!< Overall number of incident deaths due to covid in each simulation step. */
	std::vector<int> hospital_deaths;			/*!< Number of incident deaths reported at hospital due to covid in each simulation step. */
	std::vector<Compartments> ends;				/*!< Population per epidemiological state and per age group on last day. */
	std::vector<std::vector<Compartments>> 
		pop_array;								/*!< Population per epidemiological state and per age group on every day. */
};

/**
 * @brief Structure containing population counters after infection
 * 
 * Contains counters for the number of newly detected (due to testing or hospitalisation) cases, deaths and hospitalisations
 * at each time step.
 */
struct InfectionState
{
	int detected = 0;				/*!< Number of newly detected (due to testing or hospitalisation) cases at each time step. */
	int hospitalised = 0;			/*!< Number of newly detected cases due to their hospitalisation at each time step. */
	int deaths = 0;					/*!< Overall number of incident deaths due to covid at each time step.  */
	int hospital_deaths = 0;		/*!< Number of incident deaths reported at hospital due to covid at each time step. */
};

/**
 * @class KernelWindow
 * @brief Struct representing the kernel window for a given inference parameter
 */
struct KernelWindow 
{
    double kernel;
    double max;
    double min; 
}; 

/**
* @brief Supplementary input parameters used by main.cpp
*/
struct SupplementaryInputParameters
{
	seed seedlist;
	ModelStructureId model_structure = ModelStructureId::UNKNOWN;
	ModelModeId run_type = ModelModeId::UNKNOWN;
};

/**
 * @brief Common model input parameters used by all models
 */
struct CommonModelInputParameters 
{
	params paramlist;
	int herd_id;
	int totN_hcw;
};

/**
 * @brief Configuration for an inference config run
 */
struct InferenceConfig
{
	seed seedlist;
	params paramlist;
	int herd_id;
	int day_shut;
	double tau;
	double prior_pinf_shape1;
	double prior_pinf_shape2;
	double prior_phcw_shape1;
	double prior_phcw_shape2;
	double prior_chcw_mean;
	double prior_d_shape1;
	double prior_d_shape2;
	double prior_q_shape1;
	double prior_q_shape2;
	double prior_rrd_shape1;
	double prior_rrd_shape2;
	double prior_lambda_shape1;
	double prior_lambda_shape2;
	double prior_ps_shape1;
	double prior_ps_shape2;
	int nsteps;
	double kernelFactor;
	int nSim;
	int nParticleLimit;
	std::vector<double> toleranceLimit;
	ObservationsForInference observations;
};

/**
 * @brief Configuration for a prediction framework run
 */
struct PredictionConfig
{
    seed seedlist;
    int day_shut;
    int n_sim_steps;    /*! Number of steps over which to run the model */
    int n_iterations;   /*! Number of model runs to perform */
    int index;          /*! Index of selected parameters within posterior parameters file */
    std::vector<double> posterior_parameters;  /*! Set of model parameters */
    params fixedParameters; /*! Set of model fixed parameters */
};

} // namespace EERAModel
