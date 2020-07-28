#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>

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
	std::vector<int> hospital_death_outs;
	std::vector<int> death_outs;
	std::vector< std::vector<int> > end_comps;
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
    int index;          /*! Index of selected parameters within posterior parameters file */
    std::vector<double> posterior_parameters;  /*! Set of model parameters */
    params fixedParameters; /*! Set of model fixed parameters */
};

class HealthBoardData
{
	public:
		HealthBoardData(
			const std::string& HealthBoardLabel,
			const std::vector<std::vector<int>>& cases,
			const std::vector<std::vector<int>>& deaths,
			const std::vector<std::vector<double>>& age_pop,
			const std::vector<std::vector<double>>& waifw_norm, 
			const std::vector<std::vector<double>>& waifw_home, 
			const std::vector<std::vector<double>>& waifw_sdist, 
			const std::vector<std::vector<double>>& cfr_byage, 
			const std::vector<std::vector<double>>& pf_byage
		):
			health_board_label_(HealthBoardLabel),
			cases_(cases),
			deaths_(deaths),
			age_pop_(age_pop),
			waifw_norm_(waifw_norm),
			waifw_home_(waifw_home),
			waifw_sdist_(waifw_sdist),
			cfr_byage_(cfr_byage),
			pf_byage_(pf_byage) {};

		std::vector<int> DailyCases = cases_[HealthBoardMap_[health_board_label_]];
		std::vector<int> DailyDeaths = deaths_[HealthBoardMap_[health_board_label_]];

		unsigned int PopulationSize = cases_[HealthBoardMap_[health_board_label_]][0];
		double PopulationProportion_Pre_20 = 
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_Pre_20"]][HealthBoardMap_[health_board_label_]];
		double PopulationProportion_20_29 = 
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_20_29"]][HealthBoardMap_[health_board_label_]];
		double PopulationProportion_30_39 =
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_30_39"]][HealthBoardMap_[health_board_label_]];
		double PopulationProportion_40_49 =
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_40_49"]][HealthBoardMap_[health_board_label_]];;
		double PopulationProportion_50_59 =
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_50_59"]][HealthBoardMap_[health_board_label_]];;
		double PopulationProportion_60_69 =
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_60_69"]][HealthBoardMap_[health_board_label_]];;
		double PopulationProportion_Post_70 =
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_Post_70"]][HealthBoardMap_[health_board_label_]];;
		double PopulationProportion_HealthCareWorkers = 
			PopulationSize*age_pop_[AgeGroupMap_["AgeGroup_HealthCareWorkers"]][HealthBoardMap_[health_board_label_]];;

		double FullAgeContacts_ToFrom(const std::string& To, const std::string& From) {
			return waifw_norm_[AgeGroupMap_[To]][AgeGroupMap_[From]];
		}

		double HomeAgeContacts_ToFrom(const std::string& To, const std::string& From) {
			return waifw_home_[AgeGroupMap_[To]][AgeGroupMap_[From]];
		}

		double SDistAgeContacts_ToFrom(const std::string& To, const std::string& From) {
			return waifw_sdist_[AgeGroupMap_[To]][AgeGroupMap_[From]];
		}

		double ProbabilityOfHospitalisation_Pre_20 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Pre_20"]][0];
		double ProbabilityOfHospitalisation_20_29 =
			cfr_byage_[AgeGroupMap_["AgeGroup_20_29"]][0];
		double ProbabilityOfHospitalisation_30_39 =
			cfr_byage_[AgeGroupMap_["AgeGroup_30_39"]][0];
		double ProbabilityOfHospitalisation_40_49 =
			cfr_byage_[AgeGroupMap_["AgeGroup_40_49"]][0];
		double ProbabilityOfHospitalisation_50_59 =
			cfr_byage_[AgeGroupMap_["AgeGroup_50_59"]][0];
		double ProbabilityOfHospitalisation_60_69 =
			cfr_byage_[AgeGroupMap_["AgeGroup_60_69"]][0];
		double ProbabilityOfHospitalisation_Post_70 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Post_70"]][0];
		double ProbabilityOfHospitalisation_HealthCareWorkers =
			cfr_byage_[AgeGroupMap_["AgeGroup_HealthCareWorkers"]][0];

		double CaseFatalityRatio_Pre_20 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Pre_20"]][1];
		double CaseFatalityRatio_20_29 =
			cfr_byage_[AgeGroupMap_["AgeGroup_20_29"]][1];
		double CaseFatalityRatio_30_39 =
			cfr_byage_[AgeGroupMap_["AgeGroup_30_39"]][1];
		double CaseFatalityRatio_40_49 =
			cfr_byage_[AgeGroupMap_["AgeGroup_40_49"]][1];
		double CaseFatalityRatio_50_59 =
			cfr_byage_[AgeGroupMap_["AgeGroup_50_59"]][1];
		double CaseFatalityRatio_60_69 =
			cfr_byage_[AgeGroupMap_["AgeGroup_60_69"]][1];
		double CaseFatalityRatio_Post_70 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Post_70"]][1];
		double CaseFatalityRatio_HealthCareWorkers =
			cfr_byage_[AgeGroupMap_["AgeGroup_HealthCareWorkers"]][1];

		double ProbabilityOfDeath_Pre_20 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Pre_20"]][2];
		double ProbabilityOfDeath_20_29 =
			cfr_byage_[AgeGroupMap_["AgeGroup_20_29"]][2];
		double ProbabilityOfDeath_30_39 =
			cfr_byage_[AgeGroupMap_["AgeGroup_30_39"]][2];
		double ProbabilityOfDeath_40_49 =
			cfr_byage_[AgeGroupMap_["AgeGroup_40_49"]][2];
		double ProbabilityOfDeath_50_59 =
			cfr_byage_[AgeGroupMap_["AgeGroup_50_59"]][2];
		double ProbabilityOfDeath_60_69 =
			cfr_byage_[AgeGroupMap_["AgeGroup_60_69"]][2];
		double ProbabilityOfDeath_Post_70 =
			cfr_byage_[AgeGroupMap_["AgeGroup_Post_70"]][2];
		double ProbabilityOfDeath_HealthCareWorkers =
			cfr_byage_[AgeGroupMap_["AgeGroup_HealthCareWorkers"]][2];


		double SusceptibilityProbability_Pre_20 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_Pre_20"]];
		double SusceptibilityProbability_20_29 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_20_29"]];
		double SusceptibilityProbability_30_39 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_30_39"]];
		double SusceptibilityProbability_40_49 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_40_49"]];
		double SusceptibilityProbability_50_59 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_50_59"]];
		double SusceptibilityProbability_60_69 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_60_69"]];
		double SusceptibilityProbability_Post_70 =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_Post_70"]];
		double SusceptibilityProbability_HealthCareWorkers =
			pf_byage_[HealthBoardMap_[health_board_label_]][AgeGroupMap_["AgeGroup_HealthCareWorkers"]];

	private:
		std::string health_board_label_;
		std::vector<std::vector<int>> cases_;
		std::vector<std::vector<int>> deaths_;
		std::vector<std::vector<double>> age_pop_;
		std::vector<std::vector<double>> waifw_norm_;
		std::vector<std::vector<double>> waifw_home_;
		std::vector<std::vector<double>> waifw_sdist_;
		std::vector<std::vector<double>> cfr_byage_;
		std::vector<std::vector<double>> pf_byage_;
		std::map<std::string, int> HealthBoardMap_{{"Health Board 0", 0},
												   {"Health Board 1", 1},
												   {"Health Board 2", 2},
												   {"Health Board 3", 3},
												   {"Health Board 4", 4},
												   {"Health Board 5", 5},
												   {"Health Board 6", 6},
												   {"Health Board 7", 7},
												   {"Health Board 8", 8},
												   {"Health Board 9", 9},
												   {"Health Board 10", 10},
												   {"Health Board 11", 11},
												   {"Health Board 12", 12},
												   {"Health Board 13", 13},
												   {"Health Board 14", 14}};
		std::map<std::string, int> AgeGroupMap_{{"AgeGroup_Pre_20", 0},
												{"AgeGroup_20_29", 1},
												{"AgeGroup_30_39", 2},
												{"AgeGroup_40_49", 3},
												{"AgeGroup_50_59", 4},
												{"AgeGroup_60_69", 5},
												{"AgeGroup_Post_70", 6},
												{"AgeGroup_HealthCareWorkers", 7}};
};

} // namespace EERAModel
