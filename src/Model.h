#pragma once

#include "ModelTypes.h"
#include <vector>
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace Model {

void model_select(int smc, ::EERAModel::particle &outvec, std::vector<params> fixed_parameters,
	std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage, 
	std::vector<std::vector<double>> waifw_norm, std::vector<std::vector<double>> waifw_sdist,
	std::vector<std::vector<double>> waifw_home, std::vector <int> agenums, double tau,
	int duration, seed seedlist, int day_shut, int Npop, gsl_rng * r, const std::vector<int>& obsHosp,
	const std::vector<int>& obsDeaths);

void select_obs(int& Npop, int& t_index, int& duration, int& day_intro, int& day_shut, 
	std::vector<int>& obsHosp_tmp, std::vector<int>& obsDeaths_tmp, 
	std::vector<std::vector<int> > data_tmp, std::vector<std::vector<int> > death_tmp, int herd_id, 
	int time_back);

} // namespace Model
} // namespace EERAModel