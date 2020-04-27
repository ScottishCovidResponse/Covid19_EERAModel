#pragma once

#include "ModelTypes.h"
#include <vector>
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace Model {

void my_model(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				std::vector<std::vector<double>> cfr_byage, std::vector<double> pf_byage,
				std::vector<std::vector<double>> waifw_norm,
				std::vector<std::vector<double>> waifw_sdist, std::vector<std::vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, std::vector<int> agenums, 
				double tau, gsl_rng * r, std::vector<int> &sim_status,std::vector<std::vector<int>> &ends,
				std::vector<int> &death_status,std::vector<int> &deathH_status);

} // namespace Model
} // namespace EERAModel