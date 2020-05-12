#pragma once

#include "ModelTypes.h"

#include <vector>
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace FittingProcess {

void parameter_select_first_step(std::vector<double> &selected_param, std::vector<double> flag1,
	std::vector<double> flag2, gsl_rng * r, int nPar);

void parameter_select_nsteps(std::vector<double> &selected_param, int nPar, gsl_rng * r,
	std::vector<::EERAModel::particle> particleList, int pick_val, 
	const std::vector<double>& vlimitKernel, const std::vector<double>& vect_min,
	const std::vector<double>& vect_Max);

void weight_calc(int smc,int pastNpart, std::vector<::EERAModel::particle> pastPart,
	EERAModel::particle &currentPart, const std::vector<double>& vlimitKernel, int nPar);

} // namespace FittingProcess
} // namespace EERAModel