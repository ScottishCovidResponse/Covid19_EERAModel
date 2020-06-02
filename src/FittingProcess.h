#pragma once

#include "ModelTypes.h"

#include <vector>
#include <gsl/gsl_rng.h>

namespace EERAModel {
namespace FittingProcess {

void weight_calc(int smc,int pastNpart, std::vector<::EERAModel::particle> pastPart,
	EERAModel::particle &currentPart, const std::vector<double>& vlimitKernel, int nPar);

} // namespace FittingProcess
} // namespace EERAModel