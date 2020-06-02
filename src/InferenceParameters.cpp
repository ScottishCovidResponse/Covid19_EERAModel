#include "InferenceParameters.h"
#include <cmath>
#include <gsl/gsl_randist.h>

namespace EERAModel {
namespace Inference {

InferenceParameterGenerator::InferenceParameterGenerator(
    std::size_t nParameters,
    gsl_rng* r,
    const std::vector<double>& flag1,
    const std::vector<double>& flag2) 
     : nParameters_(nParameters), r_(r), flag1_(flag1), flag2_(flag2) {}

std::vector<double> InferenceParameterGenerator::GenerateInitial()
{
    std::vector<double> parameterSet(nParameters_, 0.0);

    for (std::size_t i = 0; i < nParameters_; ++i) {
        double tmpval = 0.0;

        if (i == 2) 
        {
            tmpval = (double) gsl_ran_poisson(r_, flag1_[i]);
        }
        else if (i == (nParameters_ - 1)) 
        {
            tmpval = gsl_ran_flat(r_, flag1_[i], flag2_[i]);
        } 
        else if (i == 6)
        {
            tmpval = gsl_ran_gamma(r_, flag1_[i], flag2_[i]);
        }
        else
        {
            tmpval = gsl_ran_beta(r_, flag1_[i], flag2_[i]);
        }
        
        parameterSet[i] = tmpval;
	}

	return parameterSet;
}

std::vector<double> InferenceParameterGenerator::GenerateWeighted(
    const std::vector<double>& existingSet,
    const std::vector<double>& vlimitKernel,
    const std::vector<double>& vect_Max,
    const std::vector<double>& vect_Min) 
{
    std::vector<double> parameterSet(nParameters_, 0.0);
    
    for (size_t i = 0; i < nParameters_; ++i) 
    {
        parameterSet[i] = PerturbParameter(
            existingSet[i], vlimitKernel[i], vect_Max[i], vect_Min[i]
        );
    }

    return parameterSet;
}

double InferenceParameterGenerator::PerturbParameter(double oldParam, double kernel, double max, double min)
{
    double parameter;
    
    do
    {
        parameter = oldParam + kernel * gsl_ran_flat(r_, -1, 1);
    } 
    while (std::isnan(parameter) || parameter <= min || parameter >= max);
    
    return parameter;
}

} // namespace Inference
} // namespace EERAModel
