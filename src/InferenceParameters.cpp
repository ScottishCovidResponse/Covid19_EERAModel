#include "InferenceParameters.h"
#include <cmath>
#include <gsl/gsl_randist.h>

namespace EERAModel {
namespace Inference {

InferenceParameterGenerator::InferenceParameterGenerator(
    std::size_t nParameters,
    Random::RNGInterface::Sptr rng,
    const std::vector<double>& flag1,
    const std::vector<double>& flag2) 
     : nParameters_(nParameters), rng_(rng), flag1_(flag1), flag2_(flag2) {}

std::vector<double> InferenceParameterGenerator::GenerateInitial()
{
    std::vector<double> parameterSet(nParameters_, 0.0);

    for (std::size_t i = 0; i < nParameters_; ++i) {
        double tmpval = 0.0;

        if (i == 2) 
        {
            tmpval = (double) rng_->Poisson(flag1_[i]);
        }
        else if (i == (nParameters_ - 1)) 
        {
            tmpval = rng_->Flat(flag1_[i], flag2_[i]);
        } 
        else if (i == 6)
        {
            tmpval = rng_->Gamma(flag1_[i], flag2_[i]);
        }
        else
        {
            tmpval = rng_->Beta(flag1_[i], flag2_[i]);
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
        parameter = oldParam + kernel * rng_->Flat(-1, 1);
    } 
    while (std::isnan(parameter) || parameter <= min || parameter >= max);
    
    return parameter;
}

} // namespace Inference
} // namespace EERAModel
