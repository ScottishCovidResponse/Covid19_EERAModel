#include "InferenceParameters.h"
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

#if 0
std::vector<double> InferenceParameterGenerator::GenerateWeighted(
    std::vector<particle> particleList,
    int pick_val, 
    const std::vector<double>& vlimitKernel,
    const std::vector<double>& vect_min,
    const std::vector<double>& vect_Max) 
{
	std::vector<particle> particleList, int pick_val, const std::vector<double>& vlimitKernel, 
	const std::vector<double>& vect_min, const std::vector<double>& vect_Max) {

	std::vector<double> selected_param(nPar, 0);

    //define the elements of the chosen particle
    particle perturbList = particleList[pick_val];

	//perturb all elements of the chosen particle with a kernel following a uniform distribution +/- kernelFactor

    for (int xx = 0; xx < nPar; ++xx) {
        double tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        while(std::isnan(tmpval) || tmpval<=vect_min[xx] || tmpval>=vect_Max[xx]){
        	tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        }

    	//return a vector with all selection measures
    	selected_param[xx] = tmpval;
	}

	return selected_param;
}
#endif
} // namespace Inference
} // namespace EERAModel
