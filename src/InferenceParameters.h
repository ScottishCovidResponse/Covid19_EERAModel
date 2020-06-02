#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include <gsl/gsl_rng.h>
#include "ModelTypes.h"

namespace EERAModel {
namespace Inference {

/**
 * @class InferenceParameterGenerator
 * @brief Generate inference parameter sets on demand
 * 
 * Inference parameter sets are collection of parameters used as inputs to given run of the 
 * epidemiological model within the inference framework. The parameter sets control the behaviour
 * of the epidemiological model. Repeated runs of the model with different sets of inference 
 * parameters result in different model outputs. If the outputs closely match the observations
 * of disese spread in the world, then the corresponding inference parameter set is retained as 
 * part of the output from the inference framework.
 * 
 * Inference parameter sets can be generated in one of two ways:
 *   - An initial set of parameters can be generated based on input observations/priors
 *   - Once the inference framework has built up one or more sets of accepted inference parameters,
 *     subsequent sets can be generated from these accepted sets, using a weighted generation process,
 *     where the weight attached to a previously accepted set determines its contribution to the
 *     newly generated set.
 */
class InferenceParameterGenerator
{
public:
    using Sptr = std::shared_ptr<InferenceParameterGenerator>;
    
    /**
     * @brief Constructor
     * 
     * @param nParameters Number of inference parameters to generate in each set
     * @param flag1 First data set of priors 
     * @param flag2 Second data set of priors
     */
    InferenceParameterGenerator(
        std::size_t nParameters,
        gsl_rng* r, 
        const std::vector<double>& flag1,
        const std::vector<double>& flag2
    );

    /**
     * @brief Generate an initial set of parameters
     * 
     * The initial parameter set is generated for use on the first iteration of the inference
     * framework. The parameter set is generated based on perturbations from an initial set of
     * observational priors
     * 
     * Randomly picks parameters' value from priors with gamma shape
     * 
     * @return Parameter set
     */
    std::vector<double> GenerateInitial();

    /**
     * @brief Generate a weighted set of parameters
     * 
     * @return Parameter set
     */ 
    std::vector<double> GenerateWeighted(
        std::vector<particle> acceptedParticles,
        int pick_val, 
	    const std::vector<double>& vlimitKernel,
        const std::vector<double>& vect_min,
	    const std::vector<double>& vect_Max
    ) {}

private:
    /**
     * @private
     * @brief Number of parameters to generate in each set
     */
    std::size_t nParameters_;

    /**
     * @private
     * @brief Local random number generator
     */
    gsl_rng* r_;

    /**
     * @private
     * @brief First set of priors
     */
    const std::vector<double> flag1_;

    /**
     * @private
     * @brief Second set of priors
     */
    const std::vector<double> flag2_;
};

} // namespace Inference
} // namespace EERAModel