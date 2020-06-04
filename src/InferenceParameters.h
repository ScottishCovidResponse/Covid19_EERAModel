#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include <gsl/gsl_rng.h>
#include "ModelTypes.h"
#include "Random.h"

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
 * of disease spread in the world, then the corresponding inference parameter set is retained as 
 * part of the output from the inference framework.
 * 
 * Inference parameter sets can be generated in one of two ways:
 *   - An initial set of parameters can be generated based on input observations/priors
 *   - Once the inference framework has built up one or more sets of accepted inference parameters,
 *     subsequent sets can be generated from a selected parameter from that set, by perturbing the
 *     selected parameter set.
 */
class InferenceParameterGenerator
{
public:
    using Sptr = std::shared_ptr<InferenceParameterGenerator>;
    
    /**
     * @brief Constructor
     * 
     * @param nParameters Number of inference parameters to generate in each set
     * @param rng Random number generator
     * @param flag1 First data set of priors 
     * @param flag2 Second data set of priors
     */
    InferenceParameterGenerator(
        std::size_t nParameters,
        Random::RNGInterface::Sptr rng, 
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
     * @brief Generate a set of parameters by perturbing an existing set
     * 
     * @param existingSet Existing parameter set to perturb
     * @param Kernel for the perturbation (scales the range of the perturbation)
     * @param vect_Min Minimum allowed values for the perturbed parameters
     * @param vect_Max Maximum allowed values for the perturbed parameters
     * 
     * @return Parameter set
     */ 
    std::vector<double> GenerateWeighted(
        const std::vector<double>& existingSet,
	    const std::vector<double>& vlimitKernel,
	    const std::vector<double>& vect_Max,
        const std::vector<double>& vect_Min
    );

    /**
     * @brief Generate a single parameter value by perturbing an existing value
     * 
     * Generates a single parameter value, as a perturbation of an existing parameter value. The 
     * returned value is guaranteed to lie between @p min and @p max. The perturbation is performed
     * using a uniform value in the range +/- @p kernel.
     * 
     * @param oldParam Existing parameter value
     * @param kernel Scale factor for generation
     * @param max Upper limit for generated value
     * @param min Lower limit for generated value
     * 
     * @return Generated value
     */
    inline double PerturbParameter(double oldParam, double kernel, double max, double min);

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
    Random::RNGInterface::Sptr rng_;

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