#pragma once

#include "ModelTypes.h"
#include "Random.h"
#include <cstddef>
#include <vector>

namespace EERAModel {
/**
 * @brief Namespace handling inference paramaters
 *
 * This namespace contains functions to manage and generate
 * inference parameter sets which setup the model behaviour.
 */
namespace Inference {

/**
 * @class InferenceParameterGenerator
 * @brief Generate inference parameter sets on demand
 *
 * Inference parameter sets are collection of parameters used as inputs to given
 * run of the epidemiological model within the inference framework. The
 * parameter sets control the behaviour of the epidemiological model. Repeated
 * runs of the model with different sets of inference parameters result in
 * different model outputs. If the outputs closely match the observations of
 * disease spread in the world, then the corresponding inference parameter set
 * is retained as part of the output from the inference framework.
 *
 * Inference parameter sets can be generated in one of two ways:
 *   - An initial set of parameters can be generated based on input
 * observations/priors
 *   - Once the inference framework has built up one or more sets of accepted
 * inference parameters, subsequent sets can be generated from a selected
 * parameter from that set, by perturbing the selected parameter set.
 */
class InferenceParameterGenerator {
public:
  using Sptr = std::shared_ptr<InferenceParameterGenerator>;

  /**
   * @brief Constructor
   *
   * @param rng Random number generator
   * @param flag1 First data set of priors
   * @param flag2 Second data set of priors
   */
  InferenceParameterGenerator(Random::RNGInterface::Sptr rng,
                              const std::vector<double> &flag1,
                              const std::vector<double> &flag2);

  /**
   * @brief Generate an initial set of parameters
   *
   * The initial parameter set is generated for use on the first iteration of
   * the inference framework. The parameter set is generated based on
   * perturbations from an initial set of observational priors
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
   * @param kernelWindows Kernel windows for each parameter
   *
   * @return Parameter set
   */
  std::vector<double>
  GenerateWeighted(const std::vector<double> &existingSet,
                   const std::vector<KernelWindow> &kernelWindows);

  /**
   * @brief Generate a single parameter value by perturbing an existing value
   *
   * Generates a single parameter value, as a perturbation of an existing
   * parameter value. The returned value is guaranteed to lie between @p
   * window.min and @p window.max. The perturbation is performed using a uniform
   * value in the range +/- @p window.kernel.
   *
   * @param oldParam Existing parameter value
   * @param window Kernel window for the parameter
   *
   * @return Generated value
   */
  inline double PerturbParameter(double oldParam, const KernelWindow &window);

private:
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