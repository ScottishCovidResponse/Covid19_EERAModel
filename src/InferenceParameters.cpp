#include "InferenceParameters.h"
#include "ModelCommon.h"
#include <cmath>

namespace EERAModel {
namespace Inference {

InferenceParameterGenerator::InferenceParameterGenerator(
    Random::RNGInterface::Sptr rng, const std::vector<double> &flag1,
    const std::vector<double> &flag2)
    : rng_(rng), flag1_(flag1), flag2_(flag2) {}

std::vector<double> InferenceParameterGenerator::GenerateInitial() {
  std::vector<double> parameterSet(Model::ModelParameters::NPARAMS);

  for (unsigned int i = 0; i < Model::ModelParameters::NPARAMS; ++i) {
    if (Model::ModelParameters::CHCW == i) {
      parameterSet[i] = (double)rng_->Poisson(flag1_[i]);
    } else if (Model::ModelParameters::LAMBDA == i) {
      parameterSet[i] = rng_->Flat(flag1_[i], flag2_[i]);
    } else if (Model::ModelParameters::RRD == i) {
      parameterSet[i] = rng_->Gamma(flag1_[i], flag2_[i]);
    } else {
      parameterSet[i] = rng_->Beta(flag1_[i], flag2_[i]);
    }
  }

  return parameterSet;
}

std::vector<double> InferenceParameterGenerator::GenerateWeighted(
    const std::vector<double> &existingSet,
    const std::vector<KernelWindow> &kernelWindows) {
  std::vector<double> parameterSet(Model::ModelParameters::NPARAMS);

  for (unsigned int i = 0; i < Model::ModelParameters::NPARAMS; ++i) {
    parameterSet[i] = PerturbParameter(existingSet[i], kernelWindows[i]);
  }

  return parameterSet;
}

double
InferenceParameterGenerator::PerturbParameter(double oldParam,
                                              const KernelWindow &window) {
  double parameter;

  do {
    parameter = oldParam + window.kernel * rng_->Flat(-1, 1);
  } while (parameter <= window.min || parameter >= window.max);

  return parameter;
}

} // namespace Inference
} // namespace EERAModel
