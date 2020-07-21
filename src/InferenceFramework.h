#pragma once

#include "ModelCommon.h"
#include "ModelTypes.h"
#include "Random.h"
#include "Utilities.h"
#include "InferenceParameters.h"
#include <ctime>

namespace EERAModel {
namespace Inference {

/**
 * @class InferenceParticleGenerator
 * @brief Generate inference particles on demand
 */
class InferenceParticleGenerator
{
public:
    using Sptr = std::shared_ptr<InferenceParticleGenerator>;
    
    /**
     * @brief Constructor
     * 
     * @param nInferenceParams Number of inference parameters in each particle
     * @param kernelFactor Configuration factor for the parameter kernels
     * @param rng Seeded random number generator
     */
    InferenceParticleGenerator(unsigned int nInferenceParams, double kernelFactor,
        Random::RNGInterface::Sptr rng);

    /**
     * @brief Update the generator
     * 
     * On each loop of the ABC-smc framework, a new collection of particles is accepted by the 
     * framework and recorded. This function should be called with the new set of particles, updating
     * its internal weight distribution and kernel windows for the inference parameters
     * 
     * @param particles List of particles accepted on the last ABC-smc loop
     */
    void Update(std::vector<particle> particles);

    /**
     * @brief Generate a new inference particle
     * 
     * @param smc Inference loop number
     * @param sim Simulation loop number
     * @param inferenceParameterGenerator Parameter generator for th inference framework
     * @param previousParticles Previously accepted inference particles
     */
    particle GenerateNew(int smc, int sim, InferenceParameterGenerator::Sptr parameterGenerator,
      const std::vector<particle>& previousParticles);

    /**
     * @brief Get the current kernel windows
     * 
     * @return Const reference to the current set of kerenl windows
     */
    const std::vector<KernelWindow> KernelWindows() const { return kernelWindows_; }

private:
    /**
     * @private
     * @brief NUmber of inference parameters to generate in each particle
     */
    unsigned int nInferenceParams_;

    double kernelFactor_;

    Random::RNGInterface::Sptr rng_;
    
    /**
     * @private
     * @brief Kernel windows for accepted particles
     */
    std::vector<KernelWindow> kernelWindows_;
    
    /**
     * @private
     * @brief Distribution of accepted particles
     */
    std::discrete_distribution<int> weightDistribution_;
};

/**
 * @class InferenceFramework
 * @brief Parameter inference framework
 */
class InferenceFramework
{
public:    
    /**
     * @brief Framework constructor
     *
     * @param model Model interface
     * @param inferenceConfig Inference config
     * @param rng Seeded random number generator
     * @param outDir Outputs directory path
     * @param log Logger
     */
    InferenceFramework(Model::ModelInterface::Sptr model,
        const InferenceConfig& inferenceConfig,
        Random::RNGInterface::Sptr rng,
        const std::string& outDir,
        Utilities::logging_stream::Sptr log);
    
    /**
     * @brief Calculate the time offset for the dataset
     * 
     * Determines the offset for the start of the dataset based on
     * the given parameters
     * 
     * @param inferenceConfig Inference config
     * 
     * @return integer offset value
     */
    static int GetTimeOffSet(const InferenceConfig& inferenceConfig);

    /**
     * @brief Run the model within the inference framework
     */
     void Run();

private:
    /**
     * @private
     * @brief Run the model inside the inference framework
     * 
     * Runs the model and computes the measures associated with the difference between the model 
     * outputs and the observations
     */
    void ModelSelect(EERAModel::particle& outvec, int n_sim_steps, 
        seed seedlist, int day_shut, const std::vector<int>& obsHosp, const std::vector<int>& obsDeaths);

    /**
     * @private
     * @brief Check if an inference particle passes the tolerance limits
     * 
     * A particle passes the tolerance limits if:
     *   - The sum of squared errors (SSE) between the observed and simulated hospitalisations
     * is less than the configured tolerance limit for the given @p smc step
     *   - The sum of squared errors between the observed and simulated deaths is less than the
     * configured tolerance limit for the given @p smc step
     * 
     * @param p Particle under consideration
     * @param smc NUmber of the inference step ongoing
     * 
     * @return True if the particle passes the tolerance limits; otherwise false
     */
    bool ParticlePassesTolerances(const particle& p, int smc);

    /**
     * @private
     * @brief Compute the weight of a particle
     * 
     * The weight of a given particle is the weighted sum of the weights of all previous particles which
     * are "close" to the given particle. 
     * 
     * @note If @p smc is zero, the a weight of 1.0 is returned. 
     * 
     * @param smc ABC-smc loop number
     * @param previousParticles Previously accepted particles
     * @param currentParticle Particle for which we wna tot compute the weight
     * @param kernelWindows Kernel windows for the inference parameters in the particles
     */
    double ComputeParticleWeight(int smc, const std::vector<::EERAModel::particle>& previousParticles,
	    const particle &currentParticle, const std::vector<KernelWindow> kernelWindows);

    /**
     * @private
     * @brief Determine if two particles are close to each other
     * 
     * Closeness is defined as follows. Two particles are "close" if, for each parameter in their 
     * respective parameter sets, the absolute difference between corresponding pairs of particles
     * is less than the corresponding element of the kernel. In other words, if a given particle
     * could have been generated by perturbation of another particle, the two particles are close.
     * 
     * @param first First particle of the pair
     * @param second Second particle of the pair
     * @param kernelWindows Kernel windows for the inference parameters in the particles
     * 
     * @return True if the two particles are close, otherwise false
     */
    bool ParticlesAreClose(const particle& first, const particle& second,
        const std::vector<KernelWindow> kernelWindows);

    /**
     * @brief Update the parameter kernel windows
     * 
     * Compute a new set of kernel windows, based on the accepted particles from the last ABC-smc
     * loop. Update the framework's kernel windows with the result of the computation
     */
    void UpdateKernelWindows(int smc, std::vector<particle> particles);
    
    /**
     * @private
     * @brief Model interface
     */
    Model::ModelInterface::Sptr model_;

    /**
     * @private
     * @brief Inference configuration
     */
    InferenceConfig inferenceConfig_;

    /**
     * @private
     * @brief Random number generator
     */
    Random::RNGInterface::Sptr rng_;

    /**
     * @private
     * @brief Outputs directory path
     */
    std::string outDir_;

    /**
     * @private
     * @brief Logger
     */
    Utilities::logging_stream::Sptr log_;

    /**
     * @private
     * Inference parameter generator
     */
    InferenceParameterGenerator::Sptr inferenceParameterGenerator_;

    /**
     * @private
     * @brief Inference particle generator
     */
    InferenceParticleGenerator::Sptr inferenceParticleGenerator_;
};

/**
 * @brief Compute the Kernel Windows for each inference parameter
 * 
 * Calculates the kernel window (kernel, max limit and min limit) for each inference parameter,
 * based on a set of previously accepted inference particles.
 * 
 * @param nPar Number of parameters to compute kernel for
 * @param particleList List of previously accepted particles
 * @param kernelFactor Common kernel multiplicative factor
 * 
 * @return List of computed kernel windows for each inference parameter
 */
std::vector<KernelWindow> ComputeKernelWindow(int nPar, const std::vector<particle>& particleList,
	double kernelFactor);

/**
 * @brief Compute weight distribution
 * 
 * Create the discrete distribution of the weights for the "importance sampling" process
 * 
 * @param particleList List of previously accepted particles
 * 
 * @return Weight distribution
 */
std::discrete_distribution<int> ComputeWeightDistribution(
	const std::vector<EERAModel::particle>& particleList);

} // namespace Inference
} // namespace EERAModel
