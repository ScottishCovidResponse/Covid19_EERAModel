#include "ModelCommon.h"
#include "ModelTypes.h"
#include "Random.h"
#include <memory>
#include <vector>

namespace EERAModel {
namespace Model {

class OriginalModel : public ModelInterface
{
public:
    using Sptr = std::shared_ptr<OriginalModel>;

    /**
     * @brief Constructor
     * 
     * @param rng Random number generator to be used by the model
     */
    OriginalModel(Random::RNGInterface::Sptr rng) : rng_(rng) {}

    /**
     * @brief Run the model with the given parameters and configurations
     * 
     * @param parameter_set: set of parameters that are being infered (i.e. particles)
     * @param fixed_parameters: set of fixed (known) parameters
     * @param per_age_data: age-structured known parameters (such as case fatality ratio (CFR) and probability of severe clinical outcomes )
     * @param seedlist: seeding method to initialise infection ("random": randomly allocate n infectious individuaals at the start, or "background": background transmission over a defined pre-lockdown period )
     * @param day_shut: day of the lock down
     * @param agenums: number of individuals in each age group in the study area
     * @param n_sim_steps: Number of steps to simulate
     * 
     * @return Status of model after run
     */
    Status Run(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
                    AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
                    int n_sim_steps) override;

private:
    /**
     * @brief Construct the population seed
     * 
     * The population seed is constructed according to the following criteria:
     *    - The final age category (HCW) are omitted from the seeding always
     *    - The first age category (< 20yo) are omitted 
     * 
     * @param age_nums Vector containing the number of people in each age category
     * 
     * @return Seed population
     */
    std::vector<double> BuildPopulationSeed(const std::vector<int>& age_nums);

    /**
     * @brief Construct the population array
     * 
     * Sets up an array for the population at each timestep in each age and disease category	
     * also set up the age distribution of old ages as target for disease introduction.
     * 
     * @param age_nums Vector containing the number of people in each age group
     * @param seedlist Seed object
     * 
     * @return Vector of vectors containing compartment populations
     */ 
    std::vector<Compartments> BuildPopulationArray(const std::vector<int>& age_nums,
        const seed& seedlist);

    /**
     * @brief Introduced diseased to the population
     * 
     * Compute the total number of susceptible and the number of susceptible per age class
     * 
     * @param poparray Population array to be manipulated
     * @param seedarray Population seed array to be manipulated
     * @param bkg_lambda Lambda for generating number of diseased individuals
     */
    void GenerateDiseasedPopulation(std::vector<Compartments>& poparray,
        std::vector<double>& seedarray, const double& bkg_lambda);

    /**
     * @brief Generate vector of lambda values for the age groups
     * 
     * Creates lambda values based on compartment occupancy for each age group
     * 
     * @param inf_hosp Number of hospitalised infected
     * @param parameter_set Set of model parameters
     * @param u_val
     * @param age_data Age group data set
     * @param pops Population array containing compartments for each age group
     * @param shut State of lockdown
     */
    std::vector<double> GenerateForcesOfInfection(int& inf_hosp, const std::vector<double>& parameter_set, double u_val, 
                const AgeGroupData& age_data, const std::vector<Compartments>& pops, bool shut);

    /**
     * @brief Generate an infection spread as per the structure of the original EERA model
     * 
     * Creates an infection spread state and counters number of people in different states
     * 
     * @param pop Particular population age group
     * @param n_hospitalised Number of people in hospital prior to new spread
     * @param fixed_parameters Fixed model parameters
     * @param parameter_set Variable model parameters
     * @param cfr_tab Case Fatality Ratio table
     * @param pf_val Frailty Probability
     * @param lambda Rate of spread
     */
    InfectionState GenerateInfectionSpread(Compartments& pop,
        const int& n_hospitalised, ::EERAModel::params fixed_parameters, 
        std::vector<double> parameter_set, std::vector<double> cfr_tab,
        double pf_val, double lambda);

    /**
     * @private
     * @brief Random number generator
     */
    Random::RNGInterface::Sptr rng_;
};

} // namespace Model
} // namespace EERAModel
