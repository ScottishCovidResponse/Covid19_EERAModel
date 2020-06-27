#pragma once

#include <memory>
#include <random>
#include <gsl/gsl_rng.h>

namespace EERAModel {
/**
 * @brief Namespace containing objects related to random number generation
 *
 * This namespace contains a wrapper interace to existing random number
 * generation libraries in order for them to be used within all
 * model functions.
 */
namespace Random {

/**
 * @class RNGInterface
 * @brief Abstract interface to a random number generator
 */
class RNGInterface {
public:
    using Sptr = std::shared_ptr<RNGInterface>;

    virtual ~RNGInterface() = default;

    /**
     * @brief Flat distribution
     * 
     * Return a random number from a flat distribution
     * 
     * @param a Lower limit of distribution interval
     * @param b Upper limit of distribution interval
     * 
     * @return Random number
     */
    virtual double Flat(double a, double b) = 0;

    /**
     * @brief Multinomial distribution
     * 
     * Return a random number from a multinomial distribution
     */
    virtual void Multinomial(size_t K, unsigned int N, const double p[], unsigned int n[]) = 0;

    /**
     * @brief Poisson distribution
     * 
     * Return a random number from a Poisson distribution
     * 
     * @param mu Mean of the distribution
     * 
     * @return Random number
     */
    virtual double Poisson(double mu) = 0;

    /**
     * @brief Gamma distribution
     *
     * Return a random number from a Gamma distribution
     * 
     * @param a Gamma a parameter
     * @param b Gamma b parameter
     * 
     * @return Random number
     */
    virtual double Gamma(double a, double b) = 0;

    /**
     * @brief Beta distribution
     * 
     * Return a random number from a Beta distribution
     * 
     * @param a Beta a parameter
     * @param b Beta b parameter
     * 
     * @return Random number
     */
    virtual double Beta(double a, double b) = 0;

   /**
    * @brief Binomial distribution
    * 
    * Return a random number from a Binomial distribution
    * 
    * @param p Binomial p parameter: probability
    * @param n Binomial n parameter: the population to sample
    * 
    * @return Random number
    */
   virtual unsigned int Binomial(double p, unsigned int n) = 0;

    /**
     * @brief MT19937 generator
     * 
     * Returns a reference to a std::mt19937 object
     * 
     * @return mt19937 reference
     */
    virtual std::mt19937& MT19937() = 0;
};

/**
 * @class RNG 
 * @brief Wrapper for the GNU Scientific Library random number generator
 */
class RNG : public RNGInterface
{
public:
    using Sptr = std::shared_ptr<RNG>;

    /**
     * @brief Contructor
     * 
     * Sets up a handle for the underlying GSL and STL random number generator. Seeds the randomiser
     * with the supplied seed
     * 
     * @param seed Seed value
     */
    RNG(unsigned long int seed);

    /** @brief Flat distribution override */
    virtual double Flat(double a, double b) override;

    /** @brief Multinomial distribution override */
    virtual void Multinomial(size_t K, unsigned int N, const double p[], unsigned int n[]) override;

    /** @brief Poisson distribution override */
    virtual double Poisson(double mu) override;

    /** @brief Gamma distribution override */
    virtual double Gamma(double a, double b) override;

    /** @brief Beta distribution override */
    virtual double Beta(double a, double b) override;

    /** @brief Binomial distribution override */
    virtual unsigned int Binomial(double p, unsigned int n) override;
	
    /** @brief MT19937 override */
    virtual std::mt19937& MT19937() override;
    
private:
    /**
     * @private
     * @brief Handle for the GSL random number generator
     */
    gsl_rng* r_;

    /**
     * @private
     * @brief Mersenne twister generator
     */
    std::mt19937 gen_;
};

} // namespace Random
} // namespace EERAModel