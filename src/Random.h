#pragma once

#include <memory>
#include <gsl/gsl_rng.h>

namespace EERAModel {
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
};

/**
 * @class GSLRNG 
 * @brief Wrapper for the GNU Scientific Library random number generator
 */
class GSLRNG : public RNGInterface
{
public:
    using Sptr = std::shared_ptr<GSLRNG>;

    /**
     * @brief Contructor
     * 
     * Sets up a handle for the underlying GSL random number generator. Seeds the randomiser
     * with the supplied seed
     * 
     * @param seed Seed value
     */
    GSLRNG(unsigned long int seed);

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

private:
    /**
     * @private
     * @brief Handle for the GSL random number generator
     */
    gsl_rng* r_;
};

} // namespace Random
} // namespace EERAModel