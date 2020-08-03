#include "Random.h"
#include <gsl/gsl_randist.h>

namespace EERAModel {
namespace Random {

RNG::RNG(unsigned long int seed)
 : gen_(seed)
{
    r_ = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r_, seed);
}

double RNG::Flat(double a, double b)
{
    return gsl_ran_flat(r_, a, b);
}

void RNG::Multinomial(size_t K, unsigned int N, const double p[], unsigned int n[])
{
    gsl_ran_multinomial(r_, K, N, p, n);
}

double RNG::Poisson(double mu)
{
    return gsl_ran_poisson(r_, mu);
}

double RNG::Gamma(double a, double b)
{
    return gsl_ran_gamma(r_, a, b);
}

double RNG::Beta(double a, double b)
{
    return gsl_ran_beta(r_, a, b);
}

unsigned int RNG::Binomial(double p, unsigned int n)
{
    return gsl_ran_binomial(r_, p, n);
}

std::mt19937& RNG::MT19937()
{
    return gen_;
}

} // namespace Random
} // namespace EERAModel