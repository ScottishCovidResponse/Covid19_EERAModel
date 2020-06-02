#include "Random.h"
#include <gsl/gsl_randist.h>

namespace EERAModel {
namespace Random {

GSLRNG::GSLRNG(unsigned long int seed)
{
    r_ = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r_, seed);
}

double GSLRNG::Flat(double a, double b)
{
    return gsl_ran_flat(r_, a, b);
}

void GSLRNG::Multinomial(size_t K, unsigned int N, const double p[], unsigned int n[])
{
    gsl_ran_multinomial(r_, K, N, p, n);
}

double GSLRNG::Poisson(double mu)
{
    return gsl_ran_poisson(r_, mu);
}

double GSLRNG::Gamma(double a, double b)
{
    return gsl_ran_gamma(r_, a, b);
}

double GSLRNG::Beta(double a, double b)
{
    return gsl_ran_beta(r_, a, b);
}

} // namespace Random
} // namespace EERAModel