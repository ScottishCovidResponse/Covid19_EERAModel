#include "ModelCommon.h"
#include <cmath>

namespace EERAModel {
namespace Model {

int Flow(Random::RNGInterface::Sptr rng, int pops_from_val, int pops_new_from_val, double rate)
{
	int outs = rng->Poisson(rate * static_cast<double>(pops_from_val)); //symptomatic become hospitalized
	outs = std::min(pops_new_from_val, outs);
	return outs;
}

std::vector<params> BuildFixedParameters(unsigned int size, params parameters)
{
    return std::vector<params>(size, parameters);
}

std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const InputObservations& obs) {
	std::vector<int> agenums;
	
	// define age structure of the shb of interest. the -1 is to account for difference in number of
	// rows between two datasets (age_pop does not have a row with column name)
	const auto& agedist = obs.age_pop[shb_id - 1];
	
	for (const auto& var : agedist) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(var * (Npop - N_hcw))); 

	}		
	agenums.push_back(N_hcw);

	return agenums;
}

int GetPopulationOfRegion(const InputObservations& obs, int region_id)
{
	return obs.cases[region_id][0];
}

int ComputeNumberOfHCWInRegion(int regionalPopulation, int totalHCW, const InputObservations& observations)
{
    int scotlandPopulation = 0;
	for (unsigned int region = 0; region < observations.cases.size() - 1; ++region) {
		scotlandPopulation += observations.cases[region][0];
	}
	double regionalProportion = static_cast<double>(regionalPopulation) / scotlandPopulation;
	if(observations.cases.size()-1 != 15){
		regionalProportion = 1.0;
	}
	
    return static_cast<int>(round(totalHCW * regionalProportion)); 
}

int accumulate_compartments(const Compartments& comp)
{
	int _total = 0;
	_total += comp.S + comp.E + comp.E_t + comp.I_p;
	_total += comp.I_t + comp.I1 + comp.I2 + comp.I3;
	_total += comp.I4 + comp.I_s1 + comp.I_s2 + comp.I_s3;
	_total += comp.I_s4 + comp.H + comp.R + comp.D;

	return _total;
}

} // namespace Model
} // namespace EERAModel