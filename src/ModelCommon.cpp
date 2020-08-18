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

std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const ObservationsForModels& obs) {
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

int GetPopulationOfRegion(const ObservationsForModels& obs, int region_id)
{
	return obs.cases[region_id][0];
}

int ComputeNumberOfHCWInRegion(int regionalPopulation, int totalHCW, const ObservationsForModels& obs)
{
    int scotlandPopulation = 0;

	// TODO: In the csv files, obs.cases[0] contains column headings so probably shouldn't
	// be used in this sum. The value of obs.cases[0][0] will be -1, so unlikely to show
	// much?

	// for (unsigned int region = 0; region < obs.cases.size() - 1; ++region) {
	for (unsigned int region = 1; region < obs.cases.size() - 1; ++region) {
		scotlandPopulation += obs.cases[region][0];
	}
	double regionalProportion = static_cast<double>(regionalPopulation) / scotlandPopulation;
	
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

std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec)
{
	std::vector<std::vector<int>> _temp;

	for(auto cmps : cmps_vec)
	{
		_temp.push_back({cmps.S, cmps.E, cmps.E_t, cmps.I_p,
						cmps.I_t, cmps.I1, cmps.I2, cmps.I3,
						cmps.I4, cmps.I_s1, cmps.I_s2, cmps.I_s3,
						cmps.I_s4, cmps.H, cmps.R, cmps.D});
	}

	return _temp;
}

} // namespace Model
} // namespace EERAModel