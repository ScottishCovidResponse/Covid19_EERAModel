#include "Observations.h"
#include <algorithm>
#include <numeric>

namespace EERAModel {
namespace Observations {

SimTime GetSimDuration(int time_back, const std::vector<int>& regional_cases)
{
	SimTime sim_time_pars;

	// Find the last non-zero value within the Regional Cases dataset
	auto max_time_itr = std::find_if(regional_cases.rbegin(), regional_cases.rend(), [](int n) {return n > 0;});
	const int maximum_time = std::distance(max_time_itr, regional_cases.rend());

	// Find the first positive increase in Regional Cases dataset after the start day
	auto first_case_itr = std::find_if(regional_cases.begin()+1, regional_cases.end(), [](int n) {return n > 0;});
	sim_time_pars.t_index = std::distance(regional_cases.begin(), first_case_itr);

	// Apply Offset
	sim_time_pars.day_intro = sim_time_pars.t_index - time_back;

	// Add 1 to account for Day 0
	sim_time_pars.day_intro += 1;

	// If offset results in negative start time calculate disease duration using
	// this start time, else duration is just the time of the last non-negative
	// value for number of cases
	sim_time_pars.duration = (sim_time_pars.day_intro < 0) ? maximum_time - sim_time_pars.day_intro : maximum_time;
	
	// Set start time to be zero if negative
	sim_time_pars.day_intro = (sim_time_pars.day_intro < 0) ? 0 : sim_time_pars.day_intro;

	return sim_time_pars;

}

ObsSelect SelectObservations(int& day_shut, 
	const std::vector<int>& timeStamps,
	const std::vector<int>& regionalCases,
	const std::vector<int>& regionalDeaths,
	int time_back,
	Utilities::logging_stream::Sptr log) 
{	
	ObsSelect obs_selections;

	// Calculate time under the assumption that time stamps and regional cases data match in size
	const SimTime time_info = GetSimDuration(time_back, regionalCases);
	obs_selections.sim_time = time_info;

	// Add only non-negative case data
	for(int t{1}; t < regionalCases.size(); ++t)
	{
		if(regionalCases[t] >= 0)
		{
			obs_selections.hospitalised.push_back(regionalCases[t]);
			obs_selections.deaths.push_back(regionalDeaths[t]);
		}
	}

	(*log) << "[Observations]:\n";
	(*log) << "    day first report (t_index): " << time_info.t_index << '\n';
		
	// add the extra information on the observations
	int timeSeriesLength = static_cast<int>(obs_selections.hospitalised.size());
	
	// Pad the observations with zeroes up to the duration of the simulation
	// extra_time is the number of days prior the first detection (index cases) which are
	// considered to show silent spread (spread of disease prior detection). The duration of this
	// period is informed by the parameter hrp.
	// Because we are modelling observed/detected cases and deaths (not true presence),
	// we therefore pad the observation time series with an initial run of zeros.
	if (obs_selections.sim_time.duration > timeSeriesLength){
		int extra_time = obs_selections.sim_time.duration - timeSeriesLength;
		obs_selections.hospitalised.insert(obs_selections.hospitalised.begin(), extra_time, 0);
		obs_selections.deaths.insert(obs_selections.deaths.begin(), extra_time, 0);
		day_shut = day_shut + extra_time;
	}
	
	(*log) << "    Number of days of obs cases: " <<
		obs_selections.hospitalised.size() << std::endl;
	(*log) << "    Number of days of obs deaths: " <<
		obs_selections.deaths.size() << std::endl;
	(*log) << "    Number of weeks of obs: " <<
		static_cast<double>(obs_selections.deaths.size()) / 7.0 << std::endl;

	//transform  cumulative numbers into incident cases
	std::vector<int> casesIncidence = ComputeIncidence(obs_selections.hospitalised);
	obs_selections.hospitalised = CorrectIncidence(casesIncidence, obs_selections.hospitalised);	
	
	//transform  cumulative numbers into incident deaths
	std::vector<int> deathsIncidence = ComputeIncidence(obs_selections.deaths);
	obs_selections.deaths = CorrectIncidence(deathsIncidence, obs_selections.deaths);	

	return obs_selections;
}

std::vector<int> ComputeIncidence(const std::vector<int>& timeseries) 
{	
	std::vector<int> incidence(timeseries.size());
	std::adjacent_difference(timeseries.begin(), timeseries.end(), incidence.begin());

	return incidence;
}

std::vector<int> CorrectIncidence(const std::vector<int>& originalIncidence, std::vector<int> timeseries)
{
	
	std::vector<int> correctedIncidence = originalIncidence;

	if ( std::any_of(originalIncidence.begin(), originalIncidence.end(), [](int i){ return i < 0; })) {
		for (unsigned int ii = 1; ii < originalIncidence.size(); ++ii) {
			if (originalIncidence[ii] < 0){
				int case_bef = timeseries[ii - 1];//look next day number of cases
				int case_aft = timeseries[ii + 1];//look at the previous number of cases	
				
      		  	//check if the next tot cases is bigger than tot cases before
    			//if(bigger) compute difference and remove it from day before
		        if (case_bef < case_aft) {
		          if (originalIncidence[ii - 1] >= abs(originalIncidence[ii])) {
		            timeseries[ii - 1] += originalIncidence[ii];
		          } else {
		          	timeseries[ii] = 0;
		          }
		        } else {
		          // if smaller
		          // check when the day prior today was smaller than the next day and remove 
				  // extra cases to this day
		          unsigned int counter_check = 1;
		          while (timeseries[ii - counter_check] > case_aft){
		            ++counter_check ;
		          }

		          if (originalIncidence[ii - counter_check + 1] >= abs(originalIncidence[ii])) {
				  	for ( unsigned int jj = (ii - counter_check + 1); jj < ii; ++jj) {
						timeseries[jj] += originalIncidence[ii];
					}			
		          } else {
		            timeseries[ii] = 0;
		          }
		        }
			}

		}
		
		//recompute the incident cases
		correctedIncidence = ComputeIncidence(timeseries);
	}

	return correctedIncidence;
}

} // namespace Observations
} // namespace EERAModel
