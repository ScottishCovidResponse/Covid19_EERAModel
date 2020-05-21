#include "Observations.h"
#include <algorithm>
#include <numeric>

namespace EERAModel {
namespace Observations {

void SelectObservations(int& duration,
	int& day_intro,
	int& day_shut, 
	std::vector<int>& obsHosp,
	std::vector<int>& obsDeaths, 
	const std::vector<int>& timeStamps,
	const std::vector<int>& regionalCases,
	const std::vector<int>& regionalDeaths,
	int time_back) 
{	
	// Reference value of t_index for the first run through the loop - indicates that it hasn't
	// been populated
	int t_index = -1;
	unsigned int maxTime = 0;
	//create the vector of cases (cummulative) and define time of first detection (index)
	for (unsigned int t = 1; t < timeStamps.size(); ++t) {
		
		// This check is required because the input case timeseries may hve negative values - 
		// this is a consequence of the way in which cases are re-allocated to health boards
		if (regionalCases[t] >= 0) {
			maxTime = std::max(maxTime, t);
			
			obsHosp.push_back(regionalCases[t]);
			obsDeaths.push_back(regionalDeaths[t]);
			
			//identify when the first case is detected/hospitalised 
			if (regionalCases[t] > 0 && t_index < 0) {
				t_index = t;
			}		
		}
	}

	// identify the first day of infectiousness for the index case, which will be the start of our
	// simulation, and add days in the observation, and define duration of the disease process
	int intro =  t_index + 1 - time_back;

	// define the duration of the study period and day of the incursion
	if (intro < 0) {
		duration = static_cast<int>(maxTime) - intro + 1;
		day_intro = 0;
	} else {
		duration = maxTime + 1;
		day_intro = intro;
	}
		
	// add the extra information on the observations
	int timeSeriesLength = static_cast<int>(obsHosp.size());
	
	// Pad the observations with zeroes up to the duration of the simulation
	if (duration > timeSeriesLength){
		int extra_time = duration - timeSeriesLength;
		obsHosp.insert(obsHosp.begin(), extra_time, 0);
		obsDeaths.insert(obsDeaths.begin(), extra_time, 0);		
		day_shut = day_shut + extra_time;
	}
	
	//transform  cumulative numbers into incident cases
	std::vector<int> casesIncidence = ComputeIncidence(obsHosp);
	obsHosp = CorrectIncidence(casesIncidence, obsHosp);	
	
	//transform  cumulative numbers into incident deaths
	std::vector<int> deathsIncidence = ComputeIncidence(obsDeaths);
	obsDeaths = CorrectIncidence(deathsIncidence, obsDeaths);	
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
