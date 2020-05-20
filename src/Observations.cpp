#include "Observations.h"
#include <algorithm>
#include <numeric>

namespace EERAModel {
namespace Observations {

void select_obs(int& duration,
	int& day_intro,
	int& day_shut, 
	std::vector<int>& obsHosp_tmp,
	std::vector<int>& obsDeaths_tmp, 
	const std::vector<int>& timeStamps,
	const std::vector<int>& regionalCases,
	const std::vector<int>& regionalDeaths,
	int time_back) {
			
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
			
			obsHosp_tmp.push_back(regionalCases[t]);
			obsDeaths_tmp.push_back(regionalDeaths[t]);
			
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
	int timeSeriesLength = static_cast<int>(obsHosp_tmp.size());
	
	// Pad the observations with zeroes up to the duration of the simulation
	if(duration > timeSeriesLength){
		int extra_time = duration - timeSeriesLength;
		obsHosp_tmp.insert(obsHosp_tmp.begin(), extra_time, 0);
		obsDeaths_tmp.insert(obsDeaths_tmp.begin(), extra_time, 0);		
		day_shut = day_shut + extra_time;
	}
	
	//transform  cumulative numbers into incident cases
	std::vector<int> obsHosp_tmp2(obsHosp_tmp.size());
	Observations::compute_incidence(obsHosp_tmp,obsHosp_tmp2);
	Observations::correct_incidence(obsHosp_tmp2,obsHosp_tmp);	
	obsHosp_tmp = obsHosp_tmp2;
	
	//transform  cumulative numbers into incident deaths
	std::vector<int> obsDeaths_tmp2(obsDeaths_tmp.size());
	Observations::compute_incidence(obsDeaths_tmp,obsDeaths_tmp2);
	Observations::correct_incidence(obsDeaths_tmp2,obsDeaths_tmp);	
	obsDeaths_tmp = obsDeaths_tmp2;
}

void compute_incidence(const std::vector<int>& timeseries, std::vector<int>& incidence) {
	std::adjacent_difference(timeseries.begin(), timeseries.end(), incidence.begin());
}

//correct the timeseries to avoid negative incidence records
void correct_incidence(std::vector<int>& v, std::vector<int> cumv){
	if ( std::any_of(v.begin(), v.end(), [](int i){return i<0;}) ){
		for (unsigned int ii = 1; ii < v.size(); ++ii) {
			if(v[ii]<0){
				int case_bef = cumv[ii-1];//look next day number of cases
				int case_aft = cumv[ii+1];//look at the previous number of cases	
				
      		  	//check if the next tot cases is bigger than tot cases before
    			//if(bigger) compute difference and remove it from day before
		        if(case_bef<case_aft){
		          if(v[ii-1] >= abs(v[ii])){
		            cumv[ii-1] += v[ii];
		          } else {
		          	cumv[ii]=0;
		          }
		        } else{
		          //if smaller
		          //check when the day prior today was smaller than the next day and remove extra cases to this day
		          unsigned int counter_check = 1;
		          while( cumv[ii-counter_check]>case_aft){
		            ++counter_check ;
		          }
		          if(v[ii-counter_check+1] >= abs(v[ii])){
				  	for ( unsigned int jj = (ii-counter_check+1); jj < (ii); ++jj) {
						cumv[jj] += v[ii];
					}			
		          } else {
		            cumv[ii]=0;
		          }
		        }
			}

		}
		//recompute the incident cases
		v.clear();
		compute_incidence(cumv,v);
	}
}

} // namespace Observations
} // namespace EERAModel
