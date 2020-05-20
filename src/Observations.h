#pragma once

#include <vector>

namespace EERAModel {
namespace Observations {

/**
 * @brief Select observations to use from the input observations
 * 
 */
void select_obs(
	int& duration,
	int& day_intro,
	int& day_shut, 
	std::vector<int>& obsHosp_tmp,
	std::vector<int>& obsDeaths_tmp, 
	const std::vector<int>& timeStamps,
	const std::vector<int>& regionalCases,
	const std::vector<int>& regionalDeaths,
	int time_back);

/**
 * @brief Transform the timeseries of cummulative cases into incidence
 * 
 * A timeseries of cumulative cases represents the number of cases reported to exist at any given
 * time in the timeseries. The incidence represents the change in the number of cases between two
 * successive points in the timeseries. The incidence at any given point in the timeseries is
 * computed as the difference between the cumulative number of cases at that point and the cumulative
 * number of cases at the preceding point. The incidence at the first point in the time series is
 * simply the cumulative number of cases at that first point.
 * 
 * @param timeseries Input timeseries of cumulative instances
 * @param incidence Output computed incidence
 * 
 * @return Nothing
 */
void compute_incidence(const std::vector<int>& timeseries, std::vector<int>& incidence);

/**
 * @brief To do
 */
void correct_incidence(std::vector<int>& v, std::vector<int> cumv);


} // namespace Observations
} // namespace EERAModel
