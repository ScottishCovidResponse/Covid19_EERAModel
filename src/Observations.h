#pragma once

#include "Utilities.h"
#include <vector>

namespace EERAModel {
namespace Observations {

/**
 * @brief Select observations to use from the regional input observations
 * 
 * Bsaed on the regional cases and deaths timeseries, this function computes the incidence of cases
 * and deaths, correcting for nrgative incidence records.
 * 
 * @param duration ?
 * @param day_into Determined first day of infection
 * @param day_shut ?
 * @param obsHosp Case incidence timeseries, corrected for negative incidences
 * @param obsDeaths Deaths incidence timeseries, corrected for negative incidences
 * @param timeStamps Reference times for other time series
 * @param regionalCases Observed timeseries of cases in the region
 * @param regionalDeaths Observed timeseries of deaths in the region
 * @param time_back ?
 * @param log Logger handle
 */
void SelectObservations(
	int& duration,
	int& day_intro,
	int& day_shut, 
	std::vector<int>& obsHosp,
	std::vector<int>& obsDeaths, 
	const std::vector<int>& timeStamps,
	const std::vector<int>& regionalCases,
	const std::vector<int>& regionalDeaths,
	int time_back,
	Utilities::logging_stream::Sptr log);

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
 * 
 * @return Computed incidence timeseries
 */
std::vector<int> ComputeIncidence(const std::vector<int>& timeseries);

/**
 * @brief Correct the timeseries to avoid negative incidence records
 * 
 * Negative incidence records occur where the value of a timeseries decreases between two successive
 * points in time
 * 
 * @todo Fill in with a detailed description of the algorithm
 * 
 * @param originalIncidence Original computation of incidence, with potential negative records
 * @param timeseries Timeseries of cases from which the original incidence was computed
 * 
 * @return Corrected incidence
 */
std::vector<int> CorrectIncidence(const std::vector<int>& originalIncidence, std::vector<int> timeseries);


} // namespace Observations
} // namespace EERAModel
