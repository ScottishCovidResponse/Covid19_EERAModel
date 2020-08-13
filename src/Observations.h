#pragma once

#include "Utilities.h"
#include <vector>

namespace EERAModel {
/**
 * @brief Namespace holding objects and functions relating to observations
 *
 * This namespace contains structs and functions which handle the data
 * read from data files, adjusting for different input structures.
 */
namespace Observations {

/**
 * @brief Structure to hold deduced simulation time parameters
 */
struct SimTime {
  int duration = 0;
  int day_intro = 0;
  int t_index = 0;
};

/**
 * @brief Structure to hold observation selections
 *
 * Contains information on data collection timing and number of cases/deaths
 */
struct ObsSelect {
  SimTime sim_time;
  std::vector<int> deaths;
  std::vector<int> hospitalised;
};

/**
 * @brief Calculates the duration of the current simulation
 *
 * Uses the regional cases data to deduce the time structure of the
 * observations. As the time series data may have negative incident events the
 * function first determines the non-negative start and end points. Accounting
 * for this is vital in terms of future-proofing the inference.
 *
 * @param time_back offset in start time (relative to first non-negative event)
 * @param regional_cases vector of time series data for a given region
 *
 * @return Struct containing the simulation time parameters
 */
SimTime GetSimDuration(int time_back, const std::vector<int> &regionalCases);

/**
 * @brief Select observations to use from the regional input observations
 *
 * Bsaed on the regional cases and deaths timeseries, this function computes the
 * incidence of cases and deaths, correcting for negative incidence records.
 *
 * @param day_shut ?
 * @param timeStamps Reference times for other time series
 * @param regionalCases Observed timeseries of cases in the region
 * @param regionalDeaths Observed timeseries of deaths in the region
 * @param time_back ?
 * @param log Logger handle
 *
 * @return Struct containing deduced observation selections
 */
ObsSelect SelectObservations(int &day_shut, const std::vector<int> &timeStamps,
                             const std::vector<int> &regionalCases,
                             const std::vector<int> &regionalDeaths,
                             int time_back,
                             Utilities::logging_stream::Sptr log);

/**
 * @brief Transform the timeseries of cumulative cases into incidence
 *
 * A timeseries of cumulative cases represents the number of cases reported to
 * exist at any given time in the timeseries. The incidence represents the
 * change in the number of cases between two successive points in the
 * timeseries. The incidence at any given point in the timeseries is computed as
 * the difference between the cumulative number of cases at that point and the
 * cumulative number of cases at the preceding point. The incidence at the first
 * point in the time series is simply the cumulative number of cases at that
 * first point.
 *
 * @param timeseries Input timeseries of cumulative instances
 *
 * @return Computed incidence timeseries
 */
std::vector<int> ComputeIncidence(const std::vector<int> &timeseries);

/**
 * @brief Correct the timeseries to avoid negative incidence records
 *
 * Negative incidence records occur where the value of a timeseries decreases
 * between two successive points in time
 *
 * @todo Fill in with a detailed description of the algorithm
 *
 * @param originalIncidence Original computation of incidence, with potential
 * negative records
 * @param timeseries Timeseries of cases from which the original incidence was
 * computed
 *
 * @return Corrected incidence
 */
std::vector<int> CorrectIncidence(const std::vector<int> &originalIncidence,
                                  std::vector<int> timeseries);

} // namespace Observations
} // namespace EERAModel
