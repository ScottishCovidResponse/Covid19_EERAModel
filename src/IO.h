#pragma once

#include "ModelTypes.h"

#include <string>
#include "Utilities.h"

#ifndef ROOT_DIR
#error Macro ROOT_DIR must be defined!
#endif

namespace EERAModel {
namespace IO {

/**
 * @brief Read model input parameters from an INI file
 * 
 * @param filePath Path to the INI file
 * 
 * @return Model parameters
 */
ModelInputParameters ReadParametersFromFile(const std::string& filePath, const Utilities::logging_stream::Sptr& log);


std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, const int set_selection, const Utilities::logging_stream::Sptr& log);

/**
 * @brief Read in observations
 * 
 * @param dirPath Path to the directory containing the observation files
 * 
 * @return Observations
 */
InputObservations ReadObservationsFromFiles(const Utilities::logging_stream::Sptr& log);

/**
 * @brief  Write outputs to files
 * 
 * @param smc Iteration number
 * @param herd_id Herd Id
 * @param Nparticle Number of particles
 * @param nPar Number of parameters
 * @param particleList Vector of particles to write out
 * @param outDirPath Path to the directory in which the output files should be placed
 */
void WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, 
	const std::vector<EERAModel::particle>& particleList, const std::string& outDirPath,
	const Utilities::logging_stream::Sptr& log);

void WritePredictionsToFiles(Status status, std::vector<std::vector<int>>& end_comps, 
	const std::string& outDirPath, const Utilities::logging_stream::Sptr& log);

/**
 * @brief Convert Vector of Compartments struct to a vector of integers
 * 
 * NOTE: This is a temporary function to allow compatibility
 * converts the Compartments struct to a vector of integers
 * 
 * @param cmps_vec Vector of compartments struct containing population per category
 * 
 * @return Vector of population counters
 */
std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec);

} // namespace IO
} // namespace EERAModel
