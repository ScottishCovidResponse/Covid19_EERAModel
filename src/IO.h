#pragma once

#include "ModelTypes.h"

#include <string>

namespace EERAModel {
namespace IO {

/**
 * @brief Read model input parameters from an INI file
 * 
 * @param filePath Path to the INI file
 * 
 * @return Model parameters
 */
EERAModel::ModelInputParameters ReadParametersFromFile(const std::string& filePath);

/**
 * @brief Read in observations
 * 
 * @param dirPath Path to the directory containing the observation files
 * 
 * @return Observations
 */
EERAModel::Observations ReadObservationsFromFiles();

/**
 * @brief  Write outputs to files
 * 
 * @param smc Iteration number
 * @param herd_id Herd Id
 * @param Nparticle Number of particles
 * @param nPar Number of parameters
 * @param particleList Vector of particles to write out
 */
void WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, const std::vector<EERAModel::particle>& particleList);

} // namespace IO
} // namespace EERAModel