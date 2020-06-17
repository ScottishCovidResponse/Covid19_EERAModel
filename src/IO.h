#pragma once

#include "IniFile.h"
#include "ModelTypes.h"
#include "LocalFileStructure.h"

#include <string>
#include <sstream>
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
EERAModel::ModelInputParameters ReadParametersFromFile(const DataSourcing::DataFiles& data_files, const Utilities::logging_stream::Sptr& log);


EERAModel::PriorParticleParameters ReadPriorParametersFromFile(const DataSourcing::DataFiles& data_files, const Utilities::logging_stream::Sptr& log);

/**
 * @brief Read in observations
 * 
 * @param dirPath Path to the directory containing the observation files
 * 
 * @return Observations
 */
EERAModel::InputObservations ReadObservationsFromFiles(const DataSourcing::DataFiles& data_files, const Utilities::logging_stream::Sptr& log);

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

} // namespace IO
} // namespace EERAModel
