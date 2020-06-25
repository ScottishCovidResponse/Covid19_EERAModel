#pragma once

#include "ModelTypes.h"

#include <string>
#include "Utilities.h"
#include "IniFile.h"

#ifndef ROOT_DIR
#error Macro ROOT_DIR must be defined!
#endif

namespace EERAModel {
/**
 * @brief Namespace containing functions related to input/output of data
 *
 * This namespace contains functions which handle the reading of data from
 * various input files.
 */
namespace IO {

/**
 * @brief Read model input parameters from an INI file
 * 
 * @param filePath Path to the INI file
 * 
 * @return Model parameters
 */
ModelInputParameters ReadParametersFromFile(const std::string& filePath, const Utilities::logging_stream::Sptr& log);

/**
 * @brief Read model posterior parameters from a CSV file
 * 
 * @param filePath Path to CSV file
 * @param set_selection Selection of row in CSV file for posterior parameters
 * 
 * @return Model posterior parameters
 */
std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, const int& set_selection);

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

/**
 * @brief Writes Prediction outputs to files
 * 
 * @param status Status object
 * @param end_comps Matrix to hold the end states of the simulation organised by compartments
 * @param outDirPath Path to output directory where output files will be stored
 */
void WritePredictionsToFiles(Status status, std::vector<std::vector<int>>& end_comps, 
	const std::string& outDirPath, const Utilities::logging_stream::Sptr& log);

/**
 * @brief String to number conversion
 * 
 * @param SettingName Name of setting to retrieve
 * @param SettingCategory Name of category for the setting in the INI file
 * @param filePath File path of the INI file to be read
 * 
 * @return ParseVariableType, which is either int or double
 */
template <typename ParseVariableType>
ParseVariableType StrToNumber(std::string SettingName, std::string SettingCategory, const std::string& filePath) noexcept(false);

} // namespace IO
} // namespace EERAModel
