#pragma once

#include "ModelTypes.h"
#include "IniFile.h"
#include "Utilities.h"
#include <sstream>
#include <string>

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
std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, int set_selection);

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
 * @brief Extract a numeric value from an INI file
 * 
 * @param SettingName Name of value to retrieve
 * @param SettingCategory Name of category in which the value is located in the INI file
 * @param filePath Path to the INI file
 * 
 * @return Numeric value (supported types are int or double)
 */
template <typename ParseVariableType>
ParseVariableType ReadNumberFromFile(std::string SettingName, std::string SettingCategory, const std::string& filePath) 
{
	std::string SettingValue = CIniFile::GetValue(SettingName, SettingCategory, filePath);

	char* endptr = nullptr;
	ParseVariableType Value;
	if (std::is_same<ParseVariableType, double>::value) { Value = strtod(SettingValue.c_str(), &endptr); }
	if (std::is_same<ParseVariableType, int>::value) { Value = strtol(SettingValue.c_str(), &endptr, 0); }

	/* Error Handling for strtod or strtol functions*/
	if (Value == -HUGE_VAL || Value == HUGE_VAL || endptr == SettingValue.c_str())
	{
		std::stringstream SettingParseError;
		SettingParseError << std::endl;
		SettingParseError << "Invalid value in Parameter File: " << filePath.c_str() <<  std::endl;
		SettingParseError << "Category: " << SettingCategory.c_str() << std::endl;
		SettingParseError << "Setting: " << SettingName.c_str() << std::endl;
		SettingParseError << "Value: " << SettingValue.c_str() << std::endl;
		throw std::runtime_error(SettingParseError.str());
	}
	return Value;
}

/**
 * @brief Log fixed parameters
 * 
 * Log the fixed parameters of a model
 * 
 * @param params Model input parameters
 * @param log Logger
 */
void LogFixedParameters(const ModelInputParameters& params, Utilities::logging_stream::Sptr log);

} // namespace IO
} // namespace EERAModel
