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

class IOException: public std::exception
{
public:
    IOException(const std::string& message) : message_(message) {}
    
    const char* what() const throw() 
    {
        return message_.c_str();
    }

private:
    std::string message_;
};

/**
 * @brief Read supplementary parameters used for main.cpp.
 * 
 * @param ParamsPath Path to INI file
 * @param log Logger
 * 
 * @return Supplementary parameters
 */
SupplementaryInputParameters ReadSupplementaryParameters(const std::string& ParamsPath,
    Utilities::logging_stream::Sptr log);

/**
 * @brief Read common model input parameters used by all model types
 * 
 * @param ParamsPath Path to INI file
 * 
 * @return Common model input parameters
 */ 
CommonModelInputParameters ReadCommonParameters(const std::string& ParamsPath);

/**
 * @brief Read parameters for the inference mode
 * 
 * @param configDir Directory containing the configuration and data files
 * @param log Logger
 * 
 * @return Inference parameters
 */
InferenceConfig ReadInferenceConfig(const std::string& configDir, Utilities::logging_stream::Sptr log, const CommonModelInputParameters& commonParameters);

/**
 * @brief Read prediction framework configuration from input files
 * 
 * @param configDir Directory containing the configuration and data files
 * @param index Index of the parameter set to select from the posterior parameters file
 * @param log Logger
 * 
 * @return Prediction configuration
 */
PredictionConfig ReadPredictionConfig(const std::string& configDir, int index, Utilities::logging_stream::Sptr log, const CommonModelInputParameters& commonParameters);

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
 * @brief Read model parameters from a CSV file
 * 
 * @param filePath Path to CSV file
 * @param set_selection Selection of row in CSV file for parameters
 * 
 * @return Model parameters
 */
std::vector<double> ReadPredictionParametersFromFile(const std::string& filePath, int set_selection);

/**
 * @brief Read observations needed for Inference framework
 * 
 * @param configDir Directory containing the data files
 * @param log Logger
 * 
 * @return Observations needed for Inference framework
 */
ObservationsForInference ReadInferenceObservations(const std::string& configDir, Utilities::logging_stream::Sptr log);

/**
 * @brief Read observations needed for all models
 * 
 * @param configDir Directory containing the data files
 * @param log Logger
 * 
 * @return Observations needed for all models
 */
ObservationsForModels ReadModelObservations(const std::string& configDir, Utilities::logging_stream::Sptr log);

/**
 * @brief Read seed settings from the parameters file
 * 
 * @param ParamsPath Path to INI file
 * @param log Logger
 * 
 * @return Seed settings data structure
 */
seed ReadSeedSettings(const std::string& ParamsPath, Utilities::logging_stream::Sptr log);

/**
 * @brief Read the fixed model parameters from the parameters file
 * 
 * @param ParamsPath Path to INI file
 * 
 * @return Parameters data structure
 */
params ReadFixedModelParameters(const std::string& ParamsPath);

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
 * @param log Logger
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
 * @brief Extract a Boolean value from an INI file
 * 
 * @param SettingName Name of value to retrieve
 * @param SettingCategory Name of category in which the value is located in the INI file
 * @param filePath Path to the INI file
 * 
 * @return Boolean value
 */
bool ReadBoolFromFile(std::string SettingName, std::string SettingCategory, const std::string& filePath);

/**
 * @brief Extract a String value from an INI file
 * 
 * @param SettingName Name of value to retrieve
 * @param SettingCategory Name of category in which the value is located in the INI file
 * @param filePath Path to the INI file
 * 
 * @return String value
 */
std::string ReadStringFromFile(std::string SettingName, std::string SettingCategory, const std::string& filePath);

/**
 * @brief Confirm the existence of value record in an INI file
 * 
 * @param SettingName Name of value to retrieve
 * @param SettingCategory Name of category in which the value is located in the INI file
 * @param filePath Path to the INI file
 * 
 * @return bool indicating presence or not
 */
bool ExistsInFile(std::string SettingName, std::string SettingCategory, const std::string& filePath);

/**
 * @brief Log fixed parameters
 * 
 * Log the fixed parameters of a model
 * 
 * @param paramlist Model fixed parameters
 * @param log Logger
 */
void LogFixedParameters(const params& paramlist, Utilities::logging_stream::Sptr log);

/**
 * @brief Log randomiser settings
 * 
 * Log the settings of the configured randomiser
 * 
 * @param params Supplementary input parameters
 * @param log Logger
 */
void LogRandomiserSettings(const SupplementaryInputParameters& params, unsigned long randomiser_seed, 
    Utilities::logging_stream::Sptr log);

/**
 * @brief Log model seed settings
 * 
 * Log the settings of the model disease seeding
 * 
 * @param params Model input parameters
 * @param log Logger
 */
void LogSeedSettings(const seed& params, Utilities::logging_stream::Sptr log);

/**
 * @brief Log prediction configuration
 * 
 * Log the settings of a prediction run
 * 
 * @param params Model input parameters
 * @param log Logger
 */
void LogPredictionConfig(const PredictionConfig& config, Utilities::logging_stream::Sptr log);

/**
 * @brief Log Git repository version information
 * 
 * @param log Logger
 */
void LogGitVersionInfo(Utilities::logging_stream::Sptr log);

/**
 * @brief Write formatted fixed parameters to an output logging stream
 * 
 * @param log Output logger
 * @param paramlist Fixed parameters
 * 
 * @return Output logger
 */
void OutputFixedParameters(Utilities::logging_stream::Sptr& log, const params& paramlist);

} // namespace IO
} // namespace EERAModel
