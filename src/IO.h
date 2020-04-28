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

} // namespace IO
} // namespace EERAModel