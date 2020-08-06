#pragma once

#include "ModelTypes.h"
#include "IniFile.h"
#include "Utilities.h"
#include <sstream>
#include <string>

#include "datapipeline.hh"

namespace EERAModel {
namespace IO {

/**
 * @class IOdatapipeline
 * @brief IO class to manage access to the data pipeline API
 */
class IOdatapipeline
{
public:
    /* Exclude default constructors and operators */
    IOdatapipeline() = delete;
    IOdatapipeline(const IOdatapipeline& other) = delete;
    IOdatapipeline(IOdatapipeline&& other) = delete;
    IOdatapipeline& operator=(const IOdatapipeline& other) = delete;
    IOdatapipeline& operator=(IOdatapipeline&& other) = delete;

    /**
     * @brief Constructor
     * 
     * Constructs an object suitable for interacting with the data pipeline local
     * store. Note that a local python object "pybind11::scoped_interpreter guard{};"
     * must be active and remain in scope for the durection of the use of this class.
     * 
     * @param params_path Pathname for a parameters ".ini" file
     * @param dpconfig_path Pathname for a data pipeline configuration yaml file (empty string disables data pipeline)
     */
    IOdatapipeline(string params_path, string dpconfig_path = "");

    ~IOdatapipeline() {}

    /**
     * @brief Read common model input parameters used by all model types
     * 
     * @return Common model input parameters
     */ 
    CommonModelInputParameters ReadCommonParameters();

private:
    /**
     * @private
     * @brief Python interpreter, will live for the duration of the class. Not sure if this is correct
     */
    // pybind11::scoped_interpreter guard;

    /**
     * @private
     * @brief The path of the specified parameter ".ini" file
     */
    std::string ParamsPath;

    /**
     * @private
     * @brief Flag indicating if the data pipeline has been requested and open successfully 
     */
    bool datapipelineActive = false;

    /**
     * @private
     * @brief Handle for C++ data pipeline API
     */    
    std::unique_ptr<DataPipeline> dp;

private:
    /**
     * @brief Read the fixed model parameters from the parameters file
     * 
     * @return Parameters data structure
     */
    params ReadFixedModelParameters();
};


} // namespace IO
} // namespace EERAModel