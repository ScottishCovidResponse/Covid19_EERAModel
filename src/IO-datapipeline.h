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


    IOdatapipeline(const string &params_file);
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

// /**
//  * @brief Determine if data pipeline should be used for fixed parameters
//  * 
//  * @param filePath Path to the INI file
//  * 
//  * @return Boolean value, true, yes use data pipeline, false, no
//  */
// bool UseDatapipelineFixedParameters(const params& paramlist, Utilities::logging_stream::Sptr log);



} // namespace IO
} // namespace EERAModel