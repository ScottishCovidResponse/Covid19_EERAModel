#pragma once

#include "ModelTypes.h"
#include "IniFile.h"
#include "Utilities.h"
#include <sstream>
#include <string>
#include <vector>

#include "array.hh"
#include "table.hh"
#include "datapipeline.hh"

namespace EERAModel {
namespace IO {

namespace {
    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        std::cout << std::scientific << std::setprecision(17);

        os << "(" << v.size() << ") [";
        for (auto &ve : v) {
            os << " " << ve;
        }
        os << " ]";
        return os;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& v) {
        os << "(" << v.size() << ") [\n";
        for (auto &ve : v) {
            os << "  " << ve << "\n";
        }
        os << "]";
        return os;
    }
}

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
     * @param model_config Pathname for a directory containing model configuration data
     * @param log_stream An output stream to use for logging diagnostic messages
     * @param dpconfig_path Pathname for a data pipeline configuration yaml file (empty string disables data pipeline)
     */
    IOdatapipeline(string params_path, string model_config, Utilities::logging_stream::Sptr log_stream, string dpconfig_path = "");

    ~IOdatapipeline() {}

    /**
     * @brief Read common model input parameters used by all model types
     * 
     * @return Common model input parameters
     */ 
    CommonModelInputParameters ReadCommonParameters();

    /**
     * @brief Read model observation data used by all model types
     * 
     * @return Observation data
     */ 
    ObservationsForModels ReadModelObservations();

    /**
     * @brief Read Inference config data
     * 
     * @param commonParameters Common parameters, could've been read using ReadCommonParameters above
     * @param modelObservations Model Observations, could've been read using ReadModelObservations above
     * 
     * @return Inference configuration
     */ 
    InferenceConfig ReadInferenceConfig(
        const CommonModelInputParameters& commonParameters, const ObservationsForModels& modelObservations);

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
     * @brief The path of the model data configspecified parameter ".ini" file
     */
    std::string ModelConfigDir;

    /**
     * @private
     * @brief A stream to use for logging diagnostic messages
     */
    Utilities::logging_stream::Sptr log;

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

    /**
     * @brief Read the inference observations
     * 
     * @return Inference observations data structure
     */
    ObservationsForInference ReadInferenceObservations(const ObservationsForModels& modelObservations);

    /**
     * @brief Perform consistency checks on imported data from filePath.
     * 
     * @param data_product The data pipeline data product to load
     * @param component The component within the data product
     * @param axisLength Length of axis to be checked
     * @param expectedLength Expected Length of axis to check against
     * @param axisID String holding axis identifier ("rows" or "columns")
     */
    void ImportConsistencyCheck(
        const std::string& data_product, const std::string& component,
        const unsigned int axisLength, const unsigned int expectedValue, const std::string& axisID);

    /**
     * @brief Read a 2D array from the data pipeline and copy the contents into a vector of vectors
     * 
     * @param data_product The data pipeline data product to load
     * @param component The component within the data product
     * @param result A pointer to the vector of vectors for the result
     * @param expected_rows if >=0, the number of rows expected in the data
     * @param expected_columns, if >=0, the number of columns expected in the data
     */
    template <class T>
    void dparray_to_csv(
        const std::string& data_product, const std::string& component, std::vector<std::vector<T>> *result,
        int expected_rows = -1, int expected_columns = -1)
    {
        (*log) << "\t- (data pipeline) \"" << data_product << "\", \"" << component << "\"" << std::endl;

        Array<double> input = dp->read_array(data_product, component);
        std::vector<int> array_sizes = input.size();

        if (array_sizes.size() != 2) {
            // Should complain about this... and probably should check the dimensions as matching what is expected... if that matters.
        }

        if (expected_rows >= 0) {
            ImportConsistencyCheck(data_product, component, array_sizes[1], expected_rows, "rows");
        }

        if (expected_columns >= 0) {
            ImportConsistencyCheck(data_product, component, array_sizes[0], expected_columns, "columns");
        }

        result->resize(0);

        for (int j = 0; j < array_sizes[1]; ++j) {
            // Construct a new element of isize size
            result->emplace_back(array_sizes[0]);
            auto& row = result->back();

            // Copy the data row
            for (int i = 0; i < array_sizes[0]; ++i) {
                row[i] = static_cast<T> (input(i, j));
            }
        }

        // Some checking
        // std::cout << "Original: \"" << data_product << "\" \"" << component << "\"\n";
        // std::cout << "    size() = [" << array_sizes[0] << ", " << array_sizes[1] << "]\n";
        // std::cout << "    (0,0)=" << input(0,0) << " (1,0)=" << input(1,0) << "\n";
        // std::cout << "    (0,1)=" << input(0,1) << " (1,1)=" << input(1,1) << "\n";

        // std::cout << "VoV:\n";
        // std::cout << "    size() = [" << (*result)[0].size() << ", " << result->size() << "]\n";
        // std::cout << "    (0,0)=" << (*result)[0][0] << " (1,0)=" << (*result)[0][1] << "\n";
        // std::cout << "    (0,1)=" << (*result)[1][0] << " (1,1)=" << (*result)[1][1] << "\n";
    }

    /**
     * @brief Read a table from the data pipeline and copy the contents into a vector of vectors
     * 
     * @param data_product The data pipeline data product to load
     * @param component The component within the data product
     * @param result A pointer to the vector of vectors for the result
     * @param expected_rows if >=0, the number of rows expected in the data
     * @param expected_columns, if >=0, the number of columns expected in the data
     */
    template <class T>
    void dptable_to_csv(
        const std::string& data_product, const std::string& component, std::vector<std::vector<T>> *result,
        int expected_rows = -1, int expected_columns = -1)
    {
        (*log) << "\t- (data pipeline) \"" << data_product << "\", \"" << component << "\"" << std::endl;

        Table input = dp->read_table(data_product, component);

        // std::cout << input.to_string() << "\n";

        std::vector<string> columns = input.get_column_names();
        std::size_t col_size = input.get_column_size();

        if (expected_rows >= 0) {
            ImportConsistencyCheck(data_product, component, col_size, expected_rows, "rows");
        }

        if (expected_columns >= 0) {
            ImportConsistencyCheck(data_product, component, columns.size(), expected_columns, "columns");
        }

        result->resize(0);
        result->resize(col_size);

        for (const auto& column : columns) {
            std::vector<T>& column_values = input.get_column<T>(column);

            for (std::size_t e = 0; e < col_size; ++e) {
                (*result)[e].push_back(column_values[e]);
            }
        }
    }

    void dpdistribution(
        const std::string& data_product, const std::string& component,
        std::string p1, double *a, std::string p2 = "", double *b = nullptr);
};


} // namespace IO
} // namespace EERAModel