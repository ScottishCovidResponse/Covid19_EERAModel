#pragma once

#include <string>

namespace EERAModel
{
/**
 * @brief Namespace containing data sourcing objects
 * 
 * This namespace contains definitions for data sources which
 * define the location of the various data files.
 */
namespace DataSourcing
{
    /**
     * @brief Structure containing file addresses
     * 
     * All file locations for the given data source form
     * be it local or remote are defined in this structure.
     */
    struct DataFiles
    {
        std::string parameters; 
        std::string data;
        std::string deaths;
        std::string ages;

        struct WAIFW
        {
            std::string norm;
            std::string home;
            std::string sdist;
        };
        WAIFW waifw;
        std::string cfr_byage;
        std::string frail;
        std::string prior_params;
    };

    /**
     * @brief Core data source class
     * 
     * Class defining locations of the various data inputs
     */
    class DataSource
    {
        public:
            const DataFiles data_files;
    };

    /**
     * @brief Local file source initialisation
     * 
     * Function constructs a local file source file hierarchy
     * structure for a given root address
     *
     * @param root_dir The root data directory within which the expected file structure is found
     *
     * @return Data source object with the given file definitions
     */
    DataSource Local(const std::string root_dir=std::string(ROOT_DIR));

    /**
     * @brief Remote file source initialisation
     */
    const DataSource Remote({});

    /**
     * @brief Selection of local or remote file source
     */
    enum class SourceID
    {
        REMOTE, /*!< Files are to be read through an API. */
        LOCAL   /*!< Files are to be read from local directory. */
    };

    /**
     * @brief Return the relevant source object
     * 
     * Given the specified options return the relevant data source object
     * for performing the model run.
     * 
     * @param source_id whether the data is remote or local
     * @param local_data_location if the data is local, the data root directory
     * 
     * @return an appropriate data source object for the options specified.
     */
    const DataSource getSource(SourceID source_id, const std::string local_data_location="");

};
};