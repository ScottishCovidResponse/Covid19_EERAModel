#pragma once

#include <string>

#include "ModelTypes.h"
#include "LocalFileStructure.h"
#include "IO.h"

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
     * @brief Core data source class
     * 
     * Class defining locations of the various data inputs
     */
    class DataSource
    {
        private:
            InputObservations _input_obs;
            ModelInputParameters _input_params;
            PriorParticleParameters _prior_params;
            Utilities::logging_stream::Sptr _log;
        public:
            /**
             * @brief Constructor for Data Source object
             * 
             * @param log shared pointer for logger
             */
            DataSource(Utilities::logging_stream::Sptr log) : _log(log) {}

            /**
             * @brief Fetch input observations
             * 
             * @return current read input observations
             */
            InputObservations getInputObservations() const {return _input_obs;}

            /**
             * @brief Fetch input parameters
             * 
             * @return current read input parameters
             */
            ModelInputParameters getInputParameters() const {return _input_params;}

            /**
             * @brief Fetch prior parameters
             * 
             * @return current prior parameters
             */
            PriorParticleParameters getPriorParameters() const {return _prior_params;}

            /**
             * @brief Fetch current logger
             * 
             * @return shared pointer to current logger
             */
            Utilities::logging_stream::Sptr getLogger() const {return _log;}

            /**
             * @brief Set the current input observations
             */
            void setInputObservations(const InputObservations& input_obs) {_input_obs = input_obs;}

            /**
             * @brief Set the current input parameters
             */
            void setModelInputParameters(const ModelInputParameters& model_input) {_input_params = model_input;}

            /**
             * @brief Set the current prior parameters
             */
            void setPriors(const PriorParticleParameters& priors) {_prior_params = priors;}

            /**
             * @brief Run global setup
             * 
             * Runs setup/initialisation that is relevant to all model runs independent
             * of the data source location. Must be called at appropriate time within each subclass.
             */
            void Setup();
    };

    /**
     * @brief Local data source sub-class
     * 
     * This subclass of DataSource handles the procedure for processing
     * and importing data locally.
     */
    class LocalSource : public DataSource
    {
       private:
            const DataFiles _data_files;
            void _extract_data();
       public:
            /**
             * @brief Constructor for Local Data Source object
             * 
             * @param root_dir the root directory of the data files
             * @param log shared pointer for logger
             */
            LocalSource(const std::string root_dir, Utilities::logging_stream::Sptr log) : 
                DataSource(log),
                _data_files({
                    root_dir+"/parameters.ini",
                    root_dir+"/scot_data.csv",
                    root_dir+"/scot_deaths.csv",
                    root_dir+"/scot_age.csv",
                    {root_dir+"/waifw_norm.csv",
                        root_dir+"/waifw_home.csv",
                        root_dir+"/waifw_sdist.csv"},
                    root_dir+"/cfr_byage.csv",
                    root_dir+"/scot_frail.csv",
                    root_dir+"/src/prior_particle_params.csv"})
            {
                _extract_data();
                (*log) << "[Parameters File]:\n    " << _data_files.parameters << std::endl;
                Setup();
            }

            /**
             * @brief Fetch local data files
             * 
             * @return struct containing local data file locations
             */
            DataFiles getDataFiles() const {return _data_files;}
    };

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
    const DataSource getSource(SourceID source_id, Utilities::logging_stream::Sptr log, std::string local_data_location="");

};
};