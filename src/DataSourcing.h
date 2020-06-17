#pragma once

#include <string>

#include "ModelTypes.h"
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
        private:
            InputObservations _input_obs;
            ModelInputParameters _input_params;
            PriorParticleParameters _prior_params;
            Utilities::logging_stream::Sptr _log;
        public:
            DataSource(Utilities::logging_stream::Sptr log) : _log(log) {}
            InputObservations getInputObservations() const {return _input_obs;}
            ModelInputParameters getInputParameters() const {return _input_params;}
            PriorParticleParameters getPriorParameters() const {return _prior_params;}
            Utilities::logging_stream::Sptr getLogger() const {return _log;}
            void setInputObservations(const InputObservations& input_obs) {_input_obs = input_obs;}
            void setModelInputParameters(const ModelInputParameters& model_input) {_input_params = model_input;}
            void setPriors(const PriorParticleParameters& priors) {_prior_params = priors;}
            void Setup();
    };

    class LocalSource : public DataSource
    {
       private:
            const DataFiles _data_files;
            void _extract_data();
       public:
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