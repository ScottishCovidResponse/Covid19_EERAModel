#include "DataSourcing.h"

namespace EERAModel
{
namespace DataSourcing
{
    void DataSource::Setup()
    {
        // Read prior particle parameters if run type is "Prediction"
        if(_input_params.run_type == ModelModeId::PREDICTION)
        {
            _input_params.posterior_param_list = _posterior_params;
        }
    }

    void LocalSource::_extract_data()
    {  
        setModelInputParameters(IO::ReadParametersFromFile(_data_files, getLogger()));

        // Read prior particle parameters if run type is "Prediction"
        if(getInputParameters().run_type == ModelModeId::PREDICTION)
        {
            setPosteriors(IO::ReadPosteriorParametersFromFile(_data_files, getInputParameters().posterior_parameter_select));
        }

        setInputObservations(IO::ReadObservationsFromFiles(_data_files, getLogger()));
    }

    void RemoteSource::_extract_data()
    {

        setModelInputParameters(API::ReadParametersFromAPI(_api_info, getLogger()));

        if(getInputParameters().run_type == ModelModeId::PREDICTION)
        {
            setPosteriors(API::ReadPosteriorParametersFromAPI(_api_info, getInputParameters().posterior_parameter_select));
        }

        setInputObservations(API::ReadObservationsFromAPI(_api_info, getLogger()));
        
    }

    const DataSource getSource(SourceID source_id, Utilities::logging_stream::Sptr log,  std::string data_location)
    {
        return (source_id == SourceID::LOCAL) ? DataSource(LocalSource(data_location, log)) : DataSource(RemoteSource(data_location, log));
    }

};
};