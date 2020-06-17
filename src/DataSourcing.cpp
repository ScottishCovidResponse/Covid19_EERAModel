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
        _input_params.prior_param_list = _prior_params.prior_param_list;
    }
}

void LocalSource::_extract_data()
{  
    setModelInputParameters(IO::ReadParametersFromFile(_data_files, getLogger()));

    // Read prior particle parameters if run type is "Prediction"
    if(getInputParameters().run_type == ModelModeId::PREDICTION)
    {
        setPriors(IO::ReadPriorParametersFromFile(_data_files, getLogger()));
    }

    setInputObservations(IO::ReadObservationsFromFiles(_data_files, getLogger()));
}

void Remote::_extract_data()
{

    setModelInputParameters(API::ReadParametersFromAPI(_api_info, getLogger()));

    if(getInputParameters().run_type == ModelModeId::PREDICTION)
    {
        setPriors(API::ReadPriorParametersFromAPI(_api_info, getLogger()));
    }

    setInputObservations(API::ReadObservationsFromAPI(_api_info, getLogger()));
    
}

const DataSource getSource(SourceID source_id, Utilities::logging_stream::Sptr log,  std::string data_location)
{
    return (source_id == SourceID::LOCAL) ? DataSource(LocalSource(data_location, log)) : DataSource(Remote(data_location, log));
}

void PushOutputs() {}
};
};