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

const DataSource getSource(SourceID source_id, Utilities::logging_stream::Sptr log,  std::string local_data_location)
{
    return (source_id == SourceID::LOCAL) ? LocalSource(local_data_location, log) : Remote;
}
};
};