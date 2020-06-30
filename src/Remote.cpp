#include "Remote.h"

namespace EERAModel
{
namespace API
{
    ModelInputParameters ReadParametersFromAPI(const std::string APIInfo, Utilities::logging_stream::Sptr log) {return {};}
    std::vector<double> ReadPosteriorParametersFromAPI(const std::string APIInfo, int set_selection) {return {};}
    InputObservations ReadObservationsFromAPI(const std::string APIInfo, Utilities::logging_stream::Sptr log) {return {};}
};
};