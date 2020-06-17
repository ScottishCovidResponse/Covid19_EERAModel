#pragma once

#include <string>

#include "ModelTypes.h"
#include "Utilities.h"

namespace EERAModel
{
namespace API
{
    /**
     * @brief Fetch parameter from API
     * 
     * @param param_name name of parameter
     * @param version version of parameter
     */
    template<typename T>
    T ReadParameterFromAPI(const std::string APIInfo, const std::string param_name, const std::string version="latest") {}
    ModelInputParameters ReadParametersFromAPI(const std::string APIInfo, Utilities::logging_stream::Sptr log);
    PriorParticleParameters ReadPriorParametersFromAPI(const std::string APIInfo, Utilities::logging_stream::Sptr log);
    InputObservations ReadObservationsFromAPI(const std::string APIInfo, Utilities::logging_stream::Sptr log);

    void PushOutputs();
};
};