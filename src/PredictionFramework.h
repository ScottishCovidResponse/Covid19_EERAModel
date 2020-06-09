#pragma once

#include <vector>
#include <functional>
#include "ModelTypes.h"
#include "Random.h"
#include "Model.h"

namespace EERAModel {
namespace Prediction {

class Framework
{
public:
    using ModelRunner = std::function<
        Status (
            std::vector<double>,
            const std::vector<EERAModel::params>&,
            AgeGroupData,
            seed,
            int,
            std::vector<int>,
            int,
            ModelStructureId,
            Random::RNGInterface::Sptr
        )
    >;
    
    Framework(ModelRunner modelRunner,
        ModelStructureId modelStructure,
        unsigned int size,
        AgeGroupData ageGroupData,
        const ModelInputParameters& modelInputParameters,
        int dayShut,
        std::vector<int> ageNums,
        Random::RNGInterface::Sptr rng);

    void RunModel(EERAModel::params parameterSet, int n_sim_steps);

private:
    ModelRunner modelRunner_;

    std::vector<EERAModel::params> fixed_parameters_;

    AgeGroupData ageGroupData_;

    seed seedlist_;

    int dayShut_;

    std::vector<int> ageNums_;

    ModelStructureId modelStructure_;

    Random::RNGInterface::Sptr rng_;
};

} // namespace Prediction
} // namespace EERAModel