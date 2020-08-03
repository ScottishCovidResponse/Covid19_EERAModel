#include "gtest/gtest.h"
#include "Random.h"
#include "ModelCommon.h"

using namespace EERAModel::Model;
using namespace EERAModel::Random;

TEST(TestFlowFunction, TestNotNegative)
{
    const int n_trials = 10;

    RNG::Sptr rng = std::make_shared<RNG>(time(nullptr));

    for(int i{0}; i < n_trials; ++i)
    {
        const double rate = pow(rng->Flat(1, 9), -1*rng->Flat(1, 7));
        const int orig_val = static_cast<int>(rng->Flat(0, 1E5));
        const int new_val = static_cast<int>(rng->Flat(0, orig_val)); 

        std::cout << "Running Trial: " << i << "/" << n_trials << " : ";
        std::cout << "Flow(rng, " << orig_val << ", " << new_val << ", ";
        std::cout << rate << ")" << std::endl;

        EXPECT_TRUE(Flow(rng, orig_val, new_val, rate) >= 0);
    }
}