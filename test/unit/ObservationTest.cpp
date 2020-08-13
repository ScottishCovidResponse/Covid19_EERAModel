#include "gtest/gtest.h"

#include "Observations.h"

using namespace EERAModel::Observations;

// Regression test against original code to ensure
// same behaviour
TEST(TestObservations, TestDurationCalculator) {
  std::vector<int> time_back_vals = {10, 0};

  for (int &time_back : time_back_vals) {
    std::cout << "Running Trial: time_back=" << time_back << std::endl;

    // Test based on former code
    const std::vector<int> regionalCases = {1, 0, 1, 10, 11, 11, 23, 34};

    int t_index = -1;
    unsigned int maxTime = 0;

    int duration = 0;

    for (unsigned int t = 1; t < regionalCases.size(); ++t) {

      if (regionalCases[t] >= 0) {
        maxTime = std::max(maxTime, t);
        // identify when the first case is detected/hospitalised
        if (regionalCases[t] > 0 && t_index < 0) {
          t_index = t;
        }
      }
    }

    int intro = t_index + 1 - time_back;

    // define the duration of the study period and day of the incursion
    if (intro < 0) {
      duration = static_cast<int>(maxTime) - intro + 1;
      intro = 0;
    } else {
      duration = maxTime + 1;
      intro = intro;
    }

    SimTime candidate = GetSimDuration(time_back, regionalCases);

    EXPECT_EQ(duration, candidate.duration);
    EXPECT_EQ(intro, candidate.day_intro);
    EXPECT_EQ(t_index, candidate.t_index);
  }
}