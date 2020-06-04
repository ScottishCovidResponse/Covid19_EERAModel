#include "gtest/gtest.h"

#include "Utilities.h"

using namespace EERAModel::Utilities;

TEST(TestUtilities, TestSumEveryNth)
{
    const std::vector<int> original = {1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
                                       1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7};
    
    const int weekly_index = 7;

    std::vector<int> reduced = {};

    int weekly_value = 0;
	for(unsigned int ttime =0; ttime< original.size(); ++ttime){
		weekly_value+=original[ttime];
				
		if( (ttime % weekly_index) == weekly_index-1){
			reduced.push_back(weekly_value);
			weekly_value=0;
		}
	}

    std::vector<int> test_reduced = AccumulateEveryNth(original, weekly_index);

    EXPECT_EQ(test_reduced, reduced);

}

TEST(TestUtilities, TestSumSq)
{
    const std::vector<int> vector_sim = {1, 2, 3, 4, 5, 6, 7};
    const std::vector<int> vector_obs = {7, 6, 5, 4, 3, 2, 1};

    const int expected_sum_sq = 2*(pow(6, 2)+pow(4, 2)+pow(2, 2));

    EXPECT_EQ(expected_sum_sq, sse_calc(vector_sim, vector_obs));
}