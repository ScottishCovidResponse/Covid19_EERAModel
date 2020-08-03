#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Utilities.h"

using namespace EERAModel::Utilities;
using namespace testing;

TEST(TestUtilities, TestSummEveryNReturnsEmptyIfVectorSizeLessThanN) {
  std::vector<int> data(3, 0);

  auto actual = AccumulateEveryN(data, 4);

  ASSERT_EQ(actual.size(), 0);
}

TEST(TestUtilities, TestSumEveryNReturnsEmptyIfDataIsEmpty) {
  std::vector<int> data;

  auto actual = AccumulateEveryN(data, 1);

  ASSERT_EQ(actual.size(), 0);
}

TEST(TestUtilities, TestSumEveryNReturnsSingleValueIfNEqualsDataSize) {
  std::vector<int> data{0, 1, 2, 3, 4, 5};

  auto actual = AccumulateEveryN(data, 6);

  ASSERT_THAT(actual, ElementsAre(15));
}

TEST(TestUtilities, TestSumEveryNDiscardsResidualData) {
  std::vector<int> data{0, 1, 2, 3, 4, 5, 6};

  auto actual = AccumulateEveryN(data, 2);

  ASSERT_THAT(actual, ElementsAre(1, 5, 9));
}

TEST(TestUtilities, TestSumEveryN) {
  const std::vector<int> data = {1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
                                 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7};

  auto actual = AccumulateEveryN(data, 7);

  ASSERT_THAT(actual, ElementsAre(28, 28, 28, 28));
}

TEST(TestUtilities, TestSumSq) {
  const std::vector<int> vector_sim = {1, 2, 3, 4, 5, 6, 7};
  const std::vector<int> vector_obs = {7, 6, 5, 4, 3, 2, 1};

  const int expected_sum_sq = 2 * (pow(6, 2) + pow(4, 2) + pow(2, 2));

  EXPECT_EQ(expected_sum_sq, sse_calc(vector_sim, vector_obs));
}