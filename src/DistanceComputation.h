#pragma once

#include <vector>
#include <assert.h>

namespace EERAModel {
namespace DistanceComputation {

double sse_calc_int(std::vector<int> simval, const std::vector<int>& obsval);

} // namespace DistanceComputation
} // namespace EERAModel
