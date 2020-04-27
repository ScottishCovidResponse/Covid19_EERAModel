#pragma once

#include <vector>

namespace EERAModel {
namespace DistanceComputation {

double sse_calc(int Npop, std::vector<double> simHosp, const std::vector<int>& obsHosp);

double sse_calc_int(std::vector<int> simval, std::vector<int> obsval);

double nsse_calc_int(std::vector<int> simHosp, const std::vector<int>& obsHosp);

double sst_calc(int Npop, std::vector<double> simHosp, const std::vector<int>& obsHosp);

double sst_calc_int(int Npop, std::vector<int> simHosp, const std::vector<int>& obsHosp);

} // namespace DistanceComputation
} // namespace EERAModel