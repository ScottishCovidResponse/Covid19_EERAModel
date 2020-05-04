#pragma once

#include <vector>

namespace EERAModel {
namespace Utilities {

double mean_calc_double(std::vector<double> v);

void cumsum_calc(std::vector<double> v, std::vector<double>& r_val);

double max_calc(double a,double b);

double sum_calc(std::vector<double> v);

double colsum_calc(std::vector<std::vector<double>> m, int colid);

double mean_calc_int(std::vector<int> v);

void cumsum_calc_int(std::vector<int> v, std::vector<int>& r_val);

int max_calc_int(int a,int b);

int sum_calc_int(std::vector<int> v);

int colsum_calc_int(std::vector<std::vector<int>> m, int colid);

} // namespace Utilities
} // namespace EERAModel