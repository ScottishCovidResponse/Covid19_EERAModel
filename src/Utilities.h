#pragma once

#include <vector>

namespace EERAModel {
namespace Utilities {

template<typename T>
double mean_calc( std::vector<T> v) {
	T return_value = 0.0;
	int n = v.size();
	for(std::vector<double>::iterator j=v.begin();j!=v.end();++j)  return_value += *j;
    return ( return_value / n);
}

template<typename T>
void cumsum_calc(std::vector<T> v,std::vector<T>& r_val) {
	int size_vec = v.size();
	T tmp_sum=0.0;
	for (int xx = 0; xx < size_vec; ++xx) {
		tmp_sum += v[xx];
		r_val.push_back(tmp_sum);
	}
}

template<typename T>
T max_calc(T a, T b) {
	T rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

template<typename T>
T sum_calc(std::vector<T> v) {
    T rval = 0.0;
	for(std::vector<double>::iterator j=v.begin();j!=v.end();++j)  rval += *j;
	return rval;
}

template<typename T>
T colsum_calc(std::vector<std::vector<T>> m, int colid) {
    T rval = 0;
	for (unsigned int j = 0; j < m.size(); ++j) {
		rval += m[j][colid];
	}
	return rval;
}


} // namespace Utilities
} // namespace EERAModel
