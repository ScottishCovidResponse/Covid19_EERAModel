#pragma once

#include <vector>
#include <assert.h>

namespace EERAModel {
namespace Utilities {

template<typename T>
double sse_calc(std::vector<T> simval, const std::vector<int>& obsval){
	
	//verify that the 2 vectors have the same size
	assert(simval.size() == obsval.size());

	//compute the sum of squared errors
	double sum_sq=0.0;
	for (unsigned int xx = 0; xx < obsval.size(); ++xx) {
		double error = 0.0;
		error = (double)( obsval[xx] - simval[xx] ) ;
		sum_sq += std::pow(error,2);
	}

	return sum_sq;
	
}


template<typename T>
T max_calc(T a, T b) {
	T rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

} // namespace Utilities
} // namespace EERAModel
