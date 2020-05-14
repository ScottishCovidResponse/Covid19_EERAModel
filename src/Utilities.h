#pragma once

#include <vector>
#include <assert.h>

namespace EERAModel {
namespace Utilities {

template<typename T>
double sse_calc(const std::vector<T>& simval, const std::vector<T>& obsval){
	
	//verify that the 2 vectors have the same size
	assert(simval.size() == obsval.size());

	//compute the sum of squared errors
	double sum_sq=0.0;
	for (unsigned int xx = 0; xx < obsval.size(); ++xx) {
		double error = 0.0;
		error = static_cast<double>( obsval[xx] - simval[xx] ) ;
		sum_sq += std::pow(error,2);
	}

	return sum_sq;
	
}
} // namespace Utilities
} // namespace EERAModel
