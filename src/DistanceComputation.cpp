#include "DistanceComputation.h"
#include "Utilities.h"

#include <iostream>
#include <cmath>
#include <numeric>

namespace EERAModel {
namespace DistanceComputation {


double sse_calc_int(std::vector<int> simval, const std::vector<int>& obsval){
	
	//verify that the 2 vectors have the same size
	if(simval.size() != obsval.size()){
		std::cout << "Warning!!! Simulated and observed vectors of different sizes." << std::endl;
	}

	//compute the sum of squared errors
	double sum_sq=0.0;
	for (unsigned int xx = 0; xx < obsval.size(); ++xx) {
		double error = 0.0;
		error = (double)( obsval[xx] - simval[xx] ) ;
		sum_sq += std::pow(error,2);
	}

	return sum_sq;
	
}
} // namespace DistanceComputation
} // namespace EERAModel
