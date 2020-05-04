#include "DistanceComputation.h"
#include "Utilities.h"

#include <iostream>
#include <cmath>
#include <numeric>

namespace EERAModel {
namespace DistanceComputation {

double sse_calc(int Npop, std::vector<double> simHosp, const std::vector<int>& obsHosp){
	double sum_sq=0.0;
	//verify that the 2 vectors have the same size
	if(simHosp.size() != obsHosp.size()){
		std::cout<< "Warning!!! Simulated and observed vectors of different sizes." << std::endl;
	}
	//compute the sum of squared errors
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = simHosp[xx]/Npop - obsHosp[xx]/Npop;
		sum_sq += pow(error,2);
	}
	return sum_sq;
}

double sse_calc_int(std::vector<int> simval, std::vector<int> obsval){
	
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

double nsse_calc_int(std::vector<int> simHosp, const std::vector<int>& obsHosp){
	
	//verify that the 2 vectors have the same size
	if(simHosp.size() != obsHosp.size()){
		std::cout<< "Warning!!! Simulated and observed vectors of different sizes." << std::endl;
	}
	
	//maximum numbers of cases
	int max_obs = std::accumulate(obsHosp.begin(), obsHosp.end(), 0);
	std::cout << "max_obs: " << max_obs <<'\n';
	//compute the sum of squared errors
	double sum_sq=0.0;
	//compute the sum of squared errors
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = (double)simHosp[xx] - (double)obsHosp[xx];
		sum_sq += pow(error,2);
	}
	
	//comute the normalised sum of squared error by the total number of cases
	double nsse = sum_sq / (double)max_obs;

	return nsse;
	
}

double sst_calc(int Npop, std::vector<double> simHosp, const std::vector<int>& obsHosp){
	double sum_tot=0.0;
	double mean_obs=0.0;
	//verify that the 2 vectors have the same size

	if(simHosp.size() != obsHosp.size()){
		std::cout << obsHosp.size() << "vs. " << simHosp.size() << std::endl;
		std::cout<< "Warning!!! Simulated and observed vectors of different sizes." << std::endl;
	}

	//compute the mean of the observed value
	mean_obs = ::EERAModel::Utilities::mean_calc_int(obsHosp);
	//compute the sum of squared total
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = simHosp[xx]/Npop - mean_obs/Npop;
		sum_tot += pow(error,2);
	}
	return sum_tot;
}

double sst_calc_int(int Npop, std::vector<int> simHosp, const std::vector<int>& obsHosp){
	double sum_tot=0.0;
	double mean_obs=0.0;
	//verify that the 2 vectors have the same size

	if(simHosp.size() != obsHosp.size()){
		std::cout << obsHosp.size() << "vs. " << simHosp.size() << std::endl;
		std::cout<< "Warning!!! Simulated and observed vectors of different sizes." << std::endl;
	}

	//compute the mean of the observed value
	mean_obs = ::EERAModel::Utilities::mean_calc_int(obsHosp);
	//compute the sum of squared total
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = ((double)simHosp[xx]/Npop) - (mean_obs/Npop);
		sum_tot += pow(error,2);
	}
	return sum_tot;
}


} // namespace DistanceComputation
} // namespace EERAModel