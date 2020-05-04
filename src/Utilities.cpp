#include "Utilities.h"

namespace EERAModel {
namespace Utilities {

double mean_calc_double( std::vector<double> v) {
	double return_value = 0.0;
	int n = v.size();
	for(std::vector<double>::iterator j=v.begin();j!=v.end();++j)  return_value += *j;
    return ( return_value / n);
}

double mean_calc_int( std::vector<int> v) {
	int return_value = 0;
	int n = v.size();
	for(std::vector<int>::iterator j=v.begin();j!=v.end();++j)  return_value += *j;
    return ( (double)return_value / n);
}

void cumsum_calc(std::vector<double> v,std::vector<double>& r_val) {
	int size_vec = v.size();
	double tmp_sum=0.0;
	for (int xx = 0; xx < size_vec; ++xx) {
		tmp_sum += v[xx];
		r_val.push_back(tmp_sum);
	}
}

void cumsum_calc_int(std::vector<int> v,std::vector<int>& r_val) {
	int size_vec = v.size();
	int tmp_sum=0;
	for (int xx = 0; xx < size_vec; ++xx) {
		tmp_sum += v[xx];
		r_val.push_back(tmp_sum);
	}
}

double max_calc(double a,double b) {
	double rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

int max_calc_int(int a,int b) {
	int rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

double sum_calc(std::vector<double> v) {
    double rval = 0.0;
	for(std::vector<double>::iterator j=v.begin();j!=v.end();++j)  rval += *j;
	return rval;
}

int sum_calc_int(std::vector<int> v) {
    int rval = 0;
	for(std::vector<int>::iterator j=v.begin();j!=v.end();++j)  rval += *j;
	return rval;
}

int colsum_calc_int(std::vector<std::vector<int>> m, int colid) {
    int rval = 0;
	for (unsigned int j = 0; j < m.size(); ++j) {
		rval += m[j][colid];
	}
	return rval;
}

double colsum_calc(std::vector<std::vector<double>> m, int colid) {
    double rval = 0.0;
	for (unsigned int j = 0; j < m.size(); ++j) {
		rval += m[j][colid];
	}
	return rval;
}

} // namespace Utilities
} // namespace EERAModel