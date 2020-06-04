#pragma once

#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

namespace EERAModel {
namespace Utilities {

/**
 * @brief Compute the sum of squared errors
 * 
 * Compares two vectors of equal length computing the sum of squared errors
 * 
 * @param simval First vector (e.g. simulation values)
 * @param obsval Second vector (e.g. observed values)
 * 
 * @return Sum of squared errors
 */
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

/**
 * @brief Create new vector summating every Nth element
 * 
 * Summate values at every Nth index, e.g. n = 7 would summate
 * daily data to be weekly data
 * 
 * @param data_vector	Vector of data to be reduced
 * @param n	Nth value at which summation should occur (i.e. index+1)
 * 
 * @return Vector of summation values
 */
template<typename T>
std::vector<T> AccumulateEveryNth(const std::vector<T>& data_vector, const int& n)
{
	std::vector<T> _temp = {};
	T data_val(0);

	for(int i{0}; i < data_vector.size(); ++i)
	{
		if(i % n == 0)
		{
			_temp.push_back(data_val);
			data_val = T(0);
		}
		data_val += data_vector[i];
	}

	return _temp;
}

/*! @brief  Logging Stream Class
    @details Class to send output to both cout and output log file
    @date   last modified 2020-05-18
*/
class logging_stream
{
 public:
	using Sptr = std::shared_ptr<logging_stream>;
	
	/*! Create a new logging stream sending the output as a date named file in
	within the specified directory.
	@param out_dir Output directory
	*/
	logging_stream(const std::string& out_dir)
	{
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H-%M-%S",timeinfo);
	_log_time = std::string(buffer);
	std::string _file_name = out_dir+"/logs/run_"+_log_time+".log";

	log_fstream = std::ofstream(_file_name);
	};
	/*! Send outputs specified by << operator to both cout and the log file
	*/
	template<typename T> logging_stream& operator<<(const T& output)
	{
	std::cout << output;
	log_fstream << output;
	return *this;
	}

	// Handle functions such as std::endl etc.
	typedef std::ostream& (*stream_function)(std::ostream&);
	logging_stream& operator<<(stream_function func)
	{
	func(std::cout);
	func(log_fstream);
	return *this;
	}

	std::string getLoggerTime() const {return _log_time;}
	
 private:
	std::ofstream log_fstream;
	std::string _log_time = "";
};
} // namespace Utilities
} // namespace EERAModel
