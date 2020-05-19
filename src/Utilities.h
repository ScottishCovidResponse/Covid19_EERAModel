#pragma once

#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <memory>

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
	std::string str(buffer);

	std::string _command = "mkdir -p "+out_dir+"/logs";

	system(_command.c_str());

	std::string _file_name = out_dir+"/logs/run_"+str+".log";

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
 private:
	std::ofstream log_fstream;
};
} // namespace Utilities
} // namespace EERAModel
