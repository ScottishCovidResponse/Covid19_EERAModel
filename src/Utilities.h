#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <cmath>

namespace EERAModel {
/**
 * @brief Namepace containing all utility functions
 *
 * This namepace contains functions which can be used universally
 * even outside of the model itself
 */
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
 * @brief Check if a directory exists
 * 
 * Checks if the specified address actually exists and so
 * can be written to
 * 
 * @param directory file directory
 * 
 * @return if directory exists
 */
bool directoryExists(const std::string directory);

/**
 * @brief Sum vector elements in blocks
 * 
 * Summate values at every Nth index, e.g. n = 7 would summate daily data to be weekly data. For 
 * cases where n is not a multiple of the vector size, the residual values at the end of the input
 * are dropped.
 * 
 * If n is greater than the size of @p data, or the size of @p data is zero, an empty vector is
 * returned.
 * 
 * @f$f([0,1,2,3,4,5,6,7], 3) -> f([0,1,2,3,4,5], 4) -> [3,12] @f$
 * 
 * @param data	Vector of data to be reduced
 * @param n	Nth value at which summation should occur (i.e. index+1) where @code{.cpp} n <= data.size() @endcode
 * 
 * @return Vector of summation values
 */
template<typename T>
std::vector<T> AccumulateEveryN(const std::vector<T>& data, int n)
{
    std::vector<T> _temp;
    T data_val(0);

    if (n <= data.size())
    {
        for (int i{0}; i < data.size(); ++i)
        {
            data_val += data[i];

            // Append to vector when index + 1 is a multiple of n
            if ((i + 1) % n == 0 )
            {
                _temp.push_back(data_val);
                data_val = T(0);
            }
            
        }
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

/**
 * @brief Read data from a CSV file into vectors
 * 
 * Converts data within a CSV file into a vector of vectors
 * for a particular data type.
 * 
 * @param inputfile address of input file
 * @param delimiter column separator character
 * 
 * @return vector of vectors containing read in data values
 */
template<typename T>
std::vector<std::vector<T>> read_csv(const std::string &inputfile, char delimiter)
{
	std::vector<std::vector<T> > data;
	std::ifstream infile(inputfile.c_str());
	if (infile.fail())  { std::cout << "Input file not found" << std::endl; return data; }
	std::string line;
	std::vector<T> record;
	
	while (std::getline(infile, line,'\n'))
	{
		int linepos = 0;
		int inquotes = false;
		char c;
		int linemax = line.length();
		std::string curstring;
		record.clear();
        while(line[linepos] != 0 && linepos < linemax)
        {
            c = line[linepos];
            
            if (!inquotes && curstring.length() == 0 && c == '"') { 	//beginquotechar
                 inquotes=true;
            }
            else if (inquotes && c=='"')
            {																	//quotechar
                if ( (linepos + 1 < linemax) && (line[linepos+1] == '"') )
                {
                    curstring.push_back(c);										//encountered 2 double quotes in a row (resolves to 1 double quote)
                    linepos++;
                }
                else 
                {
                    inquotes=false;
                }											//endquotechar
            }
            else if (!inquotes && c==delimiter)
            {																	//end of field
                record.push_back(static_cast<T>(atof(curstring.c_str())));
                curstring="";
            }
            else
            {
                curstring.push_back(c);
            }
            linepos++;
        }
        
        record.push_back(static_cast<T>(atof(curstring.c_str())));
        data.push_back(record);    
    }

	return data;
}

/**
 * @brief Convert a string to upper case
 * 
 * Iterates through a string converting the output to be
 * all upper case.
 * 
 * @param str string to be converted
 * 
 * @return upper case string
 */
std::string toUpper(std::string str);

} // namespace Utilities
} // namespace EERAModel
