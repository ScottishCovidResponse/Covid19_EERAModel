#include "CSV.h"

#include <iostream>
#include <fstream>
#include <string>

namespace EERAModel {
namespace CSV {

void read_csv_int(std::vector<std::vector<int> > &data, const std::string &inputfile, char delimiter)
{
	std::ifstream infile(inputfile.c_str());
	if (infile.fail())  { std::cout << "Input file not found" << std::endl; return; }
	std::string line;
	std::vector<int> record;
	
	while (std::getline(infile, line,'\n') )
	{
		int linepos=0;
		int inquotes=false;
		char c;
		int linemax=line.length();
		std::string curstring;
		record.clear();
	while(line[linepos]!=0 && linepos < linemax)
		{
			c = line[linepos];
			
			if (!inquotes && curstring.length()==0 && c=='"'){ inquotes=true;}	//beginquotechar
			else if (inquotes && c=='"')
			{																	//quotechar
				if ( (linepos+1 <linemax) && (line[linepos+1]=='"') )
				{
					curstring.push_back(c);										//encountered 2 double quotes in a row (resolves to 1 double quote)
					linepos++;
				}
				else {inquotes=false;}											//endquotechar
			}
			else if (!inquotes && c==delimiter)
			{																//end of field
				record.push_back(atoi(curstring.c_str()) );
				curstring="";
			}
/*			else if (!inquotes && (c=='\r' || c=='\n') )
			{
				record.push_back( atoi(curstring.c_str()) );
				break;
			}
*/			else
			{
				curstring.push_back(c);
			}
			linepos++;
		}
		record.push_back( atoi(curstring.c_str()) );
		data.push_back(record);
	}
}

void read_csv_double(std::vector<std::vector<double> > &data, const std::string &inputfile, char delimiter)
{
	std::ifstream infile(inputfile.c_str());
	if (infile.fail())  { std::cout << "Input file not found" << std::endl; return; }
	std::string line;
	std::vector<double> record;
	
	while (std::getline(infile, line,'\n') )
	{
		int linepos=0;
		int inquotes=false;
		char c;
		int linemax=line.length();
		std::string curstring;
		record.clear();
	while(line[linepos]!=0 && linepos < linemax)
		{
			c = line[linepos];
			
			if (!inquotes && curstring.length()==0 && c=='"'){ inquotes=true;}	//beginquotechar
			else if (inquotes && c=='"')
			{																	//quotechar
				if ( (linepos+1 <linemax) && (line[linepos+1]=='"') )
				{
					curstring.push_back(c);										//encountered 2 double quotes in a row (resolves to 1 double quote)
					linepos++;
				}
				else {inquotes=false;}											//endquotechar
			}
			else if (!inquotes && c==delimiter)
			{																//end of field
				record.push_back(atof(curstring.c_str()) );
				curstring="";
			}
/*			else if (!inquotes && (c=='\r' || c=='\n') )
			{
				record.push_back( atoi(curstring.c_str()) );
				break;
			}
*/			else
			{
				curstring.push_back(c);
			}
			linepos++;
		}
		record.push_back( atof(curstring.c_str()) );
		data.push_back(record);
	}
}

} // namespace CSV
} // namespace EERAModel