#pragma once

#include <vector>
#include <string>

namespace EERAModel {
namespace CSV {

void read_csv_int(std::vector<std::vector<int> >& data, const std::string& inputfile, char delimiter);

void read_csv_double(std::vector<std::vector<double> > &data, const std::string &inputfile, char delimiter);

} // namespace CSV
} // namespace EERAModel