#include "Utilities.h"
#include <sys/stat.h>

namespace EERAModel {
namespace Utilities {

std::string toUpper(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), toupper);

    return str;
}

bool directoryExists(const std::string& directory)
{
    struct stat buffer{};

    return (stat(directory.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode));
}

} //namespace Utilities
} //namespace EERAModel