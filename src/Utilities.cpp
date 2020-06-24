#include "Utilities.h"
#include <sys/stat.h>

namespace EERAModel
{
namespace Utilities
{
    std::string toUpper(std::string str)
    {
        std::string _temp = str;
        std::transform(_temp.begin(), _temp.end(), _temp.begin(), toupper);

        return _temp;
    }

    bool directoryExists(const std::string directory)
    {
        struct stat buffer;
        
        return (stat(directory.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode));
    }
};
};