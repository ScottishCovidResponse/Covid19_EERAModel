#include "Utilities.h"
#include <sys/stat.h>
#include <ctime>

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

bool fileExists(const std::string& file)
{
    struct stat buffer;

    return (stat(file.c_str(), &buffer) == 0 && S_ISREG(buffer.st_mode));
}

std::string dirname(std::string path)
{
    size_t slash = path.rfind('/');
    if (slash == std::string::npos)
    {
        return ".";
    }
    else if (slash == 0)
    {
        return "/";
    }
    else
    {
        return path.substr(0, slash);
    }
}

std::string appendPath(std::string path1, std::string path2)
{
    if (path2.front() == '/')
    {
        return path2;
    }
    else if (path1 == "/" || path1 == "")
    {
        return path1 + path2;
    }
    else
    {
        return path1 + "/" + path2;
    }
}


} //namespace Utilities
} //namespace EERAModel
