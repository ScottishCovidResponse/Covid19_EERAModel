#include "Utilities.h"


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
};
};