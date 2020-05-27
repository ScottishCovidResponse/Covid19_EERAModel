#include <iostream>
#include "Restful.h"

int main()
{
    std::string models_url = "http://data.scrc.uk/api/model/";
    auto j = Restful::get_json(models_url);
    std::cout << j.dump(4) << std::endl;
}