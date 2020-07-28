#include "IO-datapipeline.h"
#include "IO.h"
#include "ModelCommon.h"
#include "Utilities.h"
#include "Git.h"

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>

namespace EERAModel {
namespace IO {

namespace {

string dirname(const string& path)
{
    size_t slash = path.rfind('/');
    if (slash == string::npos)
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

string appendFilename(const string &path1, const string &path2)
{
    if (path2.front() == '/')
    {
        return path2;
    }
    else
    {
        return path1 + "/" + path2;
    }
}

}

IOdatapipeline::IOdatapipeline(const string &path)
{
    ParamsPath = path;
    datapipelineActive = false;

    if (ExistsInFile("datapipeline_config", "Settings", ParamsPath))
    {
        /* The config yaml file is assumed to have a relative path compared to the paramaters ini */ 
        string configPath = ReadStringFromFile("datapipeline_config", "Settings", ParamsPath);
        configPath = appendFilename(dirname(ParamsPath), configPath);

        const char *uri = "https://whatever"; // No idea what this is for...
        const char *git_sha = "git_sha"; // Or this...

        std::cout << "ParamsPath: " << ParamsPath << "\n";
        std::cout << "configPath: " << configPath << "\n";

        // Not sure what happens if the next call fails, is there an exception?
        dp.reset(new DataPipeline(configPath, uri, git_sha));
        datapipelineActive = true;
    }
}

CommonModelInputParameters IOdatapipeline::ReadCommonParameters()
{
    CommonModelInputParameters commonParameters;

    if (ExistsInFile("use_datapipeline", "Fixed parameters", ParamsPath) &&
        ReadBoolFromFile("use_datapipeline", "Fixed parameters", ParamsPath))
    {
        if (!datapipelineActive)
        {
            throw std::runtime_error("use_datapipeline requested, but config file not specified");
        }

        commonParameters.paramlist  = ReadFixedModelParameters();
    }
    else
    {
        commonParameters.paramlist  = IO::ReadFixedModelParameters(ParamsPath);
    }
    
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);
    commonParameters.totN_hcw   = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", ParamsPath);

    return commonParameters;
}

params IOdatapipeline::ReadFixedModelParameters()
{
    params paramlist;

    paramlist.T_lat = dp->read_estimate("fixed-parameters/T_lat", "T_lat");
    paramlist.juvp_s = dp->read_estimate("fixed-parameters/juvp_s", "juvp_s");
    paramlist.T_inf = dp->read_estimate("fixed-parameters/T_inf", "T_inf");
    paramlist.T_rec = dp->read_estimate("fixed-parameters/T_rec", "T_rec");
    paramlist.T_sym = dp->read_estimate("fixed-parameters/T_sym", "T_sym");
    paramlist.T_hos = dp->read_estimate("fixed-parameters/T_hos", "T_hos");
    paramlist.K = dp->read_estimate("fixed-parameters/K", "K");
    paramlist.inf_asym = dp->read_estimate("fixed-parameters/inf_asym", "inf_asym");

    return paramlist;
}


} // namespace IO
} // namespace EERAModel