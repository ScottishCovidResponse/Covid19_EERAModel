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

IOdatapipeline::IOdatapipeline(string params_path, string dpconfig_path)
{
    ParamsPath = params_path;
    datapipelineActive = false;

    if (dpconfig_path != "")
    {
        const char *uri = "https://whatever"; // I'm guessing this is the repo for the model
        const char *git_sha = "git_sha"; // And this the version ID, need to find these both

        std::cout << "ParamsPath: " << ParamsPath << "\n";
        std::cout << "ConfigPath: " << dpconfig_path << "\n";

        // If something goes wrong in the opening of the data pipeline an exception will get thrown
        dp.reset(new DataPipeline(dpconfig_path, uri, git_sha));
        datapipelineActive = true;
    }
}

CommonModelInputParameters IOdatapipeline::ReadCommonParameters()
{
    CommonModelInputParameters commonParameters;

    if (datapipelineActive)
    {
        commonParameters.paramlist  = ReadFixedModelParameters();
        commonParameters.totN_hcw   = dp->read_estimate("fixed-parameters/total_hcw", "total_hcw");
    }
    else
    {
        commonParameters.paramlist  = IO::ReadFixedModelParameters(ParamsPath);
        commonParameters.totN_hcw   = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", ParamsPath);
    }
    
    commonParameters.herd_id    = ReadNumberFromFile<int>("shb_id", "Settings", ParamsPath);

    return commonParameters;
}

params IOdatapipeline::ReadFixedModelParameters()
{
    std::cout << "(Datapipeline): ReadFixedModelParameters\n";
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