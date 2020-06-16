#include "DataSourcing.h"

namespace EERAModel
{
namespace DataSourcing
{
DataSource Local(const std::string root_dir)
{

    return DataSource({
        root_dir+"/data/parameters.ini",
        root_dir+"/data/scot_data.csv",
        root_dir+"/data/scot_deaths.csv",
        root_dir+"/data/scot_age.csv",
        {root_dir+"/data/waifw_norm.csv",
            root_dir+"/data/waifw_home.csv",
            root_dir+"/data/waifw_sdist.csv"},
        root_dir+"/data/cfr_byage.csv",
        root_dir+"/data/scot_frail.csv",
        root_dir+"/src/prior_particle_params.csv"});

}

const DataSource getSource(SourceID source_id, const std::string local_data_location)
{
    return (source_id == SourceID::LOCAL) ? Local(local_data_location) : Remote;
}
};
};