#pragma once

#include <string>

namespace EERAModel
{
namespace DataSourcing
{
    enum class SourceID
    {
        LOCAL,
        REMOTE
    };

    struct DataFiles
    {
        std::string parameters;
        std::string data;
        std::string deaths;
        std::string ages;

        struct WAIFW
        {
            std::string norm;
            std::string home;
            std::string sdist;
        };
        WAIFW waifw;
        std::string cfr_byage;
        std::string frail;
        std::string prior_params;
    };

    class DataSource
    {
        public:
            const DataFiles data_files;
    };

    class Local : DataSource
    {
        public:
            Local() : DataSource({
                std::string(ROOT_DIR)+"/data/parameters.ini",
                std::string(ROOT_DIR)+"/data/scot_data.csv",
                std::string(ROOT_DIR)+"/data/scot_deaths.csv",
                std::string(ROOT_DIR)+"/data/scot_age.csv",
                {std::string(ROOT_DIR)+"/data/waifw_norm.csv",
                 std::string(ROOT_DIR)+"/data/waifw_home.csv",
                 std::string(ROOT_DIR)+"/data/waifw_sdist.csv"},
                std::string(ROOT_DIR)+"/data/cfr_byage.csv",
                std::string(ROOT_DIR)+"/data/scot_frail.csv",
                std::string(ROOT_DIR)+"/src/prior_particle_params.csv"
            }
            )
    };

    class Remote : DataSource
    {

    };

    DataSource* fetchData(SourceID source_id);
};
};