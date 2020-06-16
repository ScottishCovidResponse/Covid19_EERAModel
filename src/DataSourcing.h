#pragma once

#include <string>

namespace EERAModel
{
namespace DataSourcing
{
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

    DataSource Local(const std::string root_dir=std::string(ROOT_DIR));

    const DataSource Remote;

    enum class SourceID
    {
        REMOTE,
        LOCAL
    };

    const DataSource getSource(SourceID source_id, const std::string local_data_location="");

};
};