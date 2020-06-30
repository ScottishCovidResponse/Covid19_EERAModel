#pragma once

namespace EERAModel
{
namespace DataSourcing
{
    /**
     * @brief Structure containing local file addresses
     * 
     * All file locations for local data sources are
     * defined in this structure
     */
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
        std::string posterior_params;
    };

};
};