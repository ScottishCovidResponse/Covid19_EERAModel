#pragma once

namespace EERAModel
{
namespace DataSourcing
{
/**
     * @brief Structure containing file addresses
     * 
     * All file locations for the given data source form
     * be it local or remote are defined in this structure.
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
        std::string prior_params;
    };

};
};