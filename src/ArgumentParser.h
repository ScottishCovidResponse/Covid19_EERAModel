#pragma once

#include <map>

#include <tclap/CmdLine.h>
#include <tclap/Constraint.h>

#include "Utilities.h"
#include "ModelTypes.h"

#ifndef VERSION
#error Macro VERSION must be defined!
#endif


namespace EERAModel
{
struct Arguments
{
    bool isLocal = false;
    std::string local_location=std::string(ROOT_DIR);
    ModelStructureId structure = ModelStructureId::ORIGINAL;
    ModelModeId mode = ModelModeId::INFERENCE;
    std::string output_dir = std::string(ROOT_DIR)+"/outputs";
};

/**
 * @brief Class for parsing model command line arguments
 * 
 * This class uses the TCLAP library to process command line arguments
 * relating to model structure, execution mode and input locations
 * 
 * @param argc number of command line arguments
 * @param argv argument char** strings
 */
class ArgumentParser
{
    private:
        Arguments _args;
    public:
        ArgumentParser(int argc, char** argv);

        /**
         * @brief Log argument values to output
         * 
         * Send printout of all argument values to the specified
         * logging stream
         * 
         * @param logging_stream shared pointer to model logging stream
         */
        void logArguments(Utilities::logging_stream::Sptr log);

        /**
         * @brief Append specified options to input parameters
         * 
         * Specifies options for structure, mode and input location
         * these overwrite values specified within a parameter file.
         * 
         * @param input_params model input parameter object to update
         */
        void AppendOptions(ModelInputParameters& input_params);

        bool runLocal() const {return _args.isLocal;}
        
        std::string localSourceDir() const {return _args.local_location;}

        ModelModeId runMode() const {return _args.mode;}

        ModelStructureId modelStructure() const {return _args.structure;}

        std::string outputDir() const {return _args.output_dir;}

        Arguments getArgs() const {return _args;}
};
};