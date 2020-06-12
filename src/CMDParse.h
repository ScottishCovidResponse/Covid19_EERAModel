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
            std::map<std::string, int> _int_args;
            std::map<std::string, std::string> _str_args;
            std::map<std::string, double> _dbl_args;
            std::map<std::string, bool> _bool_args;
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

            /**
             * @brief Returns the value of a given parameter by name
             * 
             * Return the value for a given parameter using the allocated
             * label.
             * 
             * @param name parameter name
             */
            template<typename T>
            T getArg(std::string name)
            {
                if(_str_args.find(name) != _str_args.end())
                {
                    return T(_str_args[name]);
                }

                else if(_bool_args.find(name) != _bool_args.end())
                {
                    return T(_bool_args[name]);
                }

                else if(_dbl_args.find(name) != _dbl_args.end())
                {
                    return T(_dbl_args.end());
                }
                
                else
                {
                    return T(0);
                }
            }
    };
};