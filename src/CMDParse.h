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
    class ArgumentParser
    {
        private:
            std::map<std::string, int> _int_args;
            std::map<std::string, std::string> _str_args;
            std::map<std::string, double> _dbl_args;
            std::map<std::string, bool> _bool_args;
        public:
            ArgumentParser(int argc, char** argv);
            void logArguments(Utilities::logging_stream::Sptr log);
            void AppendOptions(ModelInputParameters& input_params);

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