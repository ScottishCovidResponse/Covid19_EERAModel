#include "CMDParse.h"

EERAModel::ArgumentParser::ArgumentParser(int argc, char** argv)
{
    try
	{
		TCLAP::CmdLine cmd("A model structure argument and mode of running must be provided",   // Help string message
							' ',                                                                // Delimiter between arguments
							std::string(VERSION),                                               // Software version number
							true                                                                // Whether to generate the help and version flags
						);

		std::vector<std::string> running_modes = {"inference", "prediction"};
        std::vector<std::string> structures = {"original", "irish"};
        TCLAP::ValuesConstraint<std::string> allowedStructures( structures );
		TCLAP::ValuesConstraint<std::string> allowedModes( running_modes );
		TCLAP::ValueArg<std::string> structureArg("s", "structure", 
                                        "Model structure", 
                                        true, "", &allowedStructures);
		TCLAP::ValueArg<std::string> modeArg("m", "mode", "Running mode, if unset use default", 
                                        false, "", &allowedModes);
		TCLAP::ValueArg<std::string> dataLocArg("l", "local", "Location of local data repository", 
                                        false, "", "string");
        TCLAP::ValueArg<std::string> outputDirArg("d", "outdir", "Output directory of data files", false, std::string(ROOT_DIR)+"/outputs", "string");
        cmd.add( modeArg );
        cmd.add( dataLocArg );
		cmd.add( structureArg );
        cmd.add( outputDirArg );
		cmd.parse( argc, argv);
        _str_args["OutputDir"] = outputDirArg.getValue();
		_str_args["Structure"] = structureArg.getValue();
		_str_args["Mode"] = modeArg.getValue();
        _bool_args["Local"] = dataLocArg.getValue().empty();
        _str_args["Data Source"] = (_bool_args["Local"]) ? dataLocArg.getValue() : "scrc";
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}

void EERAModel::ArgumentParser::logArguments(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Arguments]:" << std::endl;
    for(auto& a : _str_args)
    {
        (*log) << "\t" << a.first << ": " << a.second << std::endl;
    }

    for(auto& a : _bool_args)
    {
        (*log) << "\t" << a.first << ": " << ((a.second) ? "True" : "False") << std::endl;
    }

    for(auto& a : _dbl_args)
    {
        (*log) << "\t" << a.first << ": " << a.second << std::endl;
    }
}

void EERAModel::ArgumentParser::AppendOptions(EERAModel::ModelInputParameters& input_params)
{
    if(Utilities::toUpper(_str_args["Structure"]) == "IRISH")
    {
        input_params.model_structure = ModelStructureId::IRISH;
    }

    else
    {
        input_params.model_structure = ModelStructureId::ORIGINAL;
    }

    if(Utilities::toUpper(_str_args["Mode"]) == "INFERENCE")
    {
        input_params.run_type = ModelModeId::INFERENCE;
    }

    else
    {
        input_params.run_type = ModelModeId::PREDICTION;
    }
}