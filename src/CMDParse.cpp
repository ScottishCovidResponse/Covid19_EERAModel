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
                                        false, "inference", &allowedModes);
		TCLAP::ValueArg<std::string> dataLocArg("l", "local", "Location of local data repository", 
                                        false, "", "string");
        TCLAP::ValueArg<std::string> outputDirArg("d", "outdir", "Output directory of data files", false, _args.output_dir, "string");
        cmd.add( modeArg );
        cmd.add( dataLocArg );
		cmd.add( structureArg );
        cmd.add( outputDirArg );
		cmd.parse( argc, argv);
        _args.output_dir = outputDirArg.getValue();
		_args.stucture = (Utilities::toUpper(structureArg.getValue()) == "IRISH") ? ModelStructureId::IRISH : ModelStructureId::ORIGINAL;
		_args.mode = (Utilities::toUpper(modeArg.getValue()) == "INFERENCE") ? ModelModeId::INFERENCE : ModelModeId::PREDICTION;
        _args.isLocal = dataLocArg.getValue().empty();
        _args.local_location = (!dataLocArg.getValue().empty()) ? dataLocArg.getValue() : std::string(ROOT_DIR)+"/data";
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}

void EERAModel::ArgumentParser::logArguments(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Arguments]:" << std::endl;
    (*log) << "\t" << "Local: " << ((_args.isLocal) ? "True" : "False") << std::endl;
    if(_args.isLocal){(*log) << "\t" << "Local Directory: " << _args.local_location << std::endl;}
    (*log) << "\t" << "Structure: " << ((_args.stucture == ModelStructureId::IRISH) ? "Irish" : "Original") << std::endl;
    (*log) << "\t" << "Mode: " << ((_args.mode == ModelModeId::INFERENCE) ? "Inference" : "Prediction") << std::endl;
    (*log) << "\t" << "Output Directory: " << _args.output_dir << std::endl;
}
void EERAModel::ArgumentParser::AppendOptions(ModelInputParameters& input_params)
{
    input_params.model_structure = _args.stucture;

    input_params.run_type = _args.mode;
}