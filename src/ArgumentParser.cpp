#include "ArgumentParser.h"

namespace EERAModel
{
ArgumentParser::ArgumentParser(int argc, char** argv)
{
    try
    {
        TCLAP::CmdLine cmd("Run model with specified set of options",                           // Help string message
                            ' ',                                                                // Delimiter between arguments
                            std::string(VERSION),                                               // Software version number
                            true                                                                // Whether to generate the help and version flags
                        );

        std::vector<std::string> running_modes = {"inference", "prediction"};
        std::vector<std::string> structures = {"original", "irish"};
        TCLAP::ValuesConstraint<std::string> allowedStructures( structures );
        TCLAP::ValuesConstraint<std::string> allowedModes( running_modes );
        TCLAP::ValueArg<std::string> structureArg("s", "structure", 
                                        "Model structure, if unset use 'original'", 
                                        false, "", &allowedStructures);
        TCLAP::ValueArg<std::string> modeArg("m", "mode", "Running mode, if unset use parameter file value else 'inference'", 
                                        false, "", &allowedModes);
        TCLAP::ValueArg<std::string> dataLocArg("l", "local", "Location of local data repository", 
                                        false, "", "string");
        TCLAP::ValueArg<std::string> outputDirArg("d", "outdir", "Output directory of data files", false, _args.output_dir, "string");
        cmd.add( modeArg );
        cmd.add( dataLocArg );
        cmd.add( structureArg );
        cmd.add( outputDirArg );
        cmd.parse( argc, argv);

        _args.output_dir = (Utilities::directoryExists(outputDirArg.getValue())) ? outputDirArg.getValue() : _args.output_dir;

        if(Utilities::toUpper(structureArg.getValue()) == "IRISH")
        {
            _args.structure = ModelStructureId::IRISH;
        }
        else if(Utilities::toUpper(structureArg.getValue()) == "ORIGINAL")
        {
            _args.structure = ModelStructureId::ORIGINAL;
        }
        else 
        {
            _args.structure = ModelStructureId::UNKNOWN;
        }

        if(Utilities::toUpper(modeArg.getValue()) == "INFERENCE")
        {
            _args.mode = ModelModeId::INFERENCE;
        }
        else if(Utilities::toUpper(modeArg.getValue()) == "PREDICTION")
        {
            _args.mode = ModelModeId::PREDICTION;
        }
        else
        {
            _args.mode = ModelModeId::UNKNOWN;
        }

        _args.isLocal = dataLocArg.getValue().empty();
        _args.local_location = (!dataLocArg.getValue().empty()) ? dataLocArg.getValue() : std::string(ROOT_DIR);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

void ArgumentParser::logArguments(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Arguments]:" << std::endl;
    (*log) << "\t" << "Local: " << ((_args.isLocal) ? "True" : "False") << std::endl;
    if(_args.isLocal){(*log) << "\t" << "Local Directory: " << _args.local_location << std::endl;}
    (*log) << "\t" << "Structure: ";
    if(_args.structure != ModelStructureId(0))
    {
        (*log) << ((_args.structure == ModelStructureId::IRISH) ? "Irish" : "Original") << std::endl;
    }
    else
    {
        (*log) << "Default" << std::endl;
    }
    if(_args.mode != ModelModeId(0))
    {
        (*log) << "\t" << "Mode: " << ((_args.mode == ModelModeId::INFERENCE) ? "Inference" : "Prediction") << std::endl;
    }
    else
    {
        (*log) << "Default" << std::endl;
    }
    (*log) << "\t" << "Output Directory: " << _args.output_dir << std::endl;
}

void ArgumentParser::AppendOptions(ModelInputParameters& input_params)
{
    // Only set model structure to default of "original" if no option specified
    // in parameters file
    if(input_params.model_structure == ModelStructureId::UNKNOWN)
    {
        input_params.model_structure = ModelStructureId::ORIGINAL;
    }

    // Only overwrite this model structure option if command line argument is not empty
    if(_args.structure != ModelStructureId::UNKNOWN)
    {
        input_params.model_structure = _args.structure;
    }

    // Only set run mode to default of "inference" if no option specified
    // in parameters file
    if(input_params.run_type == ModelModeId::UNKNOWN)
    {
        input_params.run_type = ModelModeId::INFERENCE;
    }

    // Only overwrite this run mode option if command line argument is not empty
    if(_args.mode != ModelModeId::UNKNOWN)
    {
        input_params.run_type = _args.mode;
    }
}
};