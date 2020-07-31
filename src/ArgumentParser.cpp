#include "ArgumentParser.h"

namespace EERAModel {

ArgumentParser::ArgumentParser(int argc, const char* const * argv)
{
    for (int i = 0; i < argc; ++i)
        raw_args_.push_back(std::string(argv[i]));
        
    try
    {
        TCLAP::CmdLine cmd("Run model with specified set of options",   // Help string message
                            ' ',                                        // Delimiter between arguments
                            std::string(VERSION),                       // Software version number
                            true                                        // Whether to generate the help and version flags
                        );

        std::vector<std::string> running_modes = {"inference", "prediction"};
        std::vector<std::string> structures = {"original", "irish", "irish2"};
        TCLAP::ValuesConstraint<std::string> allowedStructures(structures);
        TCLAP::ValuesConstraint<std::string> allowedModes(running_modes);
        TCLAP::ValueArg<std::string> structureArg("s", "structure", 
                                        "Model structure. Can be original, irish or irish2", 
                                        true, "", &allowedStructures);
        TCLAP::ValueArg<int> indexArg("i", "index", 
                                        "Parameter set index for forward prediction mode. Defaults to zero.", 
                                        false, 0, "integer");
        TCLAP::ValueArg<std::string> modeArg("m", "mode", "Running mode. Can be either inference or prediction.", 
                                        true, "", &allowedModes);
        TCLAP::ValueArg<std::string> dataLocArg("l", "local", "Location of local data repository", 
                                        false, "", "string");
        TCLAP::ValueArg<std::string> dpConfigArg("c", "dpconfig", "yaml file used for the data download", false, _args.datapipeline_path, "string");
        TCLAP::ValueArg<std::string> outputDirArg("d", "outdir", "Output directory of data files", false, _args.output_dir, "string");
        cmd.add( modeArg );
        cmd.add( dataLocArg );
        cmd.add( structureArg );
        cmd.add( indexArg );
        cmd.add( dpConfigArg );
        cmd.add( outputDirArg );
        cmd.parse( argc, argv);

        _args.output_dir = (Utilities::directoryExists(outputDirArg.getValue())) ? outputDirArg.getValue() : _args.output_dir;

        if (Utilities::toUpper(structureArg.getValue()) == "IRISH")
            _args.structure = ModelStructureId::IRISH;
        else if (Utilities::toUpper(structureArg.getValue()) == "IRISH2")
            _args.structure = ModelStructureId::TEMP;
        else 
            _args.structure = ModelStructureId::ORIGINAL;

        if (modeArg.getValue() == "inference")
            _args.mode = ModelModeId::INFERENCE;
        else
            _args.mode = ModelModeId::PREDICTION;

        if (! dpConfigArg.getValue().empty())
            _args.datapipeline_path = dpConfigArg.getValue();

        _args.isLocal = dataLocArg.getValue().empty();
        _args.local_location = (!dataLocArg.getValue().empty()) ? dataLocArg.getValue() : std::string(ROOT_DIR);

        _args.parameter_set_index = indexArg.getValue();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        exit(1);
    }
}

void ArgumentParser::logArguments(Utilities::logging_stream::Sptr log)
{
    (*log) << "[Arguments]:" << std::endl;
    
    (*log) << "\t" << "Command line: ";
    for (const auto& s : raw_args_) (*log) << s << " ";
    (*log) << std::endl;

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
    (*log) << "\t" << "Data pipeline config file: " << _args.datapipeline_path << std::endl;
    (*log) << "\t" << "Output Directory: " << _args.output_dir << std::endl;
    (*log) << "\t" << "Forward prediction parameter set index: " << _args.parameter_set_index << std::endl;
}

void ArgumentParser::AppendOptions(SupplementaryInputParameters& input_params)
{
    input_params.model_structure = _args.structure;
    input_params.run_type = _args.mode;
}

} // namespace EERAModel
