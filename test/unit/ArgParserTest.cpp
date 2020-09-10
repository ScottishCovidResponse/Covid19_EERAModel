#include "gtest/gtest.h"
#include "ArgumentParser.h"
#include "ModelTypes.h"

#include <iostream>

using namespace EERAModel;
using namespace testing;

TEST(AnArgumentParser, RecognisesInferenceMode)
{
    const char* _test_args_inf[] = {"exe", "-m", "inference", "-s", "original"};

    ArgumentParser parse_inf(5, _test_args_inf);
    
    EXPECT_EQ(parse_inf.getArgs().mode, ModelModeId::INFERENCE);
}

TEST(AnArgumentParser, RecognisesPredictionMode)
{
    const char* _test_args_prd[] = {"exe", "-m", "prediction", "-s", "original"};

    ArgumentParser parse_prd(5, _test_args_prd);
    
    EXPECT_EQ(parse_prd.getArgs().mode, ModelModeId::PREDICTION);
}

TEST(AnArgumentParser, TerminatesIfNoModeProvided)
{
    const char* _test_args_default[] = {"exe", "-s", "original"};

    EXPECT_EXIT( ArgumentParser(3, _test_args_default), ExitedWithCode(1), "");
}

TEST(AnArgumentParser, RecognisesOriginalModelStructure)
{
    const char* _test_args_inf[] = {"exe", "-s", "original", "-m", "inference"};

    ArgumentParser parse_inf(5, _test_args_inf);
    
    EXPECT_EQ(parse_inf.getArgs().structure, ModelStructureId::ORIGINAL);
}

TEST(AnArgumentParser, RecognisesIrishModelStructure)
{
    const char* _test_args_prd[] = {"exe", "-s", "irish", "-m", "inference"};

    ArgumentParser parse_prd(5, _test_args_prd);
    
    EXPECT_EQ(parse_prd.getArgs().structure, ModelStructureId::IRISH);
}

TEST(AnArgumentParser, RecognisesIrish2ModelStructure)
{
    const char* _test_args_prd[] = {"exe", "-s", "irish2", "-m", "inference"};

    ArgumentParser parse_prd(5, _test_args_prd);
    
    EXPECT_EQ(parse_prd.getArgs().structure, ModelStructureId::TEMP);
}

TEST(AnArgumentParser, TerminatesIfNoModelStructureProvided)
{
    const char* _test_args_default[] = {"exe", "-m", "inference"};

    EXPECT_EXIT( ArgumentParser(3, _test_args_default), ExitedWithCode(1), "");
}

TEST(AnArgumentParser, RecognisesOutputDirectory)
{
    const char* _test_args[] = {"exe", "-d", "test", "-m", "inference", "-s", "original"};

    ArgumentParser parse(7, _test_args);

    EXPECT_EQ(parse.getArgs().output_dir, "test");
}

TEST(AnArgumentParser, SwitchesToDefaultWhenNoOutputDirectoryIsSupplied)
{
    std::string default_dir(std::string(ROOT_DIR)+"/outputs");
    
    const char* _test_args[] = {"exe", "-d", "my_outputs", "-m", "inference", "-s", "original"};

    ArgumentParser parse(7, _test_args);

    EXPECT_EQ(parse.getArgs().output_dir, default_dir);
}

TEST(AnArgumentParser, DefaultsToZeroParameterSetIndex)
{    
    const char* _test_args[] = {"exe", "-m", "inference", "-s", "original"};

    ArgumentParser parse(5, _test_args);

    EXPECT_EQ(parse.getArgs().parameter_set_index, 0);
}

TEST(AnArgumentParser, RecognisesParameterSetIndex)
{    
    const char* _test_args[] = {"exe", "-m", "inference", "-s", "original", "-i", "2"};

    ArgumentParser parse(7, _test_args);

    EXPECT_EQ(parse.getArgs().parameter_set_index, 2);
}

TEST(AnArgumentParser, CheckAbsenseOfDataPipeline)
{    
    const char* _test_args[] = {"exe", "-m", "inference", "-s", "original"};

    ArgumentParser parse(5, _test_args);

    EXPECT_EQ(parse.getArgs().datapipeline_path, "");
}

TEST(AnArgumentParser, RecognisesRequestForDataPipeline)
{    
    const char* _test_args[] = {"exe", "-m", "inference", "-s", "original", "-c", "my_config.yaml"};

    ArgumentParser parse(7, _test_args);

    EXPECT_EQ(parse.getArgs().datapipeline_path, "my_config.yaml");
}

TEST(AnArgumentParser, TerminatesIfNoDataPipelinePath)
{
    const char* _test_args_default[] = {"exe", "-m", "inference", "-s", "original", "-c"};

    EXPECT_EXIT( ArgumentParser(6, _test_args_default), ExitedWithCode(1), "");
}
