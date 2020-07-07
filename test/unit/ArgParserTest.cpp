#include "gtest/gtest.h"
#include "ArgumentParser.h"
#include "ModelTypes.h"

#include <iostream>

using namespace EERAModel;

TEST(AnArgumentParser, RecognisesInferenceMode)
{
    const char* _test_args_inf[] = {"exe", "-m", "inference"};

    ArgumentParser parse_inf(3, _test_args_inf);
    
    EXPECT_EQ(parse_inf.getArgs().mode, ModelModeId::INFERENCE);
}

TEST(AnArgumentParser, RecognisesPredictionMode)
{
    const char* _test_args_prd[] = {"exe", "-m", "prediction"};

    ArgumentParser parse_prd(3, _test_args_prd);
    
    EXPECT_EQ(parse_prd.getArgs().mode, ModelModeId::PREDICTION);
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

TEST(AnArgumentParser, RecognisesNoModelStructure)
{
    const char* _test_args_default[] = {"exe", "-m", "inference"};

    ArgumentParser parse_default(3, _test_args_default);
    
    EXPECT_EQ(parse_default.getArgs().structure, ModelStructureId::UNKNOWN);
}

TEST(AnArgumentParser, RecognisesOutputDirectory)
{
    const char* _test_args[] = {"exe", "-d", "test", "-m", "inference"};

    ArgumentParser parse(5, _test_args);

    EXPECT_EQ(parse.getArgs().output_dir, "test");
}

TEST(AnArgumentParser, SwitchesToDefaultWhenNoOutputDirectoryIsSupplied)
{
    std::string default_dir(std::string(ROOT_DIR)+"/outputs");
    
    const char* _test_args[] = {"exe", "-d", "my_outputs", "-m", "inference"};

    ArgumentParser parse(5, _test_args);

    EXPECT_EQ(parse.getArgs().output_dir, default_dir);
}
