#include "gtest/gtest.h"
#include "ArgumentParser.h"
#include "ModelTypes.h"

#include <iostream>

TEST(TestArgumentParser, TestModeSwitching)
{
    char* _test_args_inf[] = {"exe", "-m", "inference", NULL};
    int argc_inf = sizeof(_test_args_inf)/sizeof(char*) - 1;
    char** _p_test_args_inf = _test_args_inf;
    char* _test_args_prd[] = {"exe", "-m", "prediction", NULL};
    int argc_prd = sizeof(_test_args_prd)/sizeof(char*) - 1;
    char** _p_test_args_prd = _test_args_prd;

    EERAModel::ArgumentParser parse_inf(argc_inf, _p_test_args_inf);
    EERAModel::ArgumentParser parse_prd(argc_prd, _p_test_args_prd);
    
    EXPECT_EQ(parse_inf.getArgs().mode, EERAModel::ModelModeId::INFERENCE);
    EXPECT_EQ(parse_prd.getArgs().mode, EERAModel::ModelModeId::PREDICTION);
}

TEST(TestArgumentParser, TestStructureSwitching)
{
    char* _test_args_inf[] = {"exe", "-s", "original", "-m", "inference", NULL};
    int argc_inf = sizeof(_test_args_inf)/sizeof(char*) - 1;
    char** _p_test_args_inf = _test_args_inf;
    char* _test_args_prd[] = {"exe", "-s", "irish", "-m", "inference", NULL};
    int argc_prd = sizeof(_test_args_prd)/sizeof(char*) - 1;
    char** _p_test_args_prd = _test_args_prd;
    char* _test_args_default[] = {"exe", "-m", "inference", NULL};
    int argc_default = sizeof(_test_args_default)/sizeof(char*) - 1;
    char** _p_test_default = _test_args_default;

    EERAModel::ArgumentParser parse_inf(argc_inf, _p_test_args_inf);
    EERAModel::ArgumentParser parse_prd(argc_prd, _p_test_args_prd);
    EERAModel::ArgumentParser parse_default(argc_default, _p_test_default);
    
    EXPECT_EQ(parse_inf.getArgs().structure, EERAModel::ModelStructureId::ORIGINAL);
    EXPECT_EQ(parse_prd.getArgs().structure, EERAModel::ModelStructureId::IRISH);
    EXPECT_EQ(parse_default.getArgs().structure, EERAModel::ModelStructureId(0));
}

TEST(TestArgumentParser, TestExistingOutputDirectory)
{
    char* _test_args[] = {"exe", "-d", "test", "-m", "inference", NULL};
    int argc = sizeof(_test_args)/sizeof(char*) - 1;
    char** _p_test_args = _test_args;

    EERAModel::ArgumentParser parse(argc, _p_test_args);

    EXPECT_EQ(parse.getArgs().output_dir, "test");
}

TEST(TestArgumentParser, TestNonExistingOutputDirectory)
{
    char* _test_args[] = {"exe", "-d", "my_outputs", "-m", "inference", NULL};
    int argc = sizeof(_test_args)/sizeof(char*) - 1;
    char** _p_test_args = _test_args;

    EERAModel::ArgumentParser parse(argc, _p_test_args);

    EXPECT_EQ(parse.getArgs().output_dir, std::string(ROOT_DIR)+"/outputs");
}