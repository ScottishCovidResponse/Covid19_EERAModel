#include "gtest/gtest.h"
#include "ModelTypes.h"
#include "IO-datapipeline.h"
#include "IO.h"
#include "Git.h"

#include <iostream>

using namespace EERAModel;
using namespace testing;

namespace {
    /* This is needed to start the python interpreter for the C++ data pipeline API.
     * It seems to fail if called once per test, but here as a single instance seems fine. */
    pybind11::scoped_interpreter guard{};
}

TEST(TestIODatapipeline, ExpectThrowForBadPath)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);
    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    EXPECT_ANY_THROW(IO::IOdatapipeline idp2(paramsFile, "", "", logger, "NoValidPath.yaml"););
}

TEST(TestIODatapipeline, CanReadFixedParameters)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    IO::IOdatapipeline idp(paramsFile, "", "", logger, datapipelineConfig);
    CommonModelInputParameters params = idp.ReadCommonParameters();

    EXPECT_EQ(params.paramlist.T_lat, 4);
    EXPECT_EQ(params.paramlist.juvp_s, 0.1);
    EXPECT_EQ(params.paramlist.T_inf, 1.5);
    EXPECT_EQ(params.paramlist.T_rec, 11);
    EXPECT_EQ(params.paramlist.T_sym, 7.0);
    EXPECT_EQ(params.paramlist.T_hos, 5);
    EXPECT_EQ(params.paramlist.K, 2000);
    EXPECT_EQ(params.paramlist.inf_asym, 1.0);
    EXPECT_EQ(params.totN_hcw, 112974);
    EXPECT_EQ(params.day_shut, 19);
}

namespace {
    void compare_eq(int a, int b)
    { 
        EXPECT_EQ(a, b);
    }

    void compare_eq(double a, double b)
    { 
        EXPECT_DOUBLE_EQ(a, b);
    }

    template <class T>
    void compare_eq(
        const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b,
        int firstj = -1, int endi = -1)
    {
        EXPECT_EQ(a.size(), b.size());
        std::size_t sz1 = std::min(a.size(), b.size());

        if (firstj < 0) firstj = 0;

        for (int j = firstj; j < sz1; ++j) {
            std::size_t sz2 = std::min(a[j].size(), b[j].size());

            if (endi < 0) {
                EXPECT_EQ(a[j].size(), b[j].size());
            } else {
                EXPECT_EQ(sz2, endi);
            }

            for (int i = 0; i < sz2; ++i) {
                compare_eq(a[j][i], b[j][i]);
            }
        }
    }
}

TEST(TestIODatapipeline, CanReadModelData)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string configDir = std::string(ROOT_DIR)+"/test/regression/run1/data";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, "", logger, datapipelineConfig);
    ObservationsForModels dp_params = idp.ReadModelObservations();

    // Load data from regression test 1
    ObservationsForModels rg_params = IO::ReadModelObservations(configDir, logger);

    compare_eq<int>(dp_params.cases, rg_params.cases, 1);
    compare_eq<double>(dp_params.age_pop, rg_params.age_pop);
    compare_eq<double>(dp_params.waifw_norm, rg_params.waifw_norm);
    compare_eq<double>(dp_params.waifw_home, rg_params.waifw_home);
    compare_eq<double>(dp_params.waifw_sdist, rg_params.waifw_sdist);
    compare_eq<double>(dp_params.cfr_byage, rg_params.cfr_byage, -1, 3);
}

TEST(TestIODatapipeline, CanReadInferenceConfig)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string configDir = std::string(ROOT_DIR)+"/test/regression/run1/data";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, "", logger, datapipelineConfig);
    CommonModelInputParameters common_params = idp.ReadCommonParameters();
    ObservationsForModels model_obs = idp.ReadModelObservations();
    InferenceConfig dp_infconfig = idp.ReadInferenceConfig(common_params, model_obs);

    EXPECT_EQ(dp_infconfig.prior_pinf_shape1, 3.0);
    EXPECT_EQ(dp_infconfig.prior_pinf_shape2, 9.0);
    EXPECT_EQ(dp_infconfig.prior_phcw_shape1, 3.0);
    EXPECT_EQ(dp_infconfig.prior_phcw_shape2, 3.0);
    EXPECT_EQ(dp_infconfig.prior_chcw_mean, 42);
    EXPECT_EQ(dp_infconfig.prior_d_shape1, 3.0);
    EXPECT_EQ(dp_infconfig.prior_d_shape2, 3.0);
    EXPECT_EQ(dp_infconfig.prior_q_shape1, 3.0);
    EXPECT_EQ(dp_infconfig.prior_q_shape2, 3.0);
    EXPECT_EQ(dp_infconfig.prior_ps_shape1, 9.0);
    EXPECT_EQ(dp_infconfig.prior_ps_shape2, 3.0);
    EXPECT_EQ(dp_infconfig.prior_rrd_shape1, 1.0);
    EXPECT_EQ(dp_infconfig.prior_rrd_shape2, 1.0);
    EXPECT_EQ(dp_infconfig.prior_lambda_shape1, 1e-9);
    EXPECT_EQ(dp_infconfig.prior_lambda_shape2, 1e-6);

    // Load from local data
    InferenceConfig rg_infconfig = IO::ReadInferenceConfig(configDir, logger, common_params);

    compare_eq<int>(dp_infconfig.observations.cases, rg_infconfig.observations.cases, 1);
    compare_eq<int>(dp_infconfig.observations.deaths, rg_infconfig.observations.deaths, 1);
}

TEST(TestIODatapipeline, CanReadPredictionConfig)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string configDir = std::string(ROOT_DIR)+"/test/regression/run1/data";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, "", logger, datapipelineConfig);
    CommonModelInputParameters common_params = idp.ReadCommonParameters();
    PredictionConfig dp_predconfig = idp.ReadPredictionConfig(0, common_params);

    EXPECT_EQ(dp_predconfig.posterior_parameters[0], 0.153532);
    EXPECT_EQ(dp_predconfig.posterior_parameters[1], 0.60916);
    EXPECT_EQ(dp_predconfig.posterior_parameters[2], 37.9059);
    EXPECT_EQ(dp_predconfig.posterior_parameters[3], 0.525139);
    EXPECT_EQ(dp_predconfig.posterior_parameters[4], 0.313957);
    EXPECT_EQ(dp_predconfig.posterior_parameters[5], 0.787278);
    EXPECT_EQ(dp_predconfig.posterior_parameters[6], 0.516736);
    EXPECT_EQ(dp_predconfig.posterior_parameters[7], 8.50135E-07);
    EXPECT_EQ(dp_predconfig.fixedParameters.T_lat, 4);
    EXPECT_EQ(dp_predconfig.fixedParameters.juvp_s, 0.1);
    EXPECT_EQ(dp_predconfig.fixedParameters.T_inf, 1.5);
    EXPECT_EQ(dp_predconfig.fixedParameters.T_rec, 11);
    EXPECT_EQ(dp_predconfig.fixedParameters.T_sym, 7);
    EXPECT_EQ(dp_predconfig.fixedParameters.T_hos, 5);
    EXPECT_EQ(dp_predconfig.fixedParameters.K, 2000);
    EXPECT_EQ(dp_predconfig.fixedParameters.inf_asym, 1);
}

TEST(TestIODatapipeline, CanWriteInferenceData)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string configDir = std::string(ROOT_DIR)+"/test/regression/run1/data";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    std::string timeStamp = logger->getLoggerTime();
    
    {
        IO::IOdatapipeline idp(paramsFile, configDir, "", logger, datapipelineConfig, timeStamp);

        std::vector<particle> particleList = {
            particle{ .nsse_cases=1,
                .nsse_deaths=2,
                .parameter_set={ 3, 4, 5, 6 },
                .iter = 0,
                .weight = 7,
                .simu_outs = { 8, 9 },
                .hospital_death_outs = { 10, 11 },
                .death_outs = { 12, 13 },
                .end_comps = {
                    Compartments{
                        .S = 14, .E = 15, .E_t = 16, .I_p = 17, .I_t = 18,
                        .I1 = 19, .I2 = 20, .I3 = 21, .I4 = 22,
                        .I_s1 = 23, .I_s2 = 24, .I_s3 = 25, .I_s4 = 26,
                        .H = 27, .R = 28, .D = 29
                    }
                }
            }
        };

        EXPECT_NO_THROW(idp.WriteOutputsToFiles(0, 0, 1, 0, particleList, "unit_test"));
    }

    // Can't reread the data back currently via the API without upload and downloading

    // {
    //     std::string uri = GitMetadata::URL();
    //     std::string git_sha = GitMetadata::CommitSHA1();

    //     DataPipeline dp(datapipelineConfig, uri.c_str(), git_sha.c_str());

    //     Table simu_table = dp.read_table("outputs/unit_test/inference/" + timeStamp, "steps/0/simu");
    // }
}

TEST(TestIODatapipeline, CanWritePredictionData)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    std::string paramsFile = std::string(ROOT_DIR)+"/test/datapipeline/parameters.ini";
    std::string configDir = std::string(ROOT_DIR)+"/test/regression/run1/data";
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/datapipeline/config.yaml";

    std::string timeStamp = logger->getLoggerTime();
    
    {
        IO::IOdatapipeline idp(paramsFile, configDir, "", logger, datapipelineConfig, timeStamp);

        std::vector<Status> statuses = {
            Status{
                .simulation = { 0, 1 },
                .deaths = { 2, 3 },
                .hospital_deaths = { 4, 5},

                .ends = {
                    Compartments{
                        .S = 1014, .E = 1015, .E_t = 1016, .I_p = 1017, .I_t = 1018,
                        .I1 = 1019, .I2 = 1020, .I3 = 1021, .I4 = 1022,
                        .I_s1 = 1023, .I_s2 = 1024, .I_s3 = 1025, .I_s4 = 1026,
                        .H = 1027, .R = 1028, .D = 1029 } },

                .pop_array = { {
                    Compartments{
                        .S = 2014, .E = 2015, .E_t = 2016, .I_p = 2017, .I_t = 2018,
                        .I1 = 2019, .I2 = 2020, .I3 = 2021, .I4 = 2022,
                        .I_s1 = 2023, .I_s2 = 2024, .I_s3 = 2025, .I_s4 = 2026,
                        .H = 2027, .R = 2028, .D = 2029 } } }
            }
        };

        EXPECT_NO_THROW(idp.WritePredictionsToFiles(statuses, "unit_test"));
    }
}