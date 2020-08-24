#include "gtest/gtest.h"
#include "ModelTypes.h"
#include "IO-datapipeline.h"
#include "IO.h"

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

    EXPECT_ANY_THROW(IO::IOdatapipeline idp2("../test/datapipeline/parameters.ini", "", logger, "NoValidPath.yaml"););
}

TEST(TestIODatapipeline, CanReadFixedParameters)
{
    const std::string out_dir = std::string(ROOT_DIR)+"/outputs";
    Utilities::logging_stream::Sptr logger = std::make_shared<Utilities::logging_stream>(out_dir);

    IO::IOdatapipeline idp("../test/datapipeline/parameters.ini", "", logger, "../test/datapipeline/config.yaml");
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
    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
    {
        os << "(" << v.size() << ") [";
        for (auto &ve : v) {
            os << " " << ve;
        }
        os << " ]";
        return os;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& v)
    {
        os << "(" << v.size() << ") [\n";
        for (auto &ve : v) {
            os << "  " << ve << "\n";
        }
        os << "]";
        return os;
    }

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
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/test/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, logger, datapipelineConfig);
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
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/test/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, logger, datapipelineConfig);
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
    std::string datapipelineConfig = std::string(ROOT_DIR)+"/test/datapipeline/config.yaml";

    // Load data from data pipeline store
    IO::IOdatapipeline idp(paramsFile, configDir, logger, datapipelineConfig);
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