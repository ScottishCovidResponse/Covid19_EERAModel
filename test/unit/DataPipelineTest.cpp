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
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        os << "(" << v.size() << ") [";
        for (auto &ve : v) {
            os << " " << ve;
        }
        os << " ]";
        return os;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& v) {
        os << "(" << v.size() << ") [\n";
        for (auto &ve : v) {
            os << "  " << ve << "\n";
        }
        os << "]";
        return os;
    }

    void compare_double_eq(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b) {
        EXPECT_EQ(a.size(), b.size());
        std::size_t sz1 = std::min(a.size(), b.size());

        for (int j = 0; j < sz1; ++j) {
            EXPECT_EQ(a[j].size(), b[j].size());
            std::size_t sz2 = std::min(a[j].size(), b[j].size());

            for (int i = 0; i < sz2; ++i) {
                EXPECT_DOUBLE_EQ(a[j][i], b[j][i]);
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

    EXPECT_EQ(dp_params.cases, rg_params.cases);
    compare_double_eq(dp_params.age_pop, rg_params.age_pop);
    compare_double_eq(dp_params.waifw_norm, rg_params.waifw_norm);
    compare_double_eq(dp_params.waifw_home, rg_params.waifw_home);
    compare_double_eq(dp_params.waifw_sdist, rg_params.waifw_sdist);
    compare_double_eq(dp_params.cfr_byage, rg_params.cfr_byage);
}
