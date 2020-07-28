#include "gtest/gtest.h"
#include "ModelTypes.h"
#include "IO-datapipeline.h"

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
    EXPECT_ANY_THROW(IO::IOdatapipeline idp{"NoValidPath.ini"};);
}

TEST(TestIODatapipeline, CanReadFixedParameters)
{
    IO::IOdatapipeline idp{"../test/datapipeline/parameters.ini"};
    CommonModelInputParameters params = idp.ReadCommonParameters();

    EXPECT_EQ(params.paramlist.T_lat, 4);
    EXPECT_EQ(params.paramlist.juvp_s, 0.1);
    EXPECT_EQ(params.paramlist.T_inf, 1.5);
    EXPECT_EQ(params.paramlist.T_rec, 11);
    EXPECT_EQ(params.paramlist.T_sym, 7.0);
    EXPECT_EQ(params.paramlist.T_hos, 5);
    EXPECT_EQ(params.paramlist.K, 2000);
    EXPECT_EQ(params.paramlist.inf_asym, 1.0);
}
