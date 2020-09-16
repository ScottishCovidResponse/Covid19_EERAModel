#include "gtest/gtest.h"
#include "IO.h"
#include <sstream>
#include <string>

using namespace EERAModel;

TEST(IOTest, WritePredictionFullHeader)
{
    const std::string expected("iter, day, age_group, S, E, E_t, I_p, I_t,"
                                " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D\n");
    
    std::ostringstream os;
    IO::WritePredictionFullHeader(os);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WritePredictionFullRow)
{
    int iter = 0;
    int day = 1;
    int age_group = 3;
    Compartments comp;
    comp.S = 1; comp.E = 2; comp.E_t = 3; comp.I_p = 4; comp.I_t = 5; comp.I1 = 6; comp.I2 = 7;
    comp.I3 = 8; comp.I4 = 9; comp.I_s1 = 10; comp.I_s2 = 11; comp.I_s3 = 12; comp.I_s4 = 13;
    comp.H = 14; comp.R = 15; comp.D = 16;
    const std::string expected("0, 1, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16\n");
    
    std::ostringstream os;
    IO::WritePredictionFullRow(os, iter, day, age_group, comp);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WriteSimuHeader)
{
    const std::string expected("iter, day, inc_case, inc_death_hospital, inc_death\n");

    std::ostringstream os;
    IO::WriteSimuHeader(os);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WriteSimuRow)
{
    int iter = 0;
    int day = 1;
    int inc_case = 2;
    int inc_hospital = 3;
    int inc_death = 4;
    const std::string expected("0, 1, 2, 3, 4\n");

    std::ostringstream os;
    IO::WriteSimuRow(os, iter, day, inc_case, inc_hospital, inc_death);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WriteInferenceEndsHeader)
{
    const std::string expected("iter, age_group, S, E, E_t, I_p, I_t,"
                                " I1, I2, I3, I4, I_s1, I_s2, I_s3, I_s4, H, R, D\n");
    
    std::ostringstream os;
    IO::WriteInferenceEndsHeader(os);

    EXPECT_EQ(expected, os.str());
}


TEST(IOTest, WriteInferenceParticlesHeader)
{
    const std::string expected("iter, nsse_cases, nsse_deaths, p_inf, p_hcw, c_hcw, d, q, p_s,"
        " rrd, lambda, weight\n");

    std::ostringstream os;
    IO::WriteInferenceParticlesHeader(os);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WriteInferenceParticlesRow)
{
    int iter = 0;
    particle p;
    p.nsse_cases = 1;
    p.nsse_deaths = 2;
    p.parameter_set = std::vector<double>{3, 4, 5, 6, 7, 8, 9, 10};
    p.weight  = 0.5;
    const std::string expected("0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0.5\n");

    std::ostringstream os;
    IO::WriteInferenceParticlesRow(os, iter, p);

    EXPECT_EQ(expected, os.str());
}

TEST(IOTest, WriteInferenceEndsRow)
{
    int iter = 0;
    int age_group = 1;
    Compartments comp;
    comp.S = 1; comp.E = 2; comp.E_t = 3; comp.I_p = 4; comp.I_t = 5; comp.I1 = 6; comp.I2 = 7;
    comp.I3 = 8; comp.I4 = 9; comp.I_s1 = 10; comp.I_s2 = 11; comp.I_s3 = 12; comp.I_s4 = 13;
    comp.H = 14; comp.R = 15; comp.D = 16;
    const std::string expected("0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16\n");
    
    std::ostringstream os;
    IO::WriteInferenceEndsRow(os, iter, age_group, comp);

    EXPECT_EQ(expected, os.str());
}