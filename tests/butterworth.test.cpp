#include <gtest/gtest.h>

#include "butterworth_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

// --- Case 1: N=2, fc=100Hz, fs=1000Hz ----------------------------------------

TEST(Butterworth, N2_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(Butterworth, N2_fc100_fs1000_Batch)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

    // Compile-time batch run + verification.
    static constexpr auto STEP32 = make_step<double, 32>();
    static constexpr auto out = constfilt::batch_filter(filt, STEP32.data);
    static_assert(all_within_tol(out, Ref::step, 1e-7), "step response");

    // Runtime verification.
    double step[32]{};
    for (unsigned int i = 0; i < 32u; ++i)
        step[i] = 1.0;
    double y[32]{};
    filt(step, y);
    for (unsigned int i = 0; i < 32u; ++i)
    {
        EXPECT_NEAR(y[i], Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

TEST(Butterworth, N2_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Case 2: N=4, fc=100Hz, fs=1000Hz ----------------------------------------

TEST(Butterworth, N4_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 4> filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_b()[4], Ref::b[4], 1e-9), "b[4]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");
    static_assert(withinTol(filt.coeffs_a()[4], Ref::a[4], 1e-9), "a[4]");

    for (unsigned int i = 0; i <= 4u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(Butterworth, N4_fc100_fs1000_Batch)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 4> filt(100.0, 1000.0);

    // Compile-time batch run + verification.
    static constexpr auto STEP32 = make_step<double, 32>();
    static constexpr auto out = constfilt::batch_filter(filt, STEP32.data);
    static_assert(all_within_tol(out, Ref::step, 1e-7), "step response");

    // Runtime verification.
    double step[32]{};
    for (unsigned int i = 0; i < 32u; ++i)
        step[i] = 1.0;
    double y[32]{};
    filt(step, y);
    for (unsigned int i = 0; i < 32u; ++i)
    {
        EXPECT_NEAR(y[i], Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

TEST(Butterworth, N4_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Case 3: N=2, fc=500Hz, fs=8000Hz ----------------------------------------

TEST(Butterworth, N2_fc500_fs8000_Coefficients)
{
    using Ref = bw_ref::case_2_500Hz_8000Hz;
    static constexpr constfilt::Butterworth<double, 2> filt(500.0, 8000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(Butterworth, N2_fc500_fs8000_RealTime)
{
    using Ref = bw_ref::case_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2> filt(500.0, 8000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Case 4: N=3, fc=200Hz, fs=4000Hz ----------------------------------------

TEST(Butterworth, N3_fc200_fs4000_Coefficients)
{
    using Ref = bw_ref::case_3_200Hz_4000Hz;
    static constexpr constfilt::Butterworth<double, 3> filt(200.0, 4000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");

    for (unsigned int i = 0; i <= 3u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(Butterworth, N3_fc200_fs4000_RealTime)
{
    using Ref = bw_ref::case_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3> filt(200.0, 4000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HPF: N=2, fc=100Hz, fs=1000Hz -------------------------------------------

TEST(ButterworthHPF, N2_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_hp_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2, constfilt::ZOH,
                                            constfilt::HighPass>
        filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthHPF, N2_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_hp_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HPF: N=4, fc=100Hz, fs=1000Hz -------------------------------------------

TEST(ButterworthHPF, N4_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_hp_4_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 4, constfilt::ZOH,
                                            constfilt::HighPass>
        filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_b()[4], Ref::b[4], 1e-9), "b[4]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");
    static_assert(withinTol(filt.coeffs_a()[4], Ref::a[4], 1e-9), "a[4]");

    for (unsigned int i = 0; i <= 4u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthHPF, N4_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_hp_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HPF: N=2, fc=500Hz, fs=8000Hz -------------------------------------------

TEST(ButterworthHPF, N2_fc500_fs8000_Coefficients)
{
    using Ref = bw_ref::case_hp_2_500Hz_8000Hz;
    static constexpr constfilt::Butterworth<double, 2, constfilt::ZOH,
                                            constfilt::HighPass>
        filt(500.0, 8000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthHPF, N2_fc500_fs8000_RealTime)
{
    using Ref = bw_ref::case_hp_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        500.0, 8000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HPF: N=3, fc=200Hz, fs=4000Hz -------------------------------------------

TEST(ButterworthHPF, N3_fc200_fs4000_Coefficients)
{
    using Ref = bw_ref::case_hp_3_200Hz_4000Hz;
    static constexpr constfilt::Butterworth<double, 3, constfilt::ZOH,
                                            constfilt::HighPass>
        filt(200.0, 4000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");

    for (unsigned int i = 0; i <= 3u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthHPF, N3_fc200_fs4000_RealTime)
{
    using Ref = bw_ref::case_hp_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3, constfilt::ZOH, constfilt::HighPass> filt(
        200.0, 4000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Matched-Z: N=2, fc=100Hz, fs=1000Hz -------------------------------------

TEST(ButterworthMatchedZ, N2_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_mz_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2, constfilt::MatchedZ>
        filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthMatchedZ, N2_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_mz_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Matched-Z: N=4, fc=100Hz, fs=1000Hz -------------------------------------

TEST(ButterworthMatchedZ, N4_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_mz_4_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 4, constfilt::MatchedZ>
        filt(100.0, 1000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_b()[4], Ref::b[4], 1e-9), "b[4]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");
    static_assert(withinTol(filt.coeffs_a()[4], Ref::a[4], 1e-9), "a[4]");

    for (unsigned int i = 0; i <= 4u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthMatchedZ, N4_fc100_fs1000_RealTime)
{
    using Ref = bw_ref::case_mz_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4, constfilt::MatchedZ> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Matched-Z: N=2, fc=500Hz, fs=8000Hz -------------------------------------

TEST(ButterworthMatchedZ, N2_fc500_fs8000_Coefficients)
{
    using Ref = bw_ref::case_mz_2_500Hz_8000Hz;
    static constexpr constfilt::Butterworth<double, 2, constfilt::MatchedZ>
        filt(500.0, 8000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthMatchedZ, N2_fc500_fs8000_RealTime)
{
    using Ref = bw_ref::case_mz_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> filt(500.0, 8000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- Matched-Z: N=3, fc=200Hz, fs=4000Hz -------------------------------------

TEST(ButterworthMatchedZ, N3_fc200_fs4000_Coefficients)
{
    using Ref = bw_ref::case_mz_3_200Hz_4000Hz;
    static constexpr constfilt::Butterworth<double, 3, constfilt::MatchedZ>
        filt(200.0, 4000.0);

    static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
    static_assert(withinTol(filt.coeffs_b()[1], Ref::b[1], 1e-9), "b[1]");
    static_assert(withinTol(filt.coeffs_b()[2], Ref::b[2], 1e-9), "b[2]");
    static_assert(withinTol(filt.coeffs_b()[3], Ref::b[3], 1e-9), "b[3]");
    static_assert(withinTol(filt.coeffs_a()[0], Ref::a[0], 1e-9), "a[0]");
    static_assert(withinTol(filt.coeffs_a()[1], Ref::a[1], 1e-9), "a[1]");
    static_assert(withinTol(filt.coeffs_a()[2], Ref::a[2], 1e-9), "a[2]");
    static_assert(withinTol(filt.coeffs_a()[3], Ref::a[3], 1e-9), "a[3]");

    for (unsigned int i = 0; i <= 3u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ButterworthMatchedZ, N3_fc200_fs4000_RealTime)
{
    using Ref = bw_ref::case_mz_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3, constfilt::MatchedZ> filt(200.0, 4000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}
