#include <gtest/gtest.h>

#include "butterworth_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

// --- Case 1: N=2, fc=100Hz, fs=1000Hz ----------------------------------------

TEST(Butterworth, N2_fc100_fs1000_Coefficients)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

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

// ============================================================================
// Impulse response tests — Octave impz reference
// ============================================================================

TEST(Butterworth, N2_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(Butterworth, N4_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4> filt(100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(Butterworth, N2_fc500_fs8000_ImpulseResponse)
{
    using Ref = bw_ref::case_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2> filt(500.0, 8000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(Butterworth, N3_fc200_fs4000_ImpulseResponse)
{
    using Ref = bw_ref::case_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3> filt(200.0, 4000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthHPF, N2_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_hp_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthHPF, N4_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_hp_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthHPF, N2_fc500_fs8000_ImpulseResponse)
{
    using Ref = bw_ref::case_hp_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        500.0, 8000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthHPF, N3_fc200_fs4000_ImpulseResponse)
{
    using Ref = bw_ref::case_hp_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3, constfilt::ZOH, constfilt::HighPass> filt(
        200.0, 4000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthMatchedZ, N2_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_mz_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> filt(100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthMatchedZ, N4_fc100_fs1000_ImpulseResponse)
{
    using Ref = bw_ref::case_mz_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4, constfilt::MatchedZ> filt(100.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthMatchedZ, N2_fc500_fs8000_ImpulseResponse)
{
    using Ref = bw_ref::case_mz_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> filt(500.0, 8000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(ButterworthMatchedZ, N3_fc200_fs4000_ImpulseResponse)
{
    using Ref = bw_ref::case_mz_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3, constfilt::MatchedZ> filt(200.0, 4000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

// ============================================================================
// Chirp (frequency sweep) tests — Octave filter(b,a,chirp) reference
// ============================================================================

TEST(Butterworth, N2_fc100_fs1000_Chirp)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(Butterworth, N4_fc100_fs1000_Chirp)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    constfilt::Butterworth<double, 4> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(Butterworth, N2_fc500_fs8000_Chirp)
{
    using Ref = bw_ref::case_2_500Hz_8000Hz;
    constfilt::Butterworth<double, 2> filt(500.0, 8000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(Butterworth, N3_fc200_fs4000_Chirp)
{
    using Ref = bw_ref::case_3_200Hz_4000Hz;
    constfilt::Butterworth<double, 3> filt(200.0, 4000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(ButterworthHPF, N2_fc100_fs1000_Chirp)
{
    using Ref = bw_ref::case_hp_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(ButterworthMatchedZ, N2_fc100_fs1000_Chirp)
{
    using Ref = bw_ref::case_mz_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> filt(100.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}
