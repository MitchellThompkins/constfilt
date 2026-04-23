#include <gtest/gtest.h>

#include "elliptic_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

// --- LP: N=2, Rp=0.5, Rs=40, fc=100Hz, fs=1000Hz ----------------------------

TEST(EllipticLP, N2_Rp05_Rs40_fc100_fs1000_Coefficients)
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    static constexpr constfilt::Elliptic<double, 2> filt(100.0, 0.5, 40.0,
                                                         1000.0);

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(EllipticLP, N2_Rp05_Rs40_fc100_fs1000_Batch)
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    static constexpr constfilt::Elliptic<double, 2> filt(100.0, 0.5, 40.0,
                                                         1000.0);
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

TEST(EllipticLP, N2_Rp05_Rs40_fc100_fs1000_RealTime)
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2> filt(100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- LP: N=4, Rp=0.5, Rs=40, fc=100Hz, fs=1000Hz ----------------------------

TEST(EllipticLP, N4_Rp05_Rs40_fc100_fs1000_Coefficients)
{
    using Ref = el_ref::lp_4_5rp_40rs_100Hz_1000Hz;
    static constexpr constfilt::Elliptic<double, 4> filt(100.0, 0.5, 40.0,
                                                         1000.0);

    for (unsigned int i = 0; i <= 4u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(EllipticLP, N4_Rp05_Rs40_fc100_fs1000_Batch)
{
    using Ref = el_ref::lp_4_5rp_40rs_100Hz_1000Hz;
    static constexpr constfilt::Elliptic<double, 4> filt(100.0, 0.5, 40.0,
                                                         1000.0);
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

TEST(EllipticLP, N4_Rp05_Rs40_fc100_fs1000_RealTime)
{
    using Ref = el_ref::lp_4_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 4> filt(100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- LP: N=3, Rp=1.0, Rs=60, fc=200Hz, fs=4000Hz ----------------------------

TEST(EllipticLP, N3_Rp10_Rs60_fc200_fs4000_Coefficients)
{
    using Ref = el_ref::lp_3_10rp_60rs_200Hz_4000Hz;
    static constexpr constfilt::Elliptic<double, 3> filt(200.0, 1.0, 60.0,
                                                         4000.0);

    for (unsigned int i = 0; i <= 3u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(EllipticLP, N3_Rp10_Rs60_fc200_fs4000_RealTime)
{
    using Ref = el_ref::lp_3_10rp_60rs_200Hz_4000Hz;
    constfilt::Elliptic<double, 3> filt(200.0, 1.0, 60.0, 4000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HP: N=2, Rp=0.5, Rs=40, fc=100Hz, fs=1000Hz ----------------------------

TEST(EllipticHP, N2_Rp05_Rs40_fc100_fs1000_Coefficients)
{
    using Ref = el_ref::hp_2_5rp_40rs_100Hz_1000Hz;
    static constexpr constfilt::Elliptic<double, 2, constfilt::ZOH,
                                         constfilt::HighPass>
        filt(100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(EllipticHP, N2_Rp05_Rs40_fc100_fs1000_RealTime)
{
    using Ref = el_ref::hp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- HP: N=3, Rp=1.0, Rs=60, fc=200Hz, fs=4000Hz ----------------------------

TEST(EllipticHP, N3_Rp10_Rs60_fc200_fs4000_Coefficients)
{
    using Ref = el_ref::hp_3_10rp_60rs_200Hz_4000Hz;
    static constexpr constfilt::Elliptic<double, 3, constfilt::ZOH,
                                         constfilt::HighPass>
        filt(200.0, 1.0, 60.0, 4000.0);

    for (unsigned int i = 0; i <= 3u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(EllipticHP, N3_Rp10_Rs60_fc200_fs4000_RealTime)
{
    using Ref = el_ref::hp_3_10rp_60rs_200Hz_4000Hz;
    constfilt::Elliptic<double, 3, constfilt::ZOH, constfilt::HighPass> filt(
        200.0, 1.0, 60.0, 4000.0);

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

TEST(EllipticLP, N2_Rp05_Rs40_fc100_fs1000_ImpulseResponse)
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2> filt(100.0, 0.5, 40.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(EllipticLP, N4_Rp05_Rs40_fc100_fs1000_ImpulseResponse)
{
    using Ref = el_ref::lp_4_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 4> filt(100.0, 0.5, 40.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(EllipticLP, N3_Rp10_Rs60_fc200_fs4000_ImpulseResponse)
{
    using Ref = el_ref::lp_3_10rp_60rs_200Hz_4000Hz;
    constfilt::Elliptic<double, 3> filt(200.0, 1.0, 60.0, 4000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(EllipticHP, N2_Rp05_Rs40_fc100_fs1000_ImpulseResponse)
{
    using Ref = el_ref::hp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 0.5, 40.0, 1000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

TEST(EllipticHP, N3_Rp10_Rs60_fc200_fs4000_ImpulseResponse)
{
    using Ref = el_ref::hp_3_10rp_60rs_200Hz_4000Hz;
    constfilt::Elliptic<double, 3, constfilt::ZOH, constfilt::HighPass> filt(
        200.0, 1.0, 60.0, 4000.0);

    EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL) << "impulse[0]";
    for (unsigned int i = 1; i < 32u; ++i)
        EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)
            << "impulse[" << i << "]";
}

// ============================================================================
// Chirp (frequency sweep) tests — Octave filter(b,a,chirp) reference
// ============================================================================

TEST(EllipticLP, N2_Rp05_Rs40_fc100_fs1000_Chirp)
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2> filt(100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(EllipticLP, N4_Rp05_Rs40_fc100_fs1000_Chirp)
{
    using Ref = el_ref::lp_4_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 4> filt(100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(EllipticLP, N3_Rp10_Rs60_fc200_fs4000_Chirp)
{
    using Ref = el_ref::lp_3_10rp_60rs_200Hz_4000Hz;
    constfilt::Elliptic<double, 3> filt(200.0, 1.0, 60.0, 4000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}

TEST(EllipticHP, N2_Rp05_Rs40_fc100_fs1000_Chirp)
{
    using Ref = el_ref::hp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> filt(
        100.0, 0.5, 40.0, 1000.0);

    for (unsigned int i = 0; i < 256u; ++i)
        EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i], CONSTFILT_STEP_TOL)
            << "chirp[" << i << "]";
}
