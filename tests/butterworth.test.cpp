#include <gtest/gtest.h>

#include "constfilt.hpp"
#include "butterworth_reference.hpp"
#include "test_tools.hpp"

// ─── Helpers ─────────────────────────────────────────────────────────────────

// Build a constexpr step-input array of length N
template <unsigned int N>
static constexpr consteig::Array<double, N> make_step()
{
    consteig::Array<double, N> s{};
    for (unsigned int i = 0; i < N; ++i)
    {
        s[i] = 1.0;
    }
    return s;
}

// ─── Case 1: N=2, fc=100Hz, fs=1000Hz ────────────────────────────────────────

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

TEST(Butterworth, N2_fc100_fs1000_BatchConstexpr)
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 2> filt(100.0, 1000.0);
    static constexpr auto step = make_step<32>();
    static constexpr auto y = filt(step);

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

// ─── Case 2: N=4, fc=100Hz, fs=1000Hz ────────────────────────────────────────

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

TEST(Butterworth, N4_fc100_fs1000_BatchConstexpr)
{
    using Ref = bw_ref::case_4_100Hz_1000Hz;
    static constexpr constfilt::Butterworth<double, 4> filt(100.0, 1000.0);
    static constexpr auto step = make_step<32>();
    static constexpr auto y = filt(step);

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

// ─── Case 3: N=2, fc=500Hz, fs=8000Hz ────────────────────────────────────────

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

// ─── Case 4: N=3, fc=200Hz, fs=4000Hz ────────────────────────────────────────

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
