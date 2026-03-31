#include <gtest/gtest.h>

#include "continuous_tf_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

namespace
{

// --- tf_to_ss: 1st-order strictly proper
// ----------------------------------------
//
// H(s) = 1/(s+1)  ->  b=[0,1], a=[1,1]
//   A=[-1], B=[1], C=[1], D=0

TEST(TfToSs, FirstOrder_StrictlyProper)
{
    constexpr double b[2] = {0.0, 1.0};
    constexpr double a[2] = {1.0, 1.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 1u>(b, a);

    EXPECT_DOUBLE_EQ(sys.A(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(sys.B(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.C(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.D, 0.0);
}

// --- tf_to_ss: 1st-order proper (D != 0)
// ----------------------------------------
//
// H(s) = (s+2)/(s+1)  ->  D=1, residual = 1/(s+1)
//   A=[-1], B=[1], C=[1], D=1

TEST(TfToSs, FirstOrder_Proper)
{
    constexpr double b[2] = {1.0, 2.0};
    constexpr double a[2] = {1.0, 1.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 1u>(b, a);

    EXPECT_DOUBLE_EQ(sys.A(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(sys.B(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.C(0, 0),
                     1.0); // e[1] = b[1]/a[0] - D*a[1]/a[0] = 2-1=1
    EXPECT_DOUBLE_EQ(sys.D, 1.0);
}

// --- tf_to_ss: 2nd-order strictly proper
// ----------------------------------------
//
// H(s) = 1/(s^2+3s+2)  ->  b=[0,0,1], a=[1,3,2]
//   A=[[0,1],[-2,-3]], B=[[0],[1]], C=[[1,0]], D=0

TEST(TfToSs, SecondOrder_StrictlyProper)
{
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, 3.0, 2.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 2u>(b, a);

    EXPECT_DOUBLE_EQ(sys.A(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(sys.A(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(sys.A(1, 0), -2.0);
    EXPECT_DOUBLE_EQ(sys.A(1, 1), -3.0);
    EXPECT_DOUBLE_EQ(sys.B(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(sys.B(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.C(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.C(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(sys.D, 0.0);
}

// --- tf_to_ss: unnormalized denominator
// -----------------------------------------
//
// H(s) = 2/(2s+4) = 1/(s+2): a[0]=2, should normalize correctly

TEST(TfToSs, UnnormalizedDenominator)
{
    constexpr double b[2] = {0.0, 2.0};
    constexpr double a[2] = {2.0, 4.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 1u>(b, a);

    EXPECT_DOUBLE_EQ(sys.A(0, 0), -2.0); // -a[1]/a[0] = -4/2
    EXPECT_DOUBLE_EQ(sys.B(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys.C(0, 0), 1.0); // e[1] = b[1]/a[0] - 0 = 1
    EXPECT_DOUBLE_EQ(sys.D, 0.0);
}

// --- Stability: stable system
// ---------------------------------------------------
//
// H(s) = 1/(s^2+3s+2): poles at -1, -2 (both left half-plane)

TEST(Stability, StableSystem)
{
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, 3.0, 2.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 2u>(b, a);
    constexpr auto stab = constfilt::check_stability(sys);

    static_assert(stab == constfilt::Stability::Stable, "should be stable");
    EXPECT_EQ(stab, constfilt::Stability::Stable);
}

// --- Stability: unstable system
// -------------------------------------------------
//
// H(s) = 1/(s^2-s+1): poles at (1 +/- j*sqrt(3))/2, Re=0.5 > 0

TEST(Stability, UnstableSystem)
{
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, -1.0, 1.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 2u>(b, a);
    constexpr auto stab = constfilt::check_stability(sys);

    static_assert(stab == constfilt::Stability::Unstable, "should be unstable");
    EXPECT_EQ(stab, constfilt::Stability::Unstable);
}

// --- Stability: marginally stable
// -----------------------------------------------
//
// H(s) = 1/(s^2+1): poles at +/-j (simple imaginary-axis poles)

TEST(Stability, MarginallyStable_SimplePoles)
{
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, 0.0, 1.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 2u>(b, a);
    constexpr auto stab = constfilt::check_stability(sys);

    static_assert(stab == constfilt::Stability::MarginallyStable,
                  "simple imaginary-axis poles → MarginallyStable");
    EXPECT_EQ(stab, constfilt::Stability::MarginallyStable);
}

// --- Stability: repeated imaginary-axis pole → Unstable
// -------------------------
//
// H(s) = 1/s^2: double pole at the origin

TEST(Stability, UnstableSystem_RepeatedImaginaryAxisPole)
{
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, 0.0, 0.0};
    constexpr auto sys = constfilt::tf_to_ss<double, 2u>(b, a);
    constexpr auto stab = constfilt::check_stability(sys);

    static_assert(stab == constfilt::Stability::Unstable,
                  "repeated pole at origin → Unstable");
    EXPECT_EQ(stab, constfilt::Stability::Unstable);
}

// --- ContinuousTF: 1st-order, ZOH, coefficients
// ---------------------------------

TEST(ContinuousTF, Case1_ZOH_Coefficients)
{
    using Ref = ctf_ref::case_1_zoh_fs10;
    static constexpr constfilt::ContinuousTF<double, 1u> filt(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);

    EXPECT_NEAR(filt.coeffs_b()[0], Ref::b[0], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_b()[1], Ref::b[1], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_a()[0], Ref::a[0], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_a()[1], Ref::a[1], CONSTFILT_COEFF_TOL);
}

TEST(ContinuousTF, Case1_ZOH_StepResponse)
{
    using Ref = ctf_ref::case_1_zoh_fs10;
    constfilt::ContinuousTF<double, 1u> filt(Ref::b_s, Ref::a_s,
                                             Ref::sample_rate_hz);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- ContinuousTF: 2nd-order, ZOH, coefficients and step response
// ---------------

TEST(ContinuousTF, Case2_ZOH_Coefficients)
{
    using Ref = ctf_ref::case_2_zoh_fs10;
    static constexpr constfilt::ContinuousTF<double, 2u> filt(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ContinuousTF, Case2_ZOH_StepResponse)
{
    using Ref = ctf_ref::case_2_zoh_fs10;
    constfilt::ContinuousTF<double, 2u> filt(Ref::b_s, Ref::a_s,
                                             Ref::sample_rate_hz);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- ContinuousTF: proper TF (D != 0), ZOH
// --------------------------------------

TEST(ContinuousTF, Case3_Proper_ZOH_Coefficients)
{
    using Ref = ctf_ref::case_3_proper_zoh_fs10;
    static constexpr constfilt::ContinuousTF<double, 1u> filt(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);

    EXPECT_NEAR(filt.coeffs_b()[0], Ref::b[0], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_b()[1], Ref::b[1], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_a()[0], Ref::a[0], CONSTFILT_COEFF_TOL);
    EXPECT_NEAR(filt.coeffs_a()[1], Ref::a[1], CONSTFILT_COEFF_TOL);
}

TEST(ContinuousTF, Case3_Proper_ZOH_StepResponse)
{
    using Ref = ctf_ref::case_3_proper_zoh_fs10;
    constfilt::ContinuousTF<double, 1u> filt(Ref::b_s, Ref::a_s,
                                             Ref::sample_rate_hz);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- ContinuousTF: 2nd-order MatchedZ
// -------------------------------------------

TEST(ContinuousTF, Case4_MatchedZ_Coefficients)
{
    using Ref = ctf_ref::case_4_mz_fs10;
    static constexpr constfilt::ContinuousTF<double, 2u, constfilt::MatchedZ>
        filt(Ref::b_s, Ref::a_s, Ref::sample_rate_hz);

    for (unsigned int i = 0; i <= 2u; ++i)
    {
        EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)
            << "b[" << i << "] mismatch";
        EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)
            << "a[" << i << "] mismatch";
    }
}

TEST(ContinuousTF, Case4_MatchedZ_StepResponse)
{
    using Ref = ctf_ref::case_4_mz_fs10;
    constfilt::ContinuousTF<double, 2u, constfilt::MatchedZ> filt(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);

    for (unsigned int i = 0; i < 32u; ++i)
    {
        double y = filt(1.0);
        EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)
            << "step[" << i << "] mismatch";
    }
}

// --- ContinuousTF: stability check disabled allows unstable filter
// ---------------

TEST(ContinuousTF, StabilityCheckDisabled_AllowsUnstableFilter)
{
    // H(s) = 1/(s^2-s+1): unstable, but CheckStab=false allows construction.
    constexpr double b[3] = {0.0, 0.0, 1.0};
    constexpr double a[3] = {1.0, -1.0, 1.0};
    constfilt::ContinuousTF<double, 2u, constfilt::ZOH, false> filt(b, a, 10.0);
    (void)filt;
}

} // namespace
