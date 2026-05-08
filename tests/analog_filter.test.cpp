#include <gtest/gtest.h>

#include "continuous_tf_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

using namespace ctf_ref;

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

// --- AnalogFilter: full matrix per case
// -------------------------------------
// Each FULL_MATRIX expansion generates 8 tests: Coefficients, Batch_Step,
// RealTime_Step, Batch_Impulse, RealTime_Impulse, Batch_Chirp, RealTime_Chirp,
// Batch_RealTime_Equivalence.

FULL_MATRIX(AnalogFilter, case_1_zoh_fs10,
            constfilt::AnalogFilter<double, 1u>(Ref::b_s, Ref::a_s,
                                                Ref::sample_rate_hz))

FULL_MATRIX(AnalogFilter, case_2_zoh_fs10,
            constfilt::AnalogFilter<double, 2u>(Ref::b_s, Ref::a_s,
                                                Ref::sample_rate_hz))

FULL_MATRIX(AnalogFilter, case_3_proper_zoh_fs10,
            constfilt::AnalogFilter<double, 1u>(Ref::b_s, Ref::a_s,
                                                Ref::sample_rate_hz))

FULL_MATRIX(AnalogFilter, case_4_mz_fs10,
            constfilt::AnalogFilter<double, 2u, constfilt::MatchedZ>(
                Ref::b_s, Ref::a_s, Ref::sample_rate_hz))

// --- AnalogFilter: stability check disabled allows unstable filter
// ---------------

TEST(AnalogFilter, StabilityCheckDisabled_AllowsUnstableFilter)
{
    // H(s) = 1/(s^2-s+1): unstable, but CheckStab=false allows construction.
    // Arrays must have static storage duration to be usable as NTTPs.
    static constexpr double b[3] = {0.0, 0.0, 1.0};
    static constexpr double a[3] = {1.0, -1.0, 1.0};
    constfilt::AnalogFilter<double, 2u, constfilt::ZOH, false> filt(b, a, 10.0);
    (void)filt;
}

} // namespace
