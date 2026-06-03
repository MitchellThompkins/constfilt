#include <gtest/gtest.h>

#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

// --- Mode 3: construct Filter directly from known b/a coefficients -----------

// Reference: a 3-tap FIR-like IIR with known hand-computed output.
// b = [0.25, 0.5, 0.25], a = [1.0, -0.5, 0.0625]
// We verify the DF2T state machine against a step-by-step reference.

namespace
{

// Known transfer function (order 2, NB=3, NA=3)
static constexpr double REF_B[3] = {0.25, 0.5, 0.25};
static constexpr double REF_A[3] = {1.0, -0.5, 0.0625};

// A Filter subclass that exposes the constructor for testing
class TestFilter : public constfilt::Filter<double, 3u, 3u>
{
  public:
    constexpr TestFilter(const double (&b)[3u], const double (&a)[3u])
        : constfilt::Filter<double, 3u, 3u>(b, a)
    {
    }
};

// --- Coefficient accessor tests
// -----------------------------------------------

TEST(FilterCoeffs, StoresB)
{
    static constexpr TestFilter filt(REF_B, REF_A);
    static constexpr auto b = filt.coeffs_b();
    EXPECT_DOUBLE_EQ(b[0], 0.25);
    EXPECT_DOUBLE_EQ(b[1], 0.5);
    EXPECT_DOUBLE_EQ(b[2], 0.25);
}

TEST(FilterCoeffs, StoresA)
{
    static constexpr TestFilter filt(REF_B, REF_A);
    static constexpr auto a = filt.coeffs_a();
    EXPECT_DOUBLE_EQ(a[0], 1.0);
    EXPECT_DOUBLE_EQ(a[1], -0.5);
    EXPECT_DOUBLE_EQ(a[2], 0.0625);
}

// --- Batch operator() correctness --------------------------------------------

// Manual DF2T trace for a step input of 4 samples:
// M = max(3,3)-1 = 2 (state has 2 elements)
//
// n=0: x=1, y = 0.25*1 + s[0]=0.25
//      s[0] = 0.5*1 - (-0.5)*0.25 + s[1] = 0.5 + 0.125 + 0 = 0.625
//      s[1] = 0.25*1 - 0.0625*0.25 = 0.25 - 0.015625 = 0.234375
// n=1: x=1, y = 0.25 + 0.625 = 0.875
//      s[0] = 0.5 - (-0.5)*0.875 + 0.234375 = 0.5 + 0.4375 + 0.234375
//      = 1.171875 s[1] = 0.25 - 0.0625*0.875 = 0.25 - 0.05468750 = 0.195312500
// n=2: x=1, y = 0.25 + 1.171875 = 1.421875
//      s[0] = 0.5 - (-0.5)*1.421875 + 0.195312500 = 0.5+0.710937500+0.195312500
//      = 1.406250 s[1] = 0.25 - 0.0625*1.421875 = 0.25 - 0.088867188 =
//      0.161132812
// n=3: x=1, y = 0.25 + 1.406250 = 1.656250

static constexpr double STEP4[4] = {1.0, 1.0, 1.0, 1.0};

TEST(FilterBatch, StepResponse4)
{
    static constexpr TestFilter filt(REF_B, REF_A);
    double y[4]{};
    filt(STEP4, y);

    // Tolerances: these are exact rational arithmetic values
    EXPECT_NEAR(y[0], 0.25, 1e-14);
    EXPECT_NEAR(y[1], 0.875, 1e-14);
    EXPECT_NEAR(y[2], 1.421875, 1e-14);
    EXPECT_NEAR(y[3], 1.656250, 1e-13);
}

// --- Real-time operator() matches batch operator() ---------------------------

TEST(FilterRealTime, MatchesBatch)
{
    static constexpr double step8[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    // batch path
    static constexpr TestFilter cfilt(REF_B, REF_A);
    double batch_out[8]{};
    cfilt(step8, batch_out);

    // real-time (mutating state)
    TestFilter rt_filt(REF_B, REF_A);
    for (int i = 0; i < 8; ++i)
    {
        double y = rt_filt(1.0);
        EXPECT_NEAR(y, batch_out[static_cast<unsigned int>(i)], 1e-13)
            << "mismatch at sample " << i;
    }
}

// --- Reset clears state
// -------------------------------------------------------

TEST(FilterReset, StateIsZeroedAfterReset)
{
    TestFilter filt(REF_B, REF_A);
    // Drive it for a few samples
    filt(1.0);
    filt(1.0);
    filt(1.0);
    filt.reset();
    // First sample after reset should equal first sample from fresh filter
    TestFilter fresh(REF_B, REF_A);
    EXPECT_DOUBLE_EQ(filt(1.0), fresh(1.0));
}

} // namespace
