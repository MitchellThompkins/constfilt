#include <gtest/gtest.h>

#include "constfilt.hpp"
#include "test_tools.hpp"

// ─── ZOH discretization: 1st-order system H(s) = 1/(s+a) ────────────────────
//
// Analytic ZOH result:
//   Ad = exp(-a*T)
//   Bd = (1 - exp(-a*T)) / a
//   Cd = 1,  Dd = 0
//
// Continuous SS:  A = [-a],  B = [1],  C = [1],  D = 0

namespace
{

// ─── a=1, T=0.1 ──────────────────────────────────────────────────────────────

TEST(ZOH, FirstOrder_a1_T0p1)
{
    // Matrix<double,1,1> init: Array<Array<double,1>,1> _data → {{{{val}}}}
    constexpr constfilt::StateSpace<double, 1u> sys_c{
        {{{{-1.0}}}}, // A = [-1]
        {{{{1.0}}}},  // B = [1]
        {{{{1.0}}}},  // C = [1]
        0.0           // D = 0
    };

    constexpr double Ts = 0.1;
    constexpr auto sys_d = constfilt::zoh_discretize(sys_c, Ts, constfilt::ZOH{});

    // Analytic: Ad = exp(-1*0.1) = exp(-0.1) ≈ 0.90483741803595957
    // Bd = (1 - exp(-0.1)) / 1   ≈ 0.09516258196404043
    EXPECT_NEAR(sys_d.A(0, 0), 0.90483741803595957, 1e-10);
    EXPECT_NEAR(sys_d.B(0, 0), 0.09516258196404043, 1e-10);
    EXPECT_DOUBLE_EQ(sys_d.C(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sys_d.D, 0.0);
}

// ─── a=5, T=0.01 ─────────────────────────────────────────────────────────────

TEST(ZOH, FirstOrder_a5_T0p01)
{
    constexpr constfilt::StateSpace<double, 1u> sys_c{
        {{{{-5.0}}}},
        {{{{1.0}}}},
        {{{{1.0}}}},
        0.0
    };

    constexpr double Ts = 0.01;
    constexpr auto sys_d = constfilt::zoh_discretize(sys_c, Ts, constfilt::ZOH{});

    // exp(-5*0.01) = exp(-0.05)
    constexpr double Ad_ref = 0.95122942450071403;
    constexpr double Bd_ref = (1.0 - 0.95122942450071403) / 5.0;

    EXPECT_NEAR(sys_d.A(0, 0), Ad_ref, 1e-10);
    EXPECT_NEAR(sys_d.B(0, 0), Bd_ref, 1e-10);
}

// ─── expm: 1×1 diagonal case ─────────────────────────────────────────────────

TEST(Expm, Scalar)
{
    // expm([-2]) should be [exp(-2)] ≈ 0.13533528323661270
    constexpr consteig::Matrix<double, 1u, 1u> A{{{{-2.0}}}};
    constexpr auto eA = constfilt::expm(A);
    EXPECT_NEAR(eA(0, 0), 0.13533528323661270, 1e-10);
}

// ─── Faddeev-LeVerrier: 2×2 diagonal matrix ──────────────────────────────────

TEST(FaddeevLeVerrier, DiagonalMatrix)
{
    // A = diag(2, 3), char poly = (λ-2)(λ-3) = λ^2 - 5λ + 6
    // Expected: [1, -5, 6]
    constexpr consteig::Matrix<double, 2u, 2u> A{{{{2.0, 0.0}, {0.0, 3.0}}}};
    constexpr auto coeffs = constfilt::faddeev_leverrier(A);

    EXPECT_NEAR(coeffs[0],  1.0, 1e-12);
    EXPECT_NEAR(coeffs[1], -5.0, 1e-12);
    EXPECT_NEAR(coeffs[2],  6.0, 1e-12);
}

// ─── ss_to_tf: 1st-order TF check ────────────────────────────────────────────

TEST(SsToTf, FirstOrder)
{
    // Discrete 1st-order: A=[0.9], B=[0.1], C=[1], D=0
    // Denominator: char poly of [0.9] = z - 0.9 → [1, -0.9]
    // h[0]=0, h[1]=C*B=0.1
    // b[0]=0, b[1]=0.1
    constexpr constfilt::StateSpace<double, 1u> sys{
        {{{{0.9}}}},
        {{{{0.1}}}},
        {{{{1.0}}}},
        0.0
    };
    constexpr auto tf = constfilt::ss_to_tf(sys);

    EXPECT_NEAR(tf.a[0],  1.0, 1e-12);
    EXPECT_NEAR(tf.a[1], -0.9, 1e-12);
    EXPECT_NEAR(tf.b[0],  0.0, 1e-12);
    EXPECT_NEAR(tf.b[1],  0.1, 1e-12);
}

} // namespace
