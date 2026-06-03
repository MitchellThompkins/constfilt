#ifndef CONSTFILT_TEST_TOOLS_HPP
#define CONSTFILT_TEST_TOOLS_HPP

// Absolute tolerance for comparing b/a transfer-function coefficients
// against the Octave ZOH reference. Raised to 1e-7 to accommodate the
// accumulated FP error in N=8 coefficient computation (~3e-8 observed).
#ifndef CONSTFILT_COEFF_TOL
#define CONSTFILT_COEFF_TOL 1e-7
#endif

// Absolute tolerance for comparing step-response samples against
// Octave's filter() output.
#ifndef CONSTFILT_STEP_TOL
#define CONSTFILT_STEP_TOL 1e-7
#endif

#include <cmath>

// Absolute-value comparison helper (no std::abs in constexpr context).
template <typename T> static constexpr bool withinTol(T a, T b, T tol)
{
    T diff = a - b;
    return (diff < static_cast<T>(0) ? -diff : diff) < tol;
}

// =============================================================================
// FULL_MATRIX test macros
//
// Generates the canonical eight-test matrix for one filter case:
//   _Coefficients, _Batch_Step,    _RealTime_Step,
//                  _Batch_Impulse, _RealTime_Impulse,
//                  _Batch_Chirp,   _RealTime_Chirp,
//                  _Batch_RealTime_Equivalence
//
// Arguments:
//   suite - gtest suite name (e.g. EllipticLP)
//   ref   - reference struct identifier (e.g. lp_2_5rp_40rs_100Hz_1000Hz)
//   ...   - prvalue expression that constructs the filter, e.g.
//             constfilt::Butterworth<double, 2>(100.0, 1000.0)
//           Variadic so commas in template args / ctor args are forwarded
//           verbatim. May reference `Ref::...` (the macros expose `using Ref =
//           ref;`). Must be side-effect-free: each sub-test constructs its own
//           filter instance (real-time tests mutate state), so the expression
//           is re-expanded once per test (8x by FULL_MATRIX, 2x within
//           BATCH_REALTIME_EQUIVALENCE_TEST).
//
// The expression is used in two contexts:
//   static constexpr auto filt = <expr>;   // constexpr ctor + prvalue elision
//   auto filt = <expr>;                    // runtime ctor for mutable state
//
// Reference struct must expose:
//   static constexpr double a[N+1], b[N+1];
//   static constexpr double step[], impulse[];
//   static constexpr double chirp_in[], chirp[];
// =============================================================================

#define CONSTFILT_DETAIL_LEN_OF(arr) (sizeof(arr) / sizeof((arr)[0]))

#define COEFF_TEST(suite, ref, ...)                                            \
    TEST(suite, ref##_Coefficients)                                            \
    {                                                                          \
        using Ref = ref;                                                       \
        static constexpr auto filt = __VA_ARGS__;                              \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::a);              \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(filt.coeffs_b()[i], Ref::b[i], CONSTFILT_COEFF_TOL)    \
                << "b[" << i << "]";                                           \
            EXPECT_NEAR(filt.coeffs_a()[i], Ref::a[i], CONSTFILT_COEFF_TOL)    \
                << "a[" << i << "]";                                           \
        }                                                                      \
    }

#define BATCH_STEP_TEST(suite, ref, ...)                                       \
    TEST(suite, ref##_Batch_Step)                                              \
    {                                                                          \
        using Ref = ref;                                                       \
        static constexpr auto filt = __VA_ARGS__;                              \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::step);           \
        double input[LEN]{};                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            input[i] = 1.0;                                                    \
        }                                                                      \
        double output[LEN]{};                                                  \
        filt(input, output);                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(output[i], Ref::step[i], CONSTFILT_STEP_TOL)           \
                << "step[" << i << "]";                                        \
        }                                                                      \
    }

#define REALTIME_STEP_TEST(suite, ref, ...)                                    \
    TEST(suite, ref##_RealTime_Step)                                           \
    {                                                                          \
        using Ref = ref;                                                       \
        auto filt = __VA_ARGS__;                                               \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::step);           \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            double y = filt(1.0);                                              \
            EXPECT_NEAR(y, Ref::step[i], CONSTFILT_STEP_TOL)                   \
                << "step[" << i << "]";                                        \
        }                                                                      \
    }

#define BATCH_IMPULSE_TEST(suite, ref, ...)                                    \
    TEST(suite, ref##_Batch_Impulse)                                           \
    {                                                                          \
        using Ref = ref;                                                       \
        static constexpr auto filt = __VA_ARGS__;                              \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::impulse);        \
        double input[LEN]{};                                                   \
        input[0] = 1.0;                                                        \
        double output[LEN]{};                                                  \
        filt(input, output);                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(output[i], Ref::impulse[i], CONSTFILT_STEP_TOL)        \
                << "impulse[" << i << "]";                                     \
        }                                                                      \
    }

#define REALTIME_IMPULSE_TEST(suite, ref, ...)                                 \
    TEST(suite, ref##_RealTime_Impulse)                                        \
    {                                                                          \
        using Ref = ref;                                                       \
        auto filt = __VA_ARGS__;                                               \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::impulse);        \
        EXPECT_NEAR(filt(1.0), Ref::impulse[0], CONSTFILT_STEP_TOL)            \
            << "impulse[0]";                                                   \
        for (unsigned i = 1; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(filt(0.0), Ref::impulse[i], CONSTFILT_STEP_TOL)        \
                << "impulse[" << i << "]";                                     \
        }                                                                      \
    }

#define BATCH_CHIRP_TEST(suite, ref, ...)                                      \
    TEST(suite, ref##_Batch_Chirp)                                             \
    {                                                                          \
        using Ref = ref;                                                       \
        static constexpr auto filt = __VA_ARGS__;                              \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::chirp);          \
        double input[LEN]{};                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            input[i] = Ref::chirp_in[i];                                       \
        }                                                                      \
        double output[LEN]{};                                                  \
        filt(input, output);                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(output[i], Ref::chirp[i], CONSTFILT_STEP_TOL)          \
                << "chirp[" << i << "]";                                       \
        }                                                                      \
    }

#define REALTIME_CHIRP_TEST(suite, ref, ...)                                   \
    TEST(suite, ref##_RealTime_Chirp)                                          \
    {                                                                          \
        using Ref = ref;                                                       \
        auto filt = __VA_ARGS__;                                               \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::chirp);          \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            EXPECT_NEAR(filt(Ref::chirp_in[i]), Ref::chirp[i],                 \
                        CONSTFILT_STEP_TOL)                                    \
                << "chirp[" << i << "]";                                       \
        }                                                                      \
    }

// Cross-check: drive the same filter through both paths on the chirp signal
// and assert the two outputs are bit-identical. Catches state-handling
// divergence between batch and real-time that the Octave reference comparison
// would miss.
#define BATCH_REALTIME_EQUIVALENCE_TEST(suite, ref, ...)                       \
    TEST(suite, ref##_Batch_RealTime_Equivalence)                              \
    {                                                                          \
        using Ref = ref;                                                       \
        constexpr unsigned LEN = CONSTFILT_DETAIL_LEN_OF(Ref::chirp);          \
        double input[LEN]{};                                                   \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            input[i] = Ref::chirp_in[i];                                       \
        }                                                                      \
        double batch_out[LEN]{};                                               \
        {                                                                      \
            auto filt = __VA_ARGS__;                                           \
            filt(input, batch_out);                                            \
        }                                                                      \
        auto rt_filt = __VA_ARGS__;                                            \
        for (unsigned i = 0; i < LEN; ++i)                                     \
        {                                                                      \
            double rt = rt_filt(input[i]);                                     \
            EXPECT_EQ(rt, batch_out[i]) << "i=" << i;                          \
        }                                                                      \
    }

#define FULL_MATRIX(suite, ref, ...)                                           \
    COEFF_TEST(suite, ref, __VA_ARGS__)                                        \
    BATCH_STEP_TEST(suite, ref, __VA_ARGS__)                                   \
    REALTIME_STEP_TEST(suite, ref, __VA_ARGS__)                                \
    BATCH_IMPULSE_TEST(suite, ref, __VA_ARGS__)                                \
    REALTIME_IMPULSE_TEST(suite, ref, __VA_ARGS__)                             \
    BATCH_CHIRP_TEST(suite, ref, __VA_ARGS__)                                  \
    REALTIME_CHIRP_TEST(suite, ref, __VA_ARGS__)                               \
    BATCH_REALTIME_EQUIVALENCE_TEST(suite, ref, __VA_ARGS__)

#endif // CONSTFILT_TEST_TOOLS_HPP
