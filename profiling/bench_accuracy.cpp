// bench_accuracy.cpp — accuracy comparison vs Octave reference (C++17)
//
// Covers constfilt (all three methods, orders 1–12) and iir1 (Tustin/bilinear
// design, orders 1–12 Butterworth).  iir1 and KFR both use bilinear/Tustin
// internally, so they are compared against the Tustin reference.
//
// KFR accuracy lives in bench_accuracy_kfr.cpp (C++20 TU) for the same
// reason as bench_kfr.cpp: consteig vendor headers break under C++20 rules.
//
// CSV to stdout (appended to by bench_accuracy_kfr):
//   library,filter_type,order,method,max_b_err,max_a_err,max_step_err
//   (max_b_err / max_a_err are -1 for libraries that use SOS internally and
//    do not expose direct-form b/a polynomials)
// Human table to stderr.

#include <cmath>
#include <cstdio>

#include "accuracy_reference.hpp"
#include <constfilt/constfilt.hpp>

#ifdef CONSTFILT_BENCH_IIR1
#include <Iir.h>
#endif

static constexpr int kStepLen = 256;

struct AccResult
{
    double max_b_err;
    double max_a_err;
    double max_step_err;
};

// NC = number of coefficients, deduced from the reference arrays.
template <size_t NC, typename F>
static AccResult check(F &filt, const double (&ref_b)[NC],
                       const double (&ref_a)[NC],
                       const double (&ref_step)[kStepLen])
{
    double max_b = 0.0;
    for (size_t i = 0; i < NC; ++i)
    {
        double e = std::fabs(filt.coeffs_b()[i] - ref_b[i]);
        if (e > max_b)
            max_b = e;
    }
    double max_a = 0.0;
    for (size_t i = 0; i < NC; ++i)
    {
        double e = std::fabs(filt.coeffs_a()[i] - ref_a[i]);
        if (e > max_a)
            max_a = e;
    }
    double max_s = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double y = filt(1.0);
        double e = std::fabs(y - ref_step[i]);
        if (e > max_s)
            max_s = e;
    }
    return {max_b, max_a, max_s};
}

static void csv_row(const char *lib, const char *ftype, int order,
                    const char *method, const AccResult &r)
{
    // max_b_err / max_a_err are -1 when the library uses SOS internally and
    // does not expose direct-form polynomial coefficients.
    if (r.max_b_err < 0.0)
        std::printf("%s,%s,%d,%s,n/a,n/a,%.6e\n", lib, ftype, order, method,
                    r.max_step_err);
    else
        std::printf("%s,%s,%d,%s,%.6e,%.6e,%.6e\n", lib, ftype, order, method,
                    r.max_b_err, r.max_a_err, r.max_step_err);
}

static void human_row(const char *lib, const char *ftype, int order,
                      const char *method, const AccResult &r)
{
    const char *flag = r.max_step_err > 1e-6 ? " !" : "  ";
    if (r.max_b_err < 0.0)
        std::fprintf(stderr,
                     "%s %-10s %-12s N=%-2d %-10s  b=n/a       a=n/a       "
                     "step=%.2e\n",
                     flag, lib, ftype, order, method, r.max_step_err);
    else
        std::fprintf(stderr,
                     "%s %-10s %-12s N=%-2d %-10s  b=%.2e  a=%.2e  "
                     "step=%.2e\n",
                     flag, lib, ftype, order, method, r.max_b_err,
                     r.max_a_err, r.max_step_err);
}

// Step-only accuracy check for libraries that use SOS internally.
// max_b_err and max_a_err are set to -1 to signal "not applicable".
template <typename StepFn>
static AccResult check_step(StepFn run_one,
                             const double (&ref_step)[kStepLen])
{
    double max_s = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double y = run_one();
        double e = std::fabs(y - ref_step[i]);
        if (e > max_s)
            max_s = e;
    }
    return {-1.0, -1.0, max_s};
}

int main()
{
    std::puts(
        "library,filter_type,order,method,max_b_err,max_a_err,max_step_err");
    std::fputs(
        "  accuracy vs Octave reference  (!= step err > 1e-6)\n\n"
        "     library    filter_type  order method      "
        "max|b_err|   max|a_err|  max|step_err|\n"
        "  ------------------------------------------------------------------"
        "--------------------\n",
        stderr);

    // ── constfilt Butterworth ─────────────────────────────────────────────────
#define BW(N, MC, MSTR)                                                        \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##MSTR##_N##N;                                \
        constfilt::Butterworth<double, N##u, constfilt::MC> f(100.0, 1000.0); \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("constfilt", "butterworth", N, #MSTR, r);                     \
        human_row("constfilt", "butterworth", N, #MSTR, r);                   \
    } while (0)

    std::fputs("  [Butterworth ZOH]\n", stderr);
    BW(1, ZOH, zoh);
    BW(2, ZOH, zoh);
    BW(3, ZOH, zoh);
    BW(4, ZOH, zoh);
    BW(5, ZOH, zoh);
    BW(6, ZOH, zoh);
    BW(7, ZOH, zoh);
    BW(8, ZOH, zoh);
    BW(9, ZOH, zoh);
    BW(10, ZOH, zoh);
    BW(11, ZOH, zoh);
    BW(12, ZOH, zoh);

    std::fputs("  [Butterworth MatchedZ]\n", stderr);
    BW(1, MatchedZ, matchedz);
    BW(2, MatchedZ, matchedz);
    BW(3, MatchedZ, matchedz);
    BW(4, MatchedZ, matchedz);
    BW(5, MatchedZ, matchedz);
    BW(6, MatchedZ, matchedz);
    BW(7, MatchedZ, matchedz);
    BW(8, MatchedZ, matchedz);
    BW(9, MatchedZ, matchedz);
    BW(10, MatchedZ, matchedz);
    BW(11, MatchedZ, matchedz);
    BW(12, MatchedZ, matchedz);

    std::fputs("  [Butterworth Tustin]\n", stderr);
    BW(1, Tustin, tustin);
    BW(2, Tustin, tustin);
    BW(3, Tustin, tustin);
    BW(4, Tustin, tustin);
    BW(5, Tustin, tustin);
    BW(6, Tustin, tustin);
    BW(7, Tustin, tustin);
    BW(8, Tustin, tustin);
    BW(9, Tustin, tustin);
    BW(10, Tustin, tustin);
    BW(11, Tustin, tustin);
    BW(12, Tustin, tustin);

#undef BW

    // ── constfilt Elliptic ────────────────────────────────────────────────────
#define EL(N, MC, MSTR)                                                        \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##MSTR##_N##N;                                \
        constfilt::Elliptic<double, N##u, constfilt::MC> f(100.0, 0.5, 40.0,  \
                                                           1000.0);           \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("constfilt", "elliptic", N, #MSTR, r);                        \
        human_row("constfilt", "elliptic", N, #MSTR, r);                      \
    } while (0)

    std::fputs("  [Elliptic ZOH]\n", stderr);
    EL(2, ZOH, zoh);
    EL(3, ZOH, zoh);
    EL(4, ZOH, zoh);
    EL(5, ZOH, zoh);
    EL(6, ZOH, zoh);
    EL(7, ZOH, zoh);
    EL(8, ZOH, zoh);
    EL(9, ZOH, zoh);
    EL(10, ZOH, zoh);
    EL(11, ZOH, zoh);
    EL(12, ZOH, zoh);

    std::fputs("  [Elliptic MatchedZ]\n", stderr);
    EL(2, MatchedZ, matchedz);
    EL(3, MatchedZ, matchedz);
    EL(4, MatchedZ, matchedz);
    EL(5, MatchedZ, matchedz);
    EL(6, MatchedZ, matchedz);
    EL(7, MatchedZ, matchedz);
    EL(8, MatchedZ, matchedz);
    EL(9, MatchedZ, matchedz);
    EL(10, MatchedZ, matchedz);
    EL(11, MatchedZ, matchedz);
    EL(12, MatchedZ, matchedz);

    std::fputs("  [Elliptic Tustin]\n", stderr);
    EL(2, Tustin, tustin);
    EL(3, Tustin, tustin);
    EL(4, Tustin, tustin);
    EL(5, Tustin, tustin);
    EL(6, Tustin, tustin);
    EL(7, Tustin, tustin);
    EL(8, Tustin, tustin);
    EL(9, Tustin, tustin);
    EL(10, Tustin, tustin);
    EL(11, Tustin, tustin);
    EL(12, Tustin, tustin);

#undef EL

#ifdef CONSTFILT_BENCH_IIR1
    // ── iir1 Butterworth (Tustin/bilinear) ────────────────────────────────────
    // iir1 uses pre-warped bilinear (exact digital cutoff).  constfilt and
    // Octave's c2d(...,'tustin') use standard bilinear (no pre-warping).  This
    // causes a systematic ~N*1% step-response difference that grows with order
    // but does NOT represent numerical breakdown -- both implementations are
    // self-consistent.  Use iir1's errors here as a baseline for the design
    // convention delta; compare them against constfilt_tustin rows for context.
    std::fputs("  [iir1  Butterworth  tustin(bilinear) -- pre-warped, see note]\n",
               stderr);

#define IIR1_BW(N)                                                             \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_tustin_N##N;                                  \
        Iir::Butterworth::LowPass<N> f;                                        \
        f.setup(1000.0, 100.0);                                                \
        auto r = check_step([&]() { return f.filter(1.0); }, Ref::step);      \
        csv_row("iir1", "butterworth", N, "tustin", r);                       \
        human_row("iir1", "butterworth", N, "tustin", r);                     \
    } while (0)

    IIR1_BW(1);
    IIR1_BW(2);
    IIR1_BW(3);
    IIR1_BW(4);
    IIR1_BW(5);
    IIR1_BW(6);
    IIR1_BW(7);
    IIR1_BW(8);
    IIR1_BW(9);
    IIR1_BW(10);
    IIR1_BW(11);
    IIR1_BW(12);

#undef IIR1_BW
#endif // CONSTFILT_BENCH_IIR1

    return 0;
}
