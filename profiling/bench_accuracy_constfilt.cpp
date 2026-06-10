// bench_accuracy_constfilt.cpp: constfilt accuracy vs Octave reference (C++17)
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include <cmath>
#include <cstdio>

#include "accuracy_reference.hpp"
#include <constfilt/constfilt.hpp>

static constexpr int kStepLen = 256;

struct AccResult
{
    double max_b_err;
    double max_a_err;
    double max_step_err;
};

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

static void csv_row(const char *ftype, int order, const char *method,
                    const AccResult &r)
{
    std::printf("constfilt,%s,%d,%s,%.6e,%.6e,%.6e\n", ftype, order, method,
                r.max_b_err, r.max_a_err, r.max_step_err);
}

static void human_row(const char *ftype, int order, const char *method,
                      const AccResult &r)
{
    const char *flag = r.max_step_err > 1e-6 ? " !" : "  ";
    std::fprintf(stderr,
                 "%s %-10s %-12s N=%-2d %-10s  b=%.2e  a=%.2e  step=%.2e\n",
                 flag, "constfilt", ftype, order, method, r.max_b_err,
                 r.max_a_err, r.max_step_err);
}

int main()
{
#define BW(N, MC, MSTR)                                                        \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##MSTR##_N##N;                                \
        constfilt::Butterworth<double, N##u, constfilt::MC> f(100.0, 1000.0); \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("butterworth", N, #MSTR, r);                                  \
        human_row("butterworth", N, #MSTR, r);                                \
    } while (0)

    std::fputs("  [constfilt  Butterworth  ZOH]\n", stderr);
    BW(1, ZOH, zoh);   BW(2, ZOH, zoh);   BW(3, ZOH, zoh);  BW(4, ZOH, zoh);
    BW(5, ZOH, zoh);   BW(6, ZOH, zoh);   BW(7, ZOH, zoh);  BW(8, ZOH, zoh);
    BW(9, ZOH, zoh);   BW(10, ZOH, zoh);  BW(11, ZOH, zoh); BW(12, ZOH, zoh);

    std::fputs("  [constfilt  Butterworth  MatchedZ]\n", stderr);
    BW(1, MatchedZ, matchedz);   BW(2, MatchedZ, matchedz);
    BW(3, MatchedZ, matchedz);   BW(4, MatchedZ, matchedz);
    BW(5, MatchedZ, matchedz);   BW(6, MatchedZ, matchedz);
    BW(7, MatchedZ, matchedz);   BW(8, MatchedZ, matchedz);
    BW(9, MatchedZ, matchedz);   BW(10, MatchedZ, matchedz);
    BW(11, MatchedZ, matchedz);  BW(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Butterworth  Tustin]\n", stderr);
    BW(1, Tustin, tustin);   BW(2, Tustin, tustin);
    BW(3, Tustin, tustin);   BW(4, Tustin, tustin);
    BW(5, Tustin, tustin);   BW(6, Tustin, tustin);
    BW(7, Tustin, tustin);   BW(8, Tustin, tustin);
    BW(9, Tustin, tustin);   BW(10, Tustin, tustin);
    BW(11, Tustin, tustin);  BW(12, Tustin, tustin);

#undef BW

#define EL(N, MC, MSTR)                                                        \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##MSTR##_N##N;                                \
        constfilt::Elliptic<double, N##u, constfilt::MC> f(100.0, 0.5, 40.0,  \
                                                           1000.0);           \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("elliptic", N, #MSTR, r);                                     \
        human_row("elliptic", N, #MSTR, r);                                   \
    } while (0)

    std::fputs("  [constfilt  Elliptic  ZOH]\n", stderr);
    EL(2, ZOH, zoh);   EL(3, ZOH, zoh);   EL(4, ZOH, zoh);
    EL(5, ZOH, zoh);   EL(6, ZOH, zoh);   EL(7, ZOH, zoh);
    EL(8, ZOH, zoh);   EL(9, ZOH, zoh);   EL(10, ZOH, zoh);
    EL(11, ZOH, zoh);  EL(12, ZOH, zoh);

    std::fputs("  [constfilt  Elliptic  MatchedZ]\n", stderr);
    EL(2, MatchedZ, matchedz);   EL(3, MatchedZ, matchedz);
    EL(4, MatchedZ, matchedz);   EL(5, MatchedZ, matchedz);
    EL(6, MatchedZ, matchedz);   EL(7, MatchedZ, matchedz);
    EL(8, MatchedZ, matchedz);   EL(9, MatchedZ, matchedz);
    EL(10, MatchedZ, matchedz);  EL(11, MatchedZ, matchedz);
    EL(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Elliptic  Tustin]\n", stderr);
    EL(2, Tustin, tustin);   EL(3, Tustin, tustin);
    EL(4, Tustin, tustin);   EL(5, Tustin, tustin);
    EL(6, Tustin, tustin);   EL(7, Tustin, tustin);
    EL(8, Tustin, tustin);   EL(9, Tustin, tustin);
    EL(10, Tustin, tustin);  EL(11, Tustin, tustin);
    EL(12, Tustin, tustin);

#undef EL

    return 0;
}
