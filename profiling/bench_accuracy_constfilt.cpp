// bench_accuracy_constfilt.cpp: constfilt accuracy vs Octave reference (C++17)
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_accuracy_common.hpp"
#include <constfilt/constfilt.hpp>

#define RUN_BW(N, MC, mstr)                                                    \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##mstr##_N##N;                                \
        constfilt::Butterworth<double, N##u, constfilt::MC> f(100.0, 1000.0); \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("constfilt", "butterworth", N, #mstr, r);                     \
        human_row("constfilt", "butterworth", N, #mstr, r);                   \
    } while (0)

#define RUN_EL(N, MC, mstr)                                                    \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##mstr##_N##N;                                \
        constfilt::Elliptic<double, N##u, constfilt::MC> f(100.0, 0.5, 40.0,  \
                                                           1000.0);           \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                         \
        csv_row("constfilt", "elliptic", N, #mstr, r);                        \
        human_row("constfilt", "elliptic", N, #mstr, r);                      \
    } while (0)

int main()
{
    std::fputs("  [constfilt  Butterworth  ZOH]\n", stderr);
    RUN_BW(1, ZOH, zoh);
    RUN_BW(2, ZOH, zoh);
    RUN_BW(3, ZOH, zoh);
    RUN_BW(4, ZOH, zoh);
    RUN_BW(5, ZOH, zoh);
    RUN_BW(6, ZOH, zoh);
    RUN_BW(7, ZOH, zoh);
    RUN_BW(8, ZOH, zoh);
    RUN_BW(9, ZOH, zoh);
    RUN_BW(10, ZOH, zoh);
    RUN_BW(11, ZOH, zoh);
    RUN_BW(12, ZOH, zoh);

    std::fputs("  [constfilt  Butterworth  MatchedZ]\n", stderr);
    RUN_BW(1, MatchedZ, matchedz);
    RUN_BW(2, MatchedZ, matchedz);
    RUN_BW(3, MatchedZ, matchedz);
    RUN_BW(4, MatchedZ, matchedz);
    RUN_BW(5, MatchedZ, matchedz);
    RUN_BW(6, MatchedZ, matchedz);
    RUN_BW(7, MatchedZ, matchedz);
    RUN_BW(8, MatchedZ, matchedz);
    RUN_BW(9, MatchedZ, matchedz);
    RUN_BW(10, MatchedZ, matchedz);
    RUN_BW(11, MatchedZ, matchedz);
    RUN_BW(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Butterworth  Tustin]\n", stderr);
    RUN_BW(1, Tustin, tustin);
    RUN_BW(2, Tustin, tustin);
    RUN_BW(3, Tustin, tustin);
    RUN_BW(4, Tustin, tustin);
    RUN_BW(5, Tustin, tustin);
    RUN_BW(6, Tustin, tustin);
    RUN_BW(7, Tustin, tustin);
    RUN_BW(8, Tustin, tustin);
    RUN_BW(9, Tustin, tustin);
    RUN_BW(10, Tustin, tustin);
    RUN_BW(11, Tustin, tustin);
    RUN_BW(12, Tustin, tustin);

    std::fputs("  [constfilt  Elliptic  ZOH]\n", stderr);
    RUN_EL(2, ZOH, zoh);
    RUN_EL(3, ZOH, zoh);
    RUN_EL(4, ZOH, zoh);
    RUN_EL(5, ZOH, zoh);
    RUN_EL(6, ZOH, zoh);
    RUN_EL(7, ZOH, zoh);
    RUN_EL(8, ZOH, zoh);
    RUN_EL(9, ZOH, zoh);
    RUN_EL(10, ZOH, zoh);
    RUN_EL(11, ZOH, zoh);
    RUN_EL(12, ZOH, zoh);

    std::fputs("  [constfilt  Elliptic  MatchedZ]\n", stderr);
    RUN_EL(2, MatchedZ, matchedz);
    RUN_EL(3, MatchedZ, matchedz);
    RUN_EL(4, MatchedZ, matchedz);
    RUN_EL(5, MatchedZ, matchedz);
    RUN_EL(6, MatchedZ, matchedz);
    RUN_EL(7, MatchedZ, matchedz);
    RUN_EL(8, MatchedZ, matchedz);
    RUN_EL(9, MatchedZ, matchedz);
    RUN_EL(10, MatchedZ, matchedz);
    RUN_EL(11, MatchedZ, matchedz);
    RUN_EL(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Elliptic  Tustin]\n", stderr);
    RUN_EL(2, Tustin, tustin);
    RUN_EL(3, Tustin, tustin);
    RUN_EL(4, Tustin, tustin);
    RUN_EL(5, Tustin, tustin);
    RUN_EL(6, Tustin, tustin);
    RUN_EL(7, Tustin, tustin);
    RUN_EL(8, Tustin, tustin);
    RUN_EL(9, Tustin, tustin);
    RUN_EL(10, Tustin, tustin);
    RUN_EL(11, Tustin, tustin);
    RUN_EL(12, Tustin, tustin);

    return 0;
}
