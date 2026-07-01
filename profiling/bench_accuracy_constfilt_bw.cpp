// bench_accuracy_constfilt_bw.cpp: Butterworth accuracy rows.
// Compiled as a separate TU from the Elliptic rows to keep per-TU template
// instantiation count (and GCC memory usage) manageable at N=25.

#include "bench_accuracy_common.hpp"
#include <constfilt/constfilt.hpp>

#define RUN_BW(N, MC, mstr)                                                    \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##mstr##_N##N;                                 \
        constfilt::Butterworth<double, N##u, constfilt::MC,                    \
                               constfilt::LowPass, false>                      \
            f(100.0, 1000.0);                                                  \
        auto r = check(f, Ref::b, Ref::a, Ref::step_sos);                      \
        csv_row("constfilt", "butterworth", N, #mstr, r);                      \
        human_row("constfilt", "butterworth", N, #mstr, r);                    \
    } while (0)

#define RUN_BW_SOS(N, MC, mstr)                                                \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##mstr##_N##N;                                 \
        constfilt::Butterworth<double, N##u, constfilt::MC,                    \
                               constfilt::LowPass, true>                       \
            f(100.0, 1000.0);                                                  \
        auto r = check_step([&] { return f(1.0); }, Ref::step_sos);            \
        csv_row("constfilt", "butterworth", N, #mstr "_sos", r);               \
        human_row("constfilt", "butterworth", N, #mstr "_sos", r);             \
    } while (0)

void run_bw_accuracy()
{
    std::fputs("  [constfilt  Butterworth  ZOH  direct]\n", stderr);
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
    RUN_BW(13, ZOH, zoh);
    RUN_BW(14, ZOH, zoh);
    RUN_BW(15, ZOH, zoh);
    RUN_BW(16, ZOH, zoh);
    RUN_BW(17, ZOH, zoh);
    RUN_BW(18, ZOH, zoh);
    RUN_BW(19, ZOH, zoh);
    RUN_BW(20, ZOH, zoh);
    RUN_BW(21, ZOH, zoh);
    RUN_BW(22, ZOH, zoh);
    RUN_BW(23, ZOH, zoh);
    RUN_BW(24, ZOH, zoh);
    RUN_BW(25, ZOH, zoh);

    std::fputs("  [constfilt  Butterworth  MatchedZ  direct]\n", stderr);
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
    RUN_BW(13, MatchedZ, matchedz);
    RUN_BW(14, MatchedZ, matchedz);
    RUN_BW(15, MatchedZ, matchedz);
    RUN_BW(16, MatchedZ, matchedz);
    RUN_BW(17, MatchedZ, matchedz);
    RUN_BW(18, MatchedZ, matchedz);
    RUN_BW(19, MatchedZ, matchedz);
    RUN_BW(20, MatchedZ, matchedz);
    RUN_BW(21, MatchedZ, matchedz);
    RUN_BW(22, MatchedZ, matchedz);
    RUN_BW(23, MatchedZ, matchedz);
    RUN_BW(24, MatchedZ, matchedz);
    RUN_BW(25, MatchedZ, matchedz);

    std::fputs("  [constfilt  Butterworth  Tustin  direct]\n", stderr);
    RUN_BW(1, TustinPW, prewarp);
    RUN_BW(2, TustinPW, prewarp);
    RUN_BW(3, TustinPW, prewarp);
    RUN_BW(4, TustinPW, prewarp);
    RUN_BW(5, TustinPW, prewarp);
    RUN_BW(6, TustinPW, prewarp);
    RUN_BW(7, TustinPW, prewarp);
    RUN_BW(8, TustinPW, prewarp);
    RUN_BW(9, TustinPW, prewarp);
    RUN_BW(10, TustinPW, prewarp);
    RUN_BW(11, TustinPW, prewarp);
    RUN_BW(12, TustinPW, prewarp);
    RUN_BW(13, TustinPW, prewarp);
    RUN_BW(14, TustinPW, prewarp);
    RUN_BW(15, TustinPW, prewarp);
    RUN_BW(16, TustinPW, prewarp);
    RUN_BW(17, TustinPW, prewarp);
    RUN_BW(18, TustinPW, prewarp);
    RUN_BW(19, TustinPW, prewarp);
    RUN_BW(20, TustinPW, prewarp);
    RUN_BW(21, TustinPW, prewarp);
    RUN_BW(22, TustinPW, prewarp);
    RUN_BW(23, TustinPW, prewarp);
    RUN_BW(24, TustinPW, prewarp);
    RUN_BW(25, TustinPW, prewarp);

    std::fputs("  [constfilt  Butterworth  MatchedZ  SOS]\n", stderr);
    RUN_BW_SOS(1, MatchedZ, matchedz);
    RUN_BW_SOS(2, MatchedZ, matchedz);
    RUN_BW_SOS(3, MatchedZ, matchedz);
    RUN_BW_SOS(4, MatchedZ, matchedz);
    RUN_BW_SOS(5, MatchedZ, matchedz);
    RUN_BW_SOS(6, MatchedZ, matchedz);
    RUN_BW_SOS(7, MatchedZ, matchedz);
    RUN_BW_SOS(8, MatchedZ, matchedz);
    RUN_BW_SOS(9, MatchedZ, matchedz);
    RUN_BW_SOS(10, MatchedZ, matchedz);
    RUN_BW_SOS(11, MatchedZ, matchedz);
    RUN_BW_SOS(12, MatchedZ, matchedz);
    RUN_BW_SOS(13, MatchedZ, matchedz);
    RUN_BW_SOS(14, MatchedZ, matchedz);
    RUN_BW_SOS(15, MatchedZ, matchedz);
    RUN_BW_SOS(16, MatchedZ, matchedz);
    RUN_BW_SOS(17, MatchedZ, matchedz);
    RUN_BW_SOS(18, MatchedZ, matchedz);
    RUN_BW_SOS(19, MatchedZ, matchedz);
    RUN_BW_SOS(20, MatchedZ, matchedz);
    RUN_BW_SOS(21, MatchedZ, matchedz);
    RUN_BW_SOS(22, MatchedZ, matchedz);
    RUN_BW_SOS(23, MatchedZ, matchedz);
    RUN_BW_SOS(24, MatchedZ, matchedz);
    RUN_BW_SOS(25, MatchedZ, matchedz);

    std::fputs("  [constfilt  Butterworth  Tustin  SOS]\n", stderr);
    RUN_BW_SOS(1, TustinPW, prewarp);
    RUN_BW_SOS(2, TustinPW, prewarp);
    RUN_BW_SOS(3, TustinPW, prewarp);
    RUN_BW_SOS(4, TustinPW, prewarp);
    RUN_BW_SOS(5, TustinPW, prewarp);
    RUN_BW_SOS(6, TustinPW, prewarp);
    RUN_BW_SOS(7, TustinPW, prewarp);
    RUN_BW_SOS(8, TustinPW, prewarp);
    RUN_BW_SOS(9, TustinPW, prewarp);
    RUN_BW_SOS(10, TustinPW, prewarp);
    RUN_BW_SOS(11, TustinPW, prewarp);
    RUN_BW_SOS(12, TustinPW, prewarp);
    RUN_BW_SOS(13, TustinPW, prewarp);
    RUN_BW_SOS(14, TustinPW, prewarp);
    RUN_BW_SOS(15, TustinPW, prewarp);
    RUN_BW_SOS(16, TustinPW, prewarp);
    RUN_BW_SOS(17, TustinPW, prewarp);
    RUN_BW_SOS(18, TustinPW, prewarp);
    RUN_BW_SOS(19, TustinPW, prewarp);
    RUN_BW_SOS(20, TustinPW, prewarp);
    RUN_BW_SOS(21, TustinPW, prewarp);
    RUN_BW_SOS(22, TustinPW, prewarp);
    RUN_BW_SOS(23, TustinPW, prewarp);
    RUN_BW_SOS(24, TustinPW, prewarp);
    RUN_BW_SOS(25, TustinPW, prewarp);
}
