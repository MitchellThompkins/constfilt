// bench_accuracy_constfilt.cpp: constfilt accuracy vs Octave reference (C++17)
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_accuracy_common.hpp"
#include <constfilt/constfilt.hpp>

#define RUN_BW(N, MC, mstr)                                                    \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##mstr##_N##N;                                 \
        constfilt::Butterworth<double, N##u, constfilt::MC,                    \
                               constfilt::LowPass, false>                      \
            f(100.0, 1000.0);                                                  \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                          \
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

#define RUN_BW_VS_SOS(N, MC, mstr)                                             \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_##mstr##_N##N;                                 \
        constfilt::Butterworth<double, N##u, constfilt::MC,                    \
                               constfilt::LowPass, false>                      \
            f(100.0, 1000.0);                                                  \
        auto r = check(f, Ref::b, Ref::a, Ref::step_sos);                      \
        csv_row("constfilt", "butterworth", N, #mstr "_vs_sos_ref", r);        \
        human_row("constfilt", "butterworth", N, #mstr "_vs_sos_ref", r);      \
    } while (0)

#define RUN_EL(N, MC, mstr)                                                    \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##mstr##_N##N;                                 \
        constfilt::Elliptic<double, N##u, constfilt::MC, constfilt::LowPass,   \
                            false>                                             \
            f(100.0, 0.5, 40.0, 1000.0);                                       \
        auto r = check(f, Ref::b, Ref::a, Ref::step);                          \
        csv_row("constfilt", "elliptic", N, #mstr, r);                         \
        human_row("constfilt", "elliptic", N, #mstr, r);                       \
    } while (0)

#define RUN_EL_SOS(N, MC, mstr)                                                \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##mstr##_N##N;                                 \
        constfilt::Elliptic<double, N##u, constfilt::MC, constfilt::LowPass,   \
                            true>                                              \
            f(100.0, 0.5, 40.0, 1000.0);                                       \
        auto r = check_step([&] { return f(1.0); }, Ref::step_sos);            \
        csv_row("constfilt", "elliptic", N, #mstr "_sos", r);                  \
        human_row("constfilt", "elliptic", N, #mstr "_sos", r);                \
    } while (0)

#define RUN_EL_VS_SOS(N, MC, mstr)                                             \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::el_##mstr##_N##N;                                 \
        constfilt::Elliptic<double, N##u, constfilt::MC, constfilt::LowPass,   \
                            false>                                             \
            f(100.0, 0.5, 40.0, 1000.0);                                       \
        auto r = check(f, Ref::b, Ref::a, Ref::step_sos);                      \
        csv_row("constfilt", "elliptic", N, #mstr "_vs_sos_ref", r);           \
        human_row("constfilt", "elliptic", N, #mstr "_vs_sos_ref", r);         \
    } while (0)

int main()
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

    std::fputs("  [constfilt  Butterworth  MatchedZ  direct vs SOS ref]\n",
               stderr);
    RUN_BW_VS_SOS(1, MatchedZ, matchedz);
    RUN_BW_VS_SOS(2, MatchedZ, matchedz);
    RUN_BW_VS_SOS(3, MatchedZ, matchedz);
    RUN_BW_VS_SOS(4, MatchedZ, matchedz);
    RUN_BW_VS_SOS(5, MatchedZ, matchedz);
    RUN_BW_VS_SOS(6, MatchedZ, matchedz);
    RUN_BW_VS_SOS(7, MatchedZ, matchedz);
    RUN_BW_VS_SOS(8, MatchedZ, matchedz);
    RUN_BW_VS_SOS(9, MatchedZ, matchedz);
    RUN_BW_VS_SOS(10, MatchedZ, matchedz);
    RUN_BW_VS_SOS(11, MatchedZ, matchedz);
    RUN_BW_VS_SOS(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Butterworth  Tustin  direct vs SOS ref]\n",
               stderr);
    RUN_BW_VS_SOS(1, TustinPW, prewarp);
    RUN_BW_VS_SOS(2, TustinPW, prewarp);
    RUN_BW_VS_SOS(3, TustinPW, prewarp);
    RUN_BW_VS_SOS(4, TustinPW, prewarp);
    RUN_BW_VS_SOS(5, TustinPW, prewarp);
    RUN_BW_VS_SOS(6, TustinPW, prewarp);
    RUN_BW_VS_SOS(7, TustinPW, prewarp);
    RUN_BW_VS_SOS(8, TustinPW, prewarp);
    RUN_BW_VS_SOS(9, TustinPW, prewarp);
    RUN_BW_VS_SOS(10, TustinPW, prewarp);
    RUN_BW_VS_SOS(11, TustinPW, prewarp);
    RUN_BW_VS_SOS(12, TustinPW, prewarp);

    std::fputs("  [constfilt  Elliptic  ZOH  direct]\n", stderr);
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

    std::fputs("  [constfilt  Elliptic  MatchedZ  direct]\n", stderr);
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

    std::fputs("  [constfilt  Elliptic  Tustin  direct]\n", stderr);
    RUN_EL(2, TustinPW, prewarp);
    RUN_EL(3, TustinPW, prewarp);
    RUN_EL(4, TustinPW, prewarp);
    RUN_EL(5, TustinPW, prewarp);
    RUN_EL(6, TustinPW, prewarp);
    RUN_EL(7, TustinPW, prewarp);
    RUN_EL(8, TustinPW, prewarp);
    RUN_EL(9, TustinPW, prewarp);
    RUN_EL(10, TustinPW, prewarp);
    RUN_EL(11, TustinPW, prewarp);
    RUN_EL(12, TustinPW, prewarp);

    std::fputs("  [constfilt  Elliptic  MatchedZ  SOS]\n", stderr);
    RUN_EL_SOS(2, MatchedZ, matchedz);
    RUN_EL_SOS(3, MatchedZ, matchedz);
    RUN_EL_SOS(4, MatchedZ, matchedz);
    RUN_EL_SOS(5, MatchedZ, matchedz);
    RUN_EL_SOS(6, MatchedZ, matchedz);
    RUN_EL_SOS(7, MatchedZ, matchedz);
    RUN_EL_SOS(8, MatchedZ, matchedz);
    RUN_EL_SOS(9, MatchedZ, matchedz);
    RUN_EL_SOS(10, MatchedZ, matchedz);
    RUN_EL_SOS(11, MatchedZ, matchedz);
    RUN_EL_SOS(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Elliptic  Tustin  SOS]\n", stderr);
    RUN_EL_SOS(2, TustinPW, prewarp);
    RUN_EL_SOS(3, TustinPW, prewarp);
    RUN_EL_SOS(4, TustinPW, prewarp);
    RUN_EL_SOS(5, TustinPW, prewarp);
    RUN_EL_SOS(6, TustinPW, prewarp);
    RUN_EL_SOS(7, TustinPW, prewarp);
    RUN_EL_SOS(8, TustinPW, prewarp);
    RUN_EL_SOS(9, TustinPW, prewarp);
    RUN_EL_SOS(10, TustinPW, prewarp);
    RUN_EL_SOS(11, TustinPW, prewarp);
    RUN_EL_SOS(12, TustinPW, prewarp);

    std::fputs("  [constfilt  Elliptic  ZOH  direct vs SOS ref]\n", stderr);
    RUN_EL_VS_SOS(2, ZOH, zoh);
    RUN_EL_VS_SOS(3, ZOH, zoh);
    RUN_EL_VS_SOS(4, ZOH, zoh);
    RUN_EL_VS_SOS(5, ZOH, zoh);
    RUN_EL_VS_SOS(6, ZOH, zoh);
    RUN_EL_VS_SOS(7, ZOH, zoh);
    RUN_EL_VS_SOS(8, ZOH, zoh);
    RUN_EL_VS_SOS(9, ZOH, zoh);
    RUN_EL_VS_SOS(10, ZOH, zoh);
    RUN_EL_VS_SOS(11, ZOH, zoh);
    RUN_EL_VS_SOS(12, ZOH, zoh);

    std::fputs("  [constfilt  Elliptic  MatchedZ  direct vs SOS ref]\n",
               stderr);
    RUN_EL_VS_SOS(2, MatchedZ, matchedz);
    RUN_EL_VS_SOS(3, MatchedZ, matchedz);
    RUN_EL_VS_SOS(4, MatchedZ, matchedz);
    RUN_EL_VS_SOS(5, MatchedZ, matchedz);
    RUN_EL_VS_SOS(6, MatchedZ, matchedz);
    RUN_EL_VS_SOS(7, MatchedZ, matchedz);
    RUN_EL_VS_SOS(8, MatchedZ, matchedz);
    RUN_EL_VS_SOS(9, MatchedZ, matchedz);
    RUN_EL_VS_SOS(10, MatchedZ, matchedz);
    RUN_EL_VS_SOS(11, MatchedZ, matchedz);
    RUN_EL_VS_SOS(12, MatchedZ, matchedz);

    std::fputs("  [constfilt  Elliptic  Tustin  direct vs SOS ref]\n", stderr);
    RUN_EL_VS_SOS(2, TustinPW, prewarp);
    RUN_EL_VS_SOS(3, TustinPW, prewarp);
    RUN_EL_VS_SOS(4, TustinPW, prewarp);
    RUN_EL_VS_SOS(5, TustinPW, prewarp);
    RUN_EL_VS_SOS(6, TustinPW, prewarp);
    RUN_EL_VS_SOS(7, TustinPW, prewarp);
    RUN_EL_VS_SOS(8, TustinPW, prewarp);
    RUN_EL_VS_SOS(9, TustinPW, prewarp);
    RUN_EL_VS_SOS(10, TustinPW, prewarp);
    RUN_EL_VS_SOS(11, TustinPW, prewarp);
    RUN_EL_VS_SOS(12, TustinPW, prewarp);

    return 0;
}
