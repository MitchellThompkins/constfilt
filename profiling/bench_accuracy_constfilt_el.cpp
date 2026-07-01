// bench_accuracy_constfilt_el.cpp: Elliptic accuracy rows.
// Compiled as a separate TU from the Butterworth rows to keep per-TU template
// instantiation count (and GCC memory usage) manageable at N=25.

#include "bench_accuracy_common.hpp"
#include <constfilt/constfilt.hpp>

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

void run_el_accuracy()
{
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
    RUN_EL(13, ZOH, zoh);
    RUN_EL(14, ZOH, zoh);
    RUN_EL(15, ZOH, zoh);
    RUN_EL(16, ZOH, zoh);
    RUN_EL(17, ZOH, zoh);
    RUN_EL(18, ZOH, zoh);
    RUN_EL(19, ZOH, zoh);
    RUN_EL(20, ZOH, zoh);
    RUN_EL(21, ZOH, zoh);
    RUN_EL(22, ZOH, zoh);
    RUN_EL(23, ZOH, zoh);
    RUN_EL(24, ZOH, zoh);
    RUN_EL(25, ZOH, zoh);

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
    RUN_EL(13, MatchedZ, matchedz);
    RUN_EL(14, MatchedZ, matchedz);
    RUN_EL(15, MatchedZ, matchedz);
    RUN_EL(16, MatchedZ, matchedz);
    RUN_EL(17, MatchedZ, matchedz);
    RUN_EL(18, MatchedZ, matchedz);
    RUN_EL(19, MatchedZ, matchedz);
    RUN_EL(20, MatchedZ, matchedz);
    RUN_EL(21, MatchedZ, matchedz);
    RUN_EL(22, MatchedZ, matchedz);
    RUN_EL(23, MatchedZ, matchedz);
    RUN_EL(24, MatchedZ, matchedz);
    RUN_EL(25, MatchedZ, matchedz);

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
    RUN_EL(13, TustinPW, prewarp);
    RUN_EL(14, TustinPW, prewarp);
    RUN_EL(15, TustinPW, prewarp);
    RUN_EL(16, TustinPW, prewarp);
    RUN_EL(17, TustinPW, prewarp);
    RUN_EL(18, TustinPW, prewarp);
    RUN_EL(19, TustinPW, prewarp);
    RUN_EL(20, TustinPW, prewarp);
    RUN_EL(21, TustinPW, prewarp);
    RUN_EL(22, TustinPW, prewarp);
    RUN_EL(23, TustinPW, prewarp);
    RUN_EL(24, TustinPW, prewarp);
    RUN_EL(25, TustinPW, prewarp);

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
    RUN_EL_SOS(13, MatchedZ, matchedz);
    RUN_EL_SOS(14, MatchedZ, matchedz);
    RUN_EL_SOS(15, MatchedZ, matchedz);
    RUN_EL_SOS(16, MatchedZ, matchedz);
    RUN_EL_SOS(17, MatchedZ, matchedz);
    RUN_EL_SOS(18, MatchedZ, matchedz);
    RUN_EL_SOS(19, MatchedZ, matchedz);
    RUN_EL_SOS(20, MatchedZ, matchedz);
    RUN_EL_SOS(21, MatchedZ, matchedz);
    RUN_EL_SOS(22, MatchedZ, matchedz);
    RUN_EL_SOS(23, MatchedZ, matchedz);
    RUN_EL_SOS(24, MatchedZ, matchedz);
    RUN_EL_SOS(25, MatchedZ, matchedz);

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
    RUN_EL_SOS(13, TustinPW, prewarp);
    RUN_EL_SOS(14, TustinPW, prewarp);
    RUN_EL_SOS(15, TustinPW, prewarp);
    RUN_EL_SOS(16, TustinPW, prewarp);
    RUN_EL_SOS(17, TustinPW, prewarp);
    RUN_EL_SOS(18, TustinPW, prewarp);
    RUN_EL_SOS(19, TustinPW, prewarp);
    RUN_EL_SOS(20, TustinPW, prewarp);
    RUN_EL_SOS(21, TustinPW, prewarp);
    RUN_EL_SOS(22, TustinPW, prewarp);
    RUN_EL_SOS(23, TustinPW, prewarp);
    RUN_EL_SOS(24, TustinPW, prewarp);
    RUN_EL_SOS(25, TustinPW, prewarp);

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
    RUN_EL_VS_SOS(13, ZOH, zoh);
    RUN_EL_VS_SOS(14, ZOH, zoh);
    RUN_EL_VS_SOS(15, ZOH, zoh);
    RUN_EL_VS_SOS(16, ZOH, zoh);
    RUN_EL_VS_SOS(17, ZOH, zoh);
    RUN_EL_VS_SOS(18, ZOH, zoh);
    RUN_EL_VS_SOS(19, ZOH, zoh);
    RUN_EL_VS_SOS(20, ZOH, zoh);
    RUN_EL_VS_SOS(21, ZOH, zoh);
    RUN_EL_VS_SOS(22, ZOH, zoh);
    RUN_EL_VS_SOS(23, ZOH, zoh);
    RUN_EL_VS_SOS(24, ZOH, zoh);
    RUN_EL_VS_SOS(25, ZOH, zoh);

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
    RUN_EL_VS_SOS(13, MatchedZ, matchedz);
    RUN_EL_VS_SOS(14, MatchedZ, matchedz);
    RUN_EL_VS_SOS(15, MatchedZ, matchedz);
    RUN_EL_VS_SOS(16, MatchedZ, matchedz);
    RUN_EL_VS_SOS(17, MatchedZ, matchedz);
    RUN_EL_VS_SOS(18, MatchedZ, matchedz);
    RUN_EL_VS_SOS(19, MatchedZ, matchedz);
    RUN_EL_VS_SOS(20, MatchedZ, matchedz);
    RUN_EL_VS_SOS(21, MatchedZ, matchedz);
    RUN_EL_VS_SOS(22, MatchedZ, matchedz);
    RUN_EL_VS_SOS(23, MatchedZ, matchedz);
    RUN_EL_VS_SOS(24, MatchedZ, matchedz);
    RUN_EL_VS_SOS(25, MatchedZ, matchedz);

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
    RUN_EL_VS_SOS(13, TustinPW, prewarp);
    RUN_EL_VS_SOS(14, TustinPW, prewarp);
    RUN_EL_VS_SOS(15, TustinPW, prewarp);
    RUN_EL_VS_SOS(16, TustinPW, prewarp);
    RUN_EL_VS_SOS(17, TustinPW, prewarp);
    RUN_EL_VS_SOS(18, TustinPW, prewarp);
    RUN_EL_VS_SOS(19, TustinPW, prewarp);
    RUN_EL_VS_SOS(20, TustinPW, prewarp);
    RUN_EL_VS_SOS(21, TustinPW, prewarp);
    RUN_EL_VS_SOS(22, TustinPW, prewarp);
    RUN_EL_VS_SOS(23, TustinPW, prewarp);
    RUN_EL_VS_SOS(24, TustinPW, prewarp);
    RUN_EL_VS_SOS(25, TustinPW, prewarp);
}
