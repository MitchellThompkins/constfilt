// bench_constfilt.cpp: constfilt runtime throughput benchmark
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_common.hpp"
#include <constfilt/constfilt.hpp>

#define RUN_CF(ftype, N, METHOD, mstr)                                         \
    do                                                                         \
    {                                                                          \
        constfilt::ftype<double, N##u, constfilt::METHOD, constfilt::LowPass,  \
                         false>                                                \
            f(100.0, 1000.0);                                                  \
        auto r = bench(f, [](auto &f) { return f(1.0); });                     \
        csv_row("constfilt", #ftype, N, mstr, r);                              \
        human_row("constfilt", #ftype, N, mstr, r);                            \
    } while (0)

#define RUN_CF_SOS(ftype, N, METHOD, mstr)                                     \
    do                                                                         \
    {                                                                          \
        constfilt::ftype<double, N##u, constfilt::METHOD, constfilt::LowPass,  \
                         true>                                                 \
            f(100.0, 1000.0);                                                  \
        auto r = bench(f, [](auto &f) { return f(1.0); });                     \
        csv_row("constfilt", #ftype, N, mstr "_sos", r);                       \
        human_row("constfilt", #ftype, N, mstr "_sos", r);                     \
    } while (0)

#define RUN_EL(N, METHOD, mstr)                                                \
    do                                                                         \
    {                                                                          \
        constfilt::Elliptic<double, N##u, constfilt::METHOD,                   \
                            constfilt::LowPass, false>                         \
            f(100.0, 0.5, 40.0, 1000.0);                                       \
        auto r = bench(f, [](auto &f) { return f(1.0); });                     \
        csv_row("constfilt", "Elliptic", N, mstr, r);                          \
        human_row("constfilt", "Elliptic", N, mstr, r);                        \
    } while (0)

#define RUN_EL_SOS(N, METHOD, mstr)                                            \
    do                                                                         \
    {                                                                          \
        constfilt::Elliptic<double, N##u, constfilt::METHOD,                   \
                            constfilt::LowPass, true>                          \
            f(100.0, 0.5, 40.0, 1000.0);                                       \
        auto r = bench(f, [](auto &f) { return f(1.0); });                     \
        csv_row("constfilt", "Elliptic", N, mstr "_sos", r);                   \
        human_row("constfilt", "Elliptic", N, mstr "_sos", r);                 \
    } while (0)

int main()
{
    std::fputs("  [constfilt  Butterworth  TustinPW  direct]\n", stderr);
    RUN_CF(Butterworth, 1, TustinPW, "tustin");
    RUN_CF(Butterworth, 2, TustinPW, "tustin");
    RUN_CF(Butterworth, 4, TustinPW, "tustin");
    RUN_CF(Butterworth, 6, TustinPW, "tustin");
    RUN_CF(Butterworth, 8, TustinPW, "tustin");

    std::fputs("  [constfilt  Butterworth  ZOH  direct]\n", stderr);
    RUN_CF(Butterworth, 1, ZOH, "zoh");
    RUN_CF(Butterworth, 2, ZOH, "zoh");
    RUN_CF(Butterworth, 4, ZOH, "zoh");
    RUN_CF(Butterworth, 6, ZOH, "zoh");
    RUN_CF(Butterworth, 8, ZOH, "zoh");

    std::fputs("  [constfilt  Butterworth  MatchedZ  direct]\n", stderr);
    RUN_CF(Butterworth, 1, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 2, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 4, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 6, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 8, MatchedZ, "matchedz");

    std::fputs("  [constfilt  Butterworth  TustinPW  SOS]\n", stderr);
    RUN_CF_SOS(Butterworth, 1, TustinPW, "tustin");
    RUN_CF_SOS(Butterworth, 2, TustinPW, "tustin");
    RUN_CF_SOS(Butterworth, 4, TustinPW, "tustin");
    RUN_CF_SOS(Butterworth, 6, TustinPW, "tustin");
    RUN_CF_SOS(Butterworth, 8, TustinPW, "tustin");

    std::fputs("  [constfilt  Butterworth  MatchedZ  SOS]\n", stderr);
    RUN_CF_SOS(Butterworth, 1, MatchedZ, "matchedz");
    RUN_CF_SOS(Butterworth, 2, MatchedZ, "matchedz");
    RUN_CF_SOS(Butterworth, 4, MatchedZ, "matchedz");
    RUN_CF_SOS(Butterworth, 6, MatchedZ, "matchedz");
    RUN_CF_SOS(Butterworth, 8, MatchedZ, "matchedz");

    std::fputs("  [constfilt  Elliptic  TustinPW  direct]\n", stderr);
    RUN_EL(2, TustinPW, "tustin");
    RUN_EL(4, TustinPW, "tustin");
    RUN_EL(6, TustinPW, "tustin");
    RUN_EL(8, TustinPW, "tustin");

    std::fputs("  [constfilt  Elliptic  ZOH  direct]\n", stderr);
    RUN_EL(2, ZOH, "zoh");
    RUN_EL(4, ZOH, "zoh");
    RUN_EL(6, ZOH, "zoh");
    RUN_EL(8, ZOH, "zoh");

    std::fputs("  [constfilt  Elliptic  MatchedZ  direct]\n", stderr);
    RUN_EL(2, MatchedZ, "matchedz");
    RUN_EL(4, MatchedZ, "matchedz");
    RUN_EL(6, MatchedZ, "matchedz");
    RUN_EL(8, MatchedZ, "matchedz");

    std::fputs("  [constfilt  Elliptic  TustinPW  SOS]\n", stderr);
    RUN_EL_SOS(2, TustinPW, "tustin");
    RUN_EL_SOS(4, TustinPW, "tustin");
    RUN_EL_SOS(6, TustinPW, "tustin");
    RUN_EL_SOS(8, TustinPW, "tustin");

    std::fputs("  [constfilt  Elliptic  MatchedZ  SOS]\n", stderr);
    RUN_EL_SOS(2, MatchedZ, "matchedz");
    RUN_EL_SOS(4, MatchedZ, "matchedz");
    RUN_EL_SOS(6, MatchedZ, "matchedz");
    RUN_EL_SOS(8, MatchedZ, "matchedz");

    return 0;
}
