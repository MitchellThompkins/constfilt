// bench_constfilt.cpp: constfilt runtime throughput benchmark
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_common.hpp"
#include <constfilt/constfilt.hpp>

#define RUN_CF(ftype, N, METHOD, mstr)                                         \
    do                                                                         \
    {                                                                          \
        constfilt::ftype<double, N##u, constfilt::METHOD> f(100.0, 1000.0);   \
        auto r = bench(f, [](auto &f) { return f(1.0); });                    \
        csv_row("constfilt", #ftype, N, mstr, r);                             \
        human_row("constfilt", #ftype, N, mstr, r);                           \
    } while (0)

#define RUN_EL(N, METHOD, mstr)                                                \
    do                                                                         \
    {                                                                          \
        constfilt::Elliptic<double, N##u, constfilt::METHOD> f(100.0, 0.5,    \
                                                               40.0, 1000.0); \
        auto r = bench(f, [](auto &f) { return f(1.0); });                    \
        csv_row("constfilt", "Elliptic", N, mstr, r);                         \
        human_row("constfilt", "Elliptic", N, mstr, r);                       \
    } while (0)

int main()
{
    std::fputs("  [constfilt  Butterworth  Tustin]\n", stderr);
    RUN_CF(Butterworth, 1, Tustin, "tustin");
    RUN_CF(Butterworth, 2, Tustin, "tustin");
    RUN_CF(Butterworth, 4, Tustin, "tustin");
    RUN_CF(Butterworth, 6, Tustin, "tustin");
    RUN_CF(Butterworth, 8, Tustin, "tustin");

    std::fputs("  [constfilt  Butterworth  ZOH]\n", stderr);
    RUN_CF(Butterworth, 1, ZOH, "zoh");
    RUN_CF(Butterworth, 2, ZOH, "zoh");
    RUN_CF(Butterworth, 4, ZOH, "zoh");
    RUN_CF(Butterworth, 6, ZOH, "zoh");
    RUN_CF(Butterworth, 8, ZOH, "zoh");

    std::fputs("  [constfilt  Butterworth  MatchedZ]\n", stderr);
    RUN_CF(Butterworth, 1, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 2, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 4, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 6, MatchedZ, "matchedz");
    RUN_CF(Butterworth, 8, MatchedZ, "matchedz");

    std::fputs("  [constfilt  Elliptic  Tustin]\n", stderr);
    RUN_EL(2, Tustin, "tustin");
    RUN_EL(4, Tustin, "tustin");
    RUN_EL(6, Tustin, "tustin");
    RUN_EL(8, Tustin, "tustin");

    std::fputs("  [constfilt  Elliptic  ZOH]\n", stderr);
    RUN_EL(2, ZOH, "zoh");
    RUN_EL(4, ZOH, "zoh");
    RUN_EL(6, ZOH, "zoh");
    RUN_EL(8, ZOH, "zoh");

    std::fputs("  [constfilt  Elliptic  MatchedZ]\n", stderr);
    RUN_EL(2, MatchedZ, "matchedz");
    RUN_EL(4, MatchedZ, "matchedz");
    RUN_EL(6, MatchedZ, "matchedz");
    RUN_EL(8, MatchedZ, "matchedz");

    return 0;
}
