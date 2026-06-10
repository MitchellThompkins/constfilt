// bench_constfilt.cpp: constfilt runtime throughput benchmark
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include <chrono>
#include <cstdio>

#include <constfilt/constfilt.hpp>

static volatile double g_sink = 0.0;

static constexpr int kReps      = 5;
static constexpr int kSamples   = 5'000'000;
static constexpr int kDcSamples = 100'000;

using Clock = std::chrono::steady_clock;
using Ns    = std::chrono::nanoseconds;

struct BenchResult
{
    double ns_per_sample;
    double dc_gain;
};

static void csv_row(const char *ftype, int order, const char *method,
                    const BenchResult &r)
{
    std::printf("constfilt,%s,%d,%s,%.3f,%.1f,%.8f\n", ftype, order, method,
                r.ns_per_sample, 1000.0 / r.ns_per_sample, r.dc_gain);
}

static void human_row(const char *ftype, int order, const char *method,
                      const BenchResult &r)
{
    std::fprintf(stderr, "  %-10s %-12s N=%-2d %-16s  %8.2f ns/smpl"
                         "  %7.1f MSa/s  dc=%.6f\n",
                 "constfilt", ftype, order, method, r.ns_per_sample,
                 1000.0 / r.ns_per_sample, r.dc_gain);
}

template <typename F>
static BenchResult bench_single(F &f)
{
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        double acc       = 0.0;
        const auto t0    = Clock::now();
        for (int i = 0; i < kSamples; ++i)
            acc += f(1.0);
        const auto t1    = Clock::now();
        g_sink           = acc;
        const Ns elapsed = t1 - t0;
        if (elapsed < best)
            best = elapsed;
    }
    f.reset();
    double y = 0.0;
    for (int i = 0; i < kDcSamples; ++i)
        y = f(1.0);
    return {static_cast<double>(best.count()) / static_cast<double>(kSamples),
            y};
}

int main()
{
    std::fputs("  [constfilt  Butterworth  Tustin]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 1, "tustin", r);
        human_row("butterworth", 1, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 2, "tustin", r);
        human_row("butterworth", 2, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 4, "tustin", r);
        human_row("butterworth", 4, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 6, "tustin", r);
        human_row("butterworth", 6, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 8, "tustin", r);
        human_row("butterworth", 8, "tustin", r);
    }

    std::fputs("  [constfilt  Butterworth  ZOH]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 1, "zoh", r);
        human_row("butterworth", 1, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 2, "zoh", r);
        human_row("butterworth", 2, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 4, "zoh", r);
        human_row("butterworth", 4, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 6, "zoh", r);
        human_row("butterworth", 6, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 8, "zoh", r);
        human_row("butterworth", 8, "zoh", r);
    }

    std::fputs("  [constfilt  Butterworth  MatchedZ]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 1, "matchedz", r);
        human_row("butterworth", 1, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 2, "matchedz", r);
        human_row("butterworth", 2, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 4, "matchedz", r);
        human_row("butterworth", 4, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 6, "matchedz", r);
        human_row("butterworth", 6, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("butterworth", 8, "matchedz", r);
        human_row("butterworth", 8, "matchedz", r);
    }

    std::fputs("  [constfilt  Elliptic  Tustin]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 2, "tustin", r);
        human_row("elliptic", 2, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 4, "tustin", r);
        human_row("elliptic", 4, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 6, "tustin", r);
        human_row("elliptic", 6, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 8, "tustin", r);
        human_row("elliptic", 8, "tustin", r);
    }

    std::fputs("  [constfilt  Elliptic  ZOH]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 2, "zoh", r);
        human_row("elliptic", 2, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 4, "zoh", r);
        human_row("elliptic", 4, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 6, "zoh", r);
        human_row("elliptic", 6, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 8, "zoh", r);
        human_row("elliptic", 8, "zoh", r);
    }

    std::fputs("  [constfilt  Elliptic  MatchedZ]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 2, "matchedz", r);
        human_row("elliptic", 2, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 4, "matchedz", r);
        human_row("elliptic", 4, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 6, "matchedz", r);
        human_row("elliptic", 6, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("elliptic", 8, "matchedz", r);
        human_row("elliptic", 8, "matchedz", r);
    }

    return 0;
}
