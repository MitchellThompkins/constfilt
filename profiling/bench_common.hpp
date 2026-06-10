#pragma once

// Shared timing helpers for bench_constfilt.cpp, bench_iir1.cpp, bench_kfr.cpp.

#include <chrono>
#include <cstdio>

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

static void csv_row(const char *lib, const char *ftype, int order,
                    const char *method, const BenchResult &r)
{
    std::printf("%s,%s,%d,%s,%.3f,%.1f,%.8f\n", lib, ftype, order, method,
                r.ns_per_sample, 1000.0 / r.ns_per_sample, r.dc_gain);
}

static void human_row(const char *lib, const char *ftype, int order,
                      const char *method, const BenchResult &r)
{
    std::fprintf(stderr, "  %-10s %-12s N=%-2d %-16s  %8.2f ns/smpl"
                         "  %7.1f MSa/s  dc=%.6f\n",
                 lib, ftype, order, method, r.ns_per_sample,
                 1000.0 / r.ns_per_sample, r.dc_gain);
}

// Time a per-sample filter. CallFn(f) must return a double sample.
// F must provide reset().
template <typename F, typename CallFn>
static BenchResult bench(F &f, CallFn call)
{
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        double acc       = 0.0;
        const auto t0    = Clock::now();
        for (int i = 0; i < kSamples; ++i)
            acc += call(f);
        const auto t1    = Clock::now();
        g_sink           = acc;
        const Ns elapsed = t1 - t0;
        if (elapsed < best)
            best = elapsed;
    }
    f.reset();
    double y = 0.0;
    for (int i = 0; i < kDcSamples; ++i)
        y = call(f);
    return {static_cast<double>(best.count()) / static_cast<double>(kSamples),
            y};
}
