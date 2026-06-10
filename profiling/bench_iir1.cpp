// bench_iir1.cpp: iir1 runtime throughput benchmark
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include <chrono>
#include <cstdio>

#include <Iir.h>

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

static void csv_row(int order, const BenchResult &r)
{
    std::printf("iir1,butterworth,%d,runtime,%.3f,%.1f,%.8f\n", order,
                r.ns_per_sample, 1000.0 / r.ns_per_sample, r.dc_gain);
}

static void human_row(int order, const BenchResult &r)
{
    std::fprintf(stderr, "  %-10s %-12s N=%-2d %-16s  %8.2f ns/smpl"
                         "  %7.1f MSa/s  dc=%.6f\n",
                 "iir1", "butterworth", order, "runtime", r.ns_per_sample,
                 1000.0 / r.ns_per_sample, r.dc_gain);
}

template <typename F>
static BenchResult bench_iir1(F &f)
{
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        double acc       = 0.0;
        const auto t0    = Clock::now();
        for (int i = 0; i < kSamples; ++i)
            acc += f.filter(1.0);
        const auto t1    = Clock::now();
        g_sink           = acc;
        const Ns elapsed = t1 - t0;
        if (elapsed < best)
            best = elapsed;
    }
    f.reset();
    double y = 0.0;
    for (int i = 0; i < kDcSamples; ++i)
        y = f.filter(1.0);
    return {static_cast<double>(best.count()) / static_cast<double>(kSamples),
            y};
}

int main()
{
    std::fputs("  [iir1  Butterworth  runtime-design  single-sample]\n", stderr);
    {
        Iir::Butterworth::LowPass<2> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row(2, r);
        human_row(2, r);
    }
    {
        Iir::Butterworth::LowPass<4> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row(4, r);
        human_row(4, r);
    }
    {
        Iir::Butterworth::LowPass<6> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row(6, r);
        human_row(6, r);
    }
    {
        Iir::Butterworth::LowPass<8> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row(8, r);
        human_row(8, r);
    }

    return 0;
}
