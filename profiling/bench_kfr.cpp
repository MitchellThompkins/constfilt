// bench_kfr.cpp — KFR runtime throughput benchmark (C++20)
//
// Compiled as C++20 (required by KFR). Kept separate from bench.cpp (C++17)
// because consteig vendor headers use C++17 aggregate initialisation that
// breaks under C++20 rules.
//
// KFR is SIMD-accelerated; benchmarks use batch=256 to expose that advantage.
// method column is "runtime+simd" to distinguish from constfilt's compile-time
// scalar path. ns/sample is still per-sample throughput.
//
// Output: CSV rows to stdout (no header — appended after bench output).
//         Human-readable table to stderr.

#include <chrono>
#include <cstdio>

#include <kfr/dsp/iir.hpp>
#include <kfr/dsp/iir_design.hpp>

static volatile double g_sink = 0.0;

static constexpr int kReps      = 5;
static constexpr int kSamples   = 5'000'000;
static constexpr int kDcSamples = 100'000;
static constexpr int kBatch     = 256;

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

static BenchResult bench_kfr(kfr::iir_filter<double> &f)
{
    double x_buf[kBatch];
    double y_buf[kBatch];
    for (int i = 0; i < kBatch; ++i)
        x_buf[i] = 1.0;

    const int n_batches = kSamples / kBatch;
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        const auto t0 = Clock::now();
        for (int b = 0; b < n_batches; ++b)
            f.apply(y_buf, x_buf, static_cast<size_t>(kBatch));
        const auto t1 = Clock::now();
        g_sink          = y_buf[0];
        const Ns elapsed = t1 - t0;
        if (elapsed < best)
            best = elapsed;
    }
    f.reset();
    for (int i = 0; i < kDcSamples / kBatch; ++i)
        f.apply(y_buf, x_buf, static_cast<size_t>(kBatch));
    const double y = y_buf[kBatch - 1];

    const long long total =
        static_cast<long long>(n_batches) * static_cast<long long>(kBatch);
    return {static_cast<double>(best.count()) / static_cast<double>(total), y};
}

int main()
{
    std::fprintf(stderr,
                 "  [KFR  Butterworth+Elliptic  runtime+simd  batch=%d]\n",
                 kBatch);

    // ── KFR Butterworth ───────────────────────────────────────────────────────
    for (int order : {2, 4, 6, 8})
    {
        auto params =
            kfr::to_sos<double>(kfr::iir_lowpass(kfr::butterworth(order), 100.0, 1000.0));
        kfr::iir_filter<double> f(params);
        auto r = bench_kfr(f);
        csv_row("kfr", "butterworth", order, "runtime+simd", r);
        human_row("kfr", "butterworth", order, "runtime+simd", r);
    }

#ifdef KFR_HAVE_ELLIPTIC
    // ── KFR Elliptic ──────────────────────────────────────────────────────────
    for (int order : {2, 4, 6, 8})
    {
        auto params = kfr::to_sos<double>(
            kfr::iir_lowpass(kfr::elliptic(order, 0.5, 40.0), 100.0, 1000.0));
        kfr::iir_filter<double> f(params);
        auto r = bench_kfr(f);
        csv_row("kfr", "elliptic", order, "runtime+simd", r);
        human_row("kfr", "elliptic", order, "runtime+simd", r);
    }
#endif

    return 0;
}
