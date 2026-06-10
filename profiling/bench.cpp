// bench.cpp — runtime throughput benchmark (constfilt + iir1)
//
// Compiled as C++17. KFR benchmarks live in bench_kfr.cpp (C++20) to avoid
// mixing KFR headers (C++20) with consteig vendor headers (C++17 aggregate
// initialisation that breaks under C++20 rules).
//
// Build defines (set by CMake when libraries are available via FetchContent):
//   CONSTFILT_BENCH_IIR1  -- enable iir1 section
//
// CSV to stdout:  library,filter_type,order,method,ns_per_sample,msa_per_s,dc_gain
// Human table to stderr.

#include <chrono>
#include <cstdio>

#include <constfilt/constfilt.hpp>

#ifdef CONSTFILT_BENCH_IIR1
#include <Iir.h>
#endif

// ---------------------------------------------------------------------------
// Anti-DCE: all filtered output accumulates here; the volatile write forces
// the compiler to treat every loop iteration as observable.
// ---------------------------------------------------------------------------
static volatile double g_sink = 0.0;

// ---------------------------------------------------------------------------
// Timing configuration
// ---------------------------------------------------------------------------
static constexpr int kReps       = 5;
static constexpr int kSamples    = 5'000'000;
static constexpr int kDcSamples  = 100'000;

using Clock = std::chrono::steady_clock;
using Ns    = std::chrono::nanoseconds;

struct BenchResult
{
    double ns_per_sample;
    double dc_gain;
};

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// constfilt: operator()(double) per-sample path
// F requires: double operator()(double), void reset()
// ---------------------------------------------------------------------------
template <typename F>
static BenchResult bench_single(F &f)
{
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        double acc = 0.0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kSamples; ++i)
            acc += f(1.0);
        const auto t1 = Clock::now();
        g_sink        = acc;
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

#ifdef CONSTFILT_BENCH_IIR1
// ---------------------------------------------------------------------------
// iir1: filter(double) per-sample path
// F requires: double filter(double), void reset()
// ---------------------------------------------------------------------------
template <typename F>
static BenchResult bench_iir1(F &f)
{
    Ns best = Ns::max();
    for (int r = 0; r < kReps; ++r)
    {
        f.reset();
        double acc = 0.0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kSamples; ++i)
            acc += f.filter(1.0);
        const auto t1 = Clock::now();
        g_sink        = acc;
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
#endif // CONSTFILT_BENCH_IIR1


// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main()
{
    // CSV header
    std::puts("library,filter_type,order,method,ns_per_sample,msa_per_s,dc_gain");

    // Human header
    std::fputs("constfilt runtime benchmark  (best of 5 reps x 5M samples)\n\n",
               stderr);
    std::fprintf(stderr, "  %-10s %-12s %-5s %-16s  %13s  %11s  %s\n",
                 "library", "filter_type", "order", "method", "ns/sample",
                 "MSa/s", "dc_gain");
    std::fputs(
        "  -----------------------------------------------------------------------"
        "---------------\n",
        stderr);

    // ── constfilt Butterworth, Tustin ─────────────────────────────────────────
    std::fputs("  [constfilt  Butterworth  Tustin  single-sample]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 1, "tustin", r);
        human_row("constfilt", "butterworth", 1, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 2, "tustin", r);
        human_row("constfilt", "butterworth", 2, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 4, "tustin", r);
        human_row("constfilt", "butterworth", 4, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 6, "tustin", r);
        human_row("constfilt", "butterworth", 6, "tustin", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::Tustin> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 8, "tustin", r);
        human_row("constfilt", "butterworth", 8, "tustin", r);
    }

    // ── constfilt Butterworth, ZOH ────────────────────────────────────────────
    std::fputs("  [constfilt  Butterworth  ZOH  single-sample]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 1, "zoh", r);
        human_row("constfilt", "butterworth", 1, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 2, "zoh", r);
        human_row("constfilt", "butterworth", 2, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 4, "zoh", r);
        human_row("constfilt", "butterworth", 4, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 6, "zoh", r);
        human_row("constfilt", "butterworth", 6, "zoh", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::ZOH> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 8, "zoh", r);
        human_row("constfilt", "butterworth", 8, "zoh", r);
    }

    // ── constfilt Butterworth, MatchedZ ───────────────────────────────────────
    std::fputs("  [constfilt  Butterworth  MatchedZ  single-sample]\n", stderr);
    {
        constfilt::Butterworth<double, 1u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 1, "matchedz", r);
        human_row("constfilt", "butterworth", 1, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 2u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 2, "matchedz", r);
        human_row("constfilt", "butterworth", 2, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 4u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 4, "matchedz", r);
        human_row("constfilt", "butterworth", 4, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 6u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 6, "matchedz", r);
        human_row("constfilt", "butterworth", 6, "matchedz", r);
    }
    {
        constfilt::Butterworth<double, 8u, constfilt::MatchedZ> f(100.0, 1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "butterworth", 8, "matchedz", r);
        human_row("constfilt", "butterworth", 8, "matchedz", r);
    }

    // ── constfilt Elliptic, Tustin ────────────────────────────────────────────
    std::fputs("  [constfilt  Elliptic  Tustin  single-sample]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 2, "tustin", r);
        human_row("constfilt", "elliptic", 2, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 4, "tustin", r);
        human_row("constfilt", "elliptic", 4, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 6, "tustin", r);
        human_row("constfilt", "elliptic", 6, "tustin", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::Tustin> f(100.0, 0.5, 40.0,
                                                             1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 8, "tustin", r);
        human_row("constfilt", "elliptic", 8, "tustin", r);
    }

    // ── constfilt Elliptic, ZOH ───────────────────────────────────────────────
    std::fputs("  [constfilt  Elliptic  ZOH  single-sample]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 2, "zoh", r);
        human_row("constfilt", "elliptic", 2, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 4, "zoh", r);
        human_row("constfilt", "elliptic", 4, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 6, "zoh", r);
        human_row("constfilt", "elliptic", 6, "zoh", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::ZOH> f(100.0, 0.5, 40.0,
                                                          1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 8, "zoh", r);
        human_row("constfilt", "elliptic", 8, "zoh", r);
    }

    // ── constfilt Elliptic, MatchedZ ──────────────────────────────────────────
    std::fputs("  [constfilt  Elliptic  MatchedZ  single-sample]\n", stderr);
    {
        constfilt::Elliptic<double, 2u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 2, "matchedz", r);
        human_row("constfilt", "elliptic", 2, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 4u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 4, "matchedz", r);
        human_row("constfilt", "elliptic", 4, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 6u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 6, "matchedz", r);
        human_row("constfilt", "elliptic", 6, "matchedz", r);
    }
    {
        constfilt::Elliptic<double, 8u, constfilt::MatchedZ> f(100.0, 0.5, 40.0,
                                                               1000.0);
        auto r = bench_single(f);
        csv_row("constfilt", "elliptic", 8, "matchedz", r);
        human_row("constfilt", "elliptic", 8, "matchedz", r);
    }

#ifdef CONSTFILT_BENCH_IIR1
    // ── iir1 Butterworth, per-sample ──────────────────────────────────────────
    std::fputs("  [iir1  Butterworth  runtime-design  single-sample]\n", stderr);
    {
        Iir::Butterworth::LowPass<2> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row("iir1", "butterworth", 2, "runtime", r);
        human_row("iir1", "butterworth", 2, "runtime", r);
    }
    {
        Iir::Butterworth::LowPass<4> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row("iir1", "butterworth", 4, "runtime", r);
        human_row("iir1", "butterworth", 4, "runtime", r);
    }
    {
        Iir::Butterworth::LowPass<6> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row("iir1", "butterworth", 6, "runtime", r);
        human_row("iir1", "butterworth", 6, "runtime", r);
    }
    {
        Iir::Butterworth::LowPass<8> f;
        f.setup(1000.0, 100.0);
        auto r = bench_iir1(f);
        csv_row("iir1", "butterworth", 8, "runtime", r);
        human_row("iir1", "butterworth", 8, "runtime", r);
    }
#endif // CONSTFILT_BENCH_IIR1

    return 0;
}
