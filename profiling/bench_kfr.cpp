// bench_kfr.cpp: KFR runtime throughput benchmark (C++20)
//
// KFR is designed for SIMD batch throughput; a batch of 256 samples exposes
// that advantage. Results are labelled "runtime+simd,batch=256".
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_common.hpp"
#include <kfr/dsp/iir.hpp>
#include <kfr/dsp/iir_design.hpp>

static constexpr int kBatch = 256;

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
        const auto t0    = Clock::now();
        for (int b = 0; b < n_batches; ++b)
            f.apply(y_buf, x_buf, static_cast<size_t>(kBatch));
        const auto t1    = Clock::now();
        g_sink           = y_buf[0];
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

#define RUN_KFR_BW(N)                                                          \
    do                                                                         \
    {                                                                          \
        auto params = kfr::to_sos<double>(                                     \
            kfr::iir_lowpass(kfr::butterworth(N), 100.0, 1000.0));            \
        kfr::iir_filter<double> f(params);                                     \
        auto r = bench_kfr(f);                                                 \
        csv_row("kfr", "Butterworth", N, "runtime+simd", r);                  \
        human_row("kfr", "Butterworth", N, "runtime+simd", r);                \
    } while (0)

#define RUN_KFR_EL(N)                                                          \
    do                                                                         \
    {                                                                          \
        auto params = kfr::to_sos<double>(                                     \
            kfr::iir_lowpass(kfr::elliptic(N, 0.5, 40.0), 100.0, 1000.0));   \
        kfr::iir_filter<double> f(params);                                     \
        auto r = bench_kfr(f);                                                 \
        csv_row("kfr", "Elliptic", N, "runtime+simd", r);                     \
        human_row("kfr", "Elliptic", N, "runtime+simd", r);                   \
    } while (0)

int main()
{
    std::fprintf(stderr, "  [KFR  Butterworth  runtime+simd  batch=%d]\n",
                 kBatch);
    RUN_KFR_BW(2);
    RUN_KFR_BW(4);
    RUN_KFR_BW(6);
    RUN_KFR_BW(8);

#ifdef KFR_HAVE_ELLIPTIC
    std::fprintf(stderr, "  [KFR  Elliptic  runtime+simd  batch=%d]\n",
                 kBatch);
    RUN_KFR_EL(2);
    RUN_KFR_EL(4);
    RUN_KFR_EL(6);
    RUN_KFR_EL(8);
#endif

    return 0;
}
