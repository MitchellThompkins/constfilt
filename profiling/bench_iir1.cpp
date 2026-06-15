// bench_iir1.cpp: iir1 runtime throughput benchmark
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_common.hpp"
#include <Iir.h>

#define RUN_IIR1(N)                                                            \
    do                                                                         \
    {                                                                          \
        Iir::Butterworth::LowPass<N> f;                                        \
        f.setup(1000.0, 100.0);                                                \
        auto r = bench(f, [](auto &f) { return f.filter(1.0); });              \
        csv_row("iir1", "Butterworth", N, "runtime", r);                       \
        human_row("iir1", "Butterworth", N, "runtime", r);                     \
    } while (0)

int main()
{
    std::fputs("  [iir1  Butterworth  runtime-design  single-sample]\n",
               stderr);
    RUN_IIR1(2);
    RUN_IIR1(4);
    RUN_IIR1(6);
    RUN_IIR1(8);

    return 0;
}
