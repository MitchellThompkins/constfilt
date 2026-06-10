// bench_accuracy_iir1.cpp: iir1 accuracy vs Octave pre-warped Tustin reference
//
// iir1 uses pre-warped bilinear (exact digital cutoff). Compared against
// acc_ref::bw_prewarp_NX. Only step response is compared (iir1 uses SOS
// internally and does not expose direct-form b/a polynomials).
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include <cmath>
#include <cstdio>

#include "accuracy_reference.hpp"
#include <Iir.h>

static constexpr int kStepLen = 256;

static void csv_row(int order, double step_err)
{
    std::printf("iir1,butterworth,%d,prewarp,n/a,n/a,%.6e\n", order, step_err);
}

static void human_row(int order, double step_err)
{
    const char *flag = step_err > 1e-6 ? " !" : "  ";
    std::fprintf(stderr,
                 "%s %-10s %-12s N=%-2d %-10s  b=n/a       a=n/a       "
                 "step=%.2e\n",
                 flag, "iir1", "butterworth", order, "prewarp", step_err);
}

int main()
{
    std::fputs("  [iir1  Butterworth  prewarp bilinear]\n", stderr);

#define IIR1_BW(N)                                                             \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_prewarp_N##N;                                 \
        Iir::Butterworth::LowPass<N> f;                                        \
        f.setup(1000.0, 100.0);                                                \
        double max_s = 0.0;                                                    \
        for (int i = 0; i < kStepLen; ++i)                                    \
        {                                                                      \
            double e = std::fabs(f.filter(1.0) - Ref::step[i]);               \
            if (e > max_s)                                                     \
                max_s = e;                                                     \
        }                                                                      \
        csv_row(N, max_s);                                                     \
        human_row(N, max_s);                                                   \
    } while (0)

    IIR1_BW(1);
    IIR1_BW(2);
    IIR1_BW(3);
    IIR1_BW(4);
    IIR1_BW(5);
    IIR1_BW(6);
    IIR1_BW(7);
    IIR1_BW(8);
    IIR1_BW(9);
    IIR1_BW(10);
    IIR1_BW(11);
    IIR1_BW(12);

#undef IIR1_BW

    return 0;
}
