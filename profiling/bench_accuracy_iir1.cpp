// bench_accuracy_iir1.cpp: iir1 accuracy vs Octave pre-warped Tustin reference
//
// iir1 uses pre-warped bilinear (exact digital cutoff). Compared against
// acc_ref::bw_prewarp_NX. Only step response is compared (iir1 uses SOS
// internally and does not expose direct-form b/a polynomials).
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_accuracy_common.hpp"
#include <Iir.h>

#define RUN_IIR1(N)                                                            \
    do                                                                         \
    {                                                                          \
        using Ref = acc_ref::bw_prewarp_N##N;                                  \
        Iir::Butterworth::LowPass<N> f;                                        \
        f.setup(1000.0, 100.0);                                                \
        auto r = check_step([&]() { return f.filter(1.0); }, Ref::step_sos);   \
        csv_row("iir1", "butterworth", N, "prewarp", r);                       \
        human_row("iir1", "butterworth", N, "prewarp", r);                     \
    } while (0)

int main()
{
    std::fputs("  [iir1  Butterworth  prewarp bilinear]\n", stderr);
    RUN_IIR1(1);
    RUN_IIR1(2);
    RUN_IIR1(3);
    RUN_IIR1(4);
    RUN_IIR1(5);
    RUN_IIR1(6);
    RUN_IIR1(7);
    RUN_IIR1(8);
    RUN_IIR1(9);
    RUN_IIR1(10);
    RUN_IIR1(11);
    RUN_IIR1(12);
    RUN_IIR1(13);
    RUN_IIR1(14);
    RUN_IIR1(15);
    RUN_IIR1(16);
    RUN_IIR1(17);
    RUN_IIR1(18);
    RUN_IIR1(19);
    RUN_IIR1(20);
    RUN_IIR1(21);
    RUN_IIR1(22);
    RUN_IIR1(23);
    RUN_IIR1(24);
    RUN_IIR1(25);

    return 0;
}
