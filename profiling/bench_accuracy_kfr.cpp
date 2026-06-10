// bench_accuracy_kfr.cpp: KFR accuracy vs Octave reference (C++20)
//
// KFR uses pre-warped bilinear for both Butterworth and Elliptic designs.
// Compared against acc_ref::bw_prewarp_NX and acc_ref::el_prewarp_NX.
// Only step response is compared (KFR uses SOS internally).
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.

#include "bench_accuracy_common.hpp"
#include <kfr/dsp/iir.hpp>
#include <kfr/dsp/iir_design.hpp>

static double kfr_step_err(kfr::iir_filter<double> &f,
                            const double (&ref_step)[kStepLen])
{
    // KFR requires a single contiguous apply() to maintain state.
    double x_buf[kStepLen];
    double y_buf[kStepLen];
    for (int i = 0; i < kStepLen; ++i)
    {
        x_buf[i] = 1.0;
    }
    f.apply(y_buf, x_buf, static_cast<size_t>(kStepLen));

    double max_e = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double e = std::fabs(y_buf[i] - ref_step[i]);
        if (e > max_e)
        {
            max_e = e;
        }
    }
    return max_e;
}

#define RUN_KFR_BW(N)                                                          \
    do                                                                         \
    {                                                                          \
        using Ref   = acc_ref::bw_prewarp_N##N;                               \
        auto params = kfr::to_sos<double>(                                     \
            kfr::iir_lowpass(kfr::butterworth(N), 100.0, 1000.0));            \
        kfr::iir_filter<double> f(params);                                     \
        AccResult r = {-1.0, -1.0, kfr_step_err(f, Ref::step)};              \
        csv_row("kfr", "butterworth", N, "prewarp", r);                       \
        human_row("kfr", "butterworth", N, "prewarp", r);                     \
    } while (0)

#ifdef KFR_HAVE_ELLIPTIC
#define RUN_KFR_EL(N)                                                          \
    do                                                                         \
    {                                                                          \
        using Ref   = acc_ref::el_prewarp_N##N;                               \
        auto params = kfr::to_sos<double>(kfr::iir_lowpass(                   \
            kfr::elliptic(N, 0.5, 40.0), 100.0, 1000.0));                     \
        kfr::iir_filter<double> f(params);                                     \
        AccResult r = {-1.0, -1.0, kfr_step_err(f, Ref::step)};              \
        csv_row("kfr", "elliptic", N, "prewarp", r);                          \
        human_row("kfr", "elliptic", N, "prewarp", r);                        \
    } while (0)
#endif

int main()
{
    std::fputs("  [KFR  Butterworth  prewarp bilinear]\n", stderr);
    RUN_KFR_BW(1);
    RUN_KFR_BW(2);
    RUN_KFR_BW(3);
    RUN_KFR_BW(4);
    RUN_KFR_BW(5);
    RUN_KFR_BW(6);
    RUN_KFR_BW(7);
    RUN_KFR_BW(8);
    RUN_KFR_BW(9);
    RUN_KFR_BW(10);
    RUN_KFR_BW(11);
    RUN_KFR_BW(12);

#ifdef KFR_HAVE_ELLIPTIC
    std::fputs("  [KFR  Elliptic  prewarp bilinear]\n", stderr);
    RUN_KFR_EL(2);
    RUN_KFR_EL(3);
    RUN_KFR_EL(4);
    RUN_KFR_EL(5);
    RUN_KFR_EL(6);
    RUN_KFR_EL(7);
    RUN_KFR_EL(8);
    RUN_KFR_EL(9);
    RUN_KFR_EL(10);
    RUN_KFR_EL(11);
    RUN_KFR_EL(12);
#endif

    return 0;
}
