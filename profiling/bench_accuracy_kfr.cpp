// bench_accuracy_kfr.cpp — KFR accuracy vs Octave reference (C++20)
//
// KFR uses bilinear/Tustin for both Butterworth and Elliptic designs, so all
// comparisons are against acc_ref::bw_tustin_NX and acc_ref::el_tustin_NX.
//
// KFR uses SOS (biquad cascade) internally and does not expose direct-form
// b/a polynomials; only step response is compared.
//
// CSV rows appended to the accuracy CSV (no header — bench_accuracy writes it).
// Human table to stderr.

#include <cmath>
#include <cstdio>

#include "accuracy_reference.hpp"
#include <kfr/dsp/iir.hpp>
#include <kfr/dsp/iir_design.hpp>

static constexpr int kStepLen = 256;

static void csv_row(const char *lib, const char *ftype, int order,
                    const char *method, double step_err)
{
    std::printf("%s,%s,%d,%s,n/a,n/a,%.6e\n", lib, ftype, order, method,
                step_err);
}

static void human_row(const char *lib, const char *ftype, int order,
                      const char *method, double step_err)
{
    const char *flag = step_err > 1e-6 ? " !" : "  ";
    std::fprintf(stderr,
                 "%s %-10s %-12s N=%-2d %-10s  b=n/a       a=n/a       "
                 "step=%.2e\n",
                 flag, lib, ftype, order, method, step_err);
}

static double max_step_err(kfr::iir_filter<double> &f,
                            const double (&ref_step)[kStepLen])
{
    // KFR's IIR expression engine requires a single contiguous apply() call to
    // maintain state correctly across samples. Calling with size=1 repeatedly
    // does not accumulate state between calls.
    double x_buf[kStepLen];
    double y_buf[kStepLen];
    for (int i = 0; i < kStepLen; ++i)
        x_buf[i] = 1.0;
    f.apply(y_buf, x_buf, static_cast<size_t>(kStepLen));

    double max_e = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double e = std::fabs(y_buf[i] - ref_step[i]);
        if (e > max_e)
            max_e = e;
    }
    return max_e;
}

int main()
{
    // ── KFR Butterworth (Tustin/bilinear) ─────────────────────────────────────
    // KFR uses pre-warped bilinear (exact digital cutoff), same convention as
    // iir1.  Systematic ~N*1% step-response offset vs constfilt/Octave standard
    // bilinear.  Not numerical breakdown -- design convention difference.
    std::fputs("  [KFR  Butterworth  runtime+simd  tustin(bilinear) -- pre-warped]\n",
               stderr);
    for (int order = 1; order <= 12; ++order)
    {
        auto params = kfr::to_sos<double>(
            kfr::iir_lowpass(kfr::butterworth(order), 100.0, 1000.0));
        kfr::iir_filter<double> f(params);

        // Select the correct acc_ref struct by order via a switch.
#define REF_CASE(N)                                                            \
    case N:                                                                    \
        e = max_step_err(f, acc_ref::bw_tustin_N##N::step);                   \
        break

        double e = 0.0;
        switch (order)
        {
            REF_CASE(1);
            REF_CASE(2);
            REF_CASE(3);
            REF_CASE(4);
            REF_CASE(5);
            REF_CASE(6);
            REF_CASE(7);
            REF_CASE(8);
            REF_CASE(9);
            REF_CASE(10);
            REF_CASE(11);
            REF_CASE(12);
        default:
            break;
        }
#undef REF_CASE

        csv_row("kfr", "butterworth", order, "tustin", e);
        human_row("kfr", "butterworth", order, "tustin", e);
    }

#ifdef KFR_HAVE_ELLIPTIC
    // ── KFR Elliptic (Tustin/bilinear) ───────────────────────────────────────
    std::fputs("  [KFR  Elliptic  runtime+simd  tustin(bilinear)]\n", stderr);
    for (int order = 2; order <= 12; ++order)
    {
        auto params = kfr::to_sos<double>(
            kfr::iir_lowpass(kfr::elliptic(order, 0.5, 40.0), 100.0, 1000.0));
        kfr::iir_filter<double> f(params);

#define REF_CASE(N)                                                            \
    case N:                                                                    \
        e = max_step_err(f, acc_ref::el_tustin_N##N::step);                   \
        break

        double e = 0.0;
        switch (order)
        {
            REF_CASE(2);
            REF_CASE(3);
            REF_CASE(4);
            REF_CASE(5);
            REF_CASE(6);
            REF_CASE(7);
            REF_CASE(8);
            REF_CASE(9);
            REF_CASE(10);
            REF_CASE(11);
            REF_CASE(12);
        default:
            break;
        }
#undef REF_CASE

        csv_row("kfr", "elliptic", order, "tustin", e);
        human_row("kfr", "elliptic", order, "elliptic", e);
    }
#endif // KFR_HAVE_ELLIPTIC

    return 0;
}
