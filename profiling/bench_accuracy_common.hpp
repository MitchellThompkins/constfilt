#pragma once

// Shared accuracy helpers for bench_accuracy_constfilt.cpp,
// bench_accuracy_iir1.cpp, bench_accuracy_kfr.cpp.

#include <cmath>
#include <cstdio>

#include "accuracy_reference.hpp"

static constexpr int kStepLen = 256;
static constexpr double kErrThreshold = 1e-6;

struct AccResult
{
    double max_b_err; // -1 when b/a not available (SOS-internal libraries)
    double max_a_err;
    double max_step_err;
};

static void csv_row(const char *lib, const char *ftype, int order,
                    const char *method, const AccResult &r)
{
    if (r.max_b_err < 0.0)
    {
        std::printf("%s,%s,%d,%s,n/a,n/a,%.6e\n", lib, ftype, order, method,
                    r.max_step_err);
    }
    else
    {
        std::printf("%s,%s,%d,%s,%.6e,%.6e,%.6e\n", lib, ftype, order, method,
                    r.max_b_err, r.max_a_err, r.max_step_err);
    }
}

static void human_row(const char *lib, const char *ftype, int order,
                      const char *method, const AccResult &r)
{
    const char *flag = r.max_step_err > kErrThreshold ? " !" : "  ";
    if (r.max_b_err < 0.0)
    {
        std::fprintf(stderr,
                     "%s %-10s %-12s N=%-2d %-10s  b=n/a       a=n/a      "
                     " step=%.2e\n",
                     flag, lib, ftype, order, method, r.max_step_err);
    }
    else
    {
        std::fprintf(stderr,
                     "%s %-10s %-12s N=%-2d %-10s  b=%.2e  a=%.2e  "
                     "step=%.2e\n",
                     flag, lib, ftype, order, method, r.max_b_err, r.max_a_err,
                     r.max_step_err);
    }
}

// Full accuracy check for filters that expose b/a polynomials.
template <size_t NC, typename F>
static AccResult check(F &filt, const double (&ref_b)[NC],
                       const double (&ref_a)[NC],
                       const double (&ref_step)[kStepLen])
{
    double max_b = 0.0;
    for (size_t i = 0; i < NC; ++i)
    {
        double e = std::fabs(filt.coeffs_b()[i] - ref_b[i]);
        if (e > max_b)
        {
            max_b = e;
        }
    }
    double max_a = 0.0;
    for (size_t i = 0; i < NC; ++i)
    {
        double e = std::fabs(filt.coeffs_a()[i] - ref_a[i]);
        if (e > max_a)
        {
            max_a = e;
        }
    }
    double max_s = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double y = filt(1.0);
        double e = std::fabs(y - ref_step[i]);
        if (e > max_s)
        {
            max_s = e;
        }
    }
    return {max_b, max_a, max_s};
}

// Step-only check for libraries that use SOS internally (b/a not available).
// StepFn() must return the next output sample for a unit-step input.
template <typename StepFn>
static AccResult check_step(StepFn step_fn, const double (&ref_step)[kStepLen])
{
    double max_s = 0.0;
    for (int i = 0; i < kStepLen; ++i)
    {
        double e = std::fabs(step_fn() - ref_step[i]);
        if (e > max_s)
        {
            max_s = e;
        }
    }
    return {-1.0, -1.0, max_s};
}
