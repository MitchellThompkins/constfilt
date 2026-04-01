#ifndef CONSTFILT_BESSEL_HPP
#define CONSTFILT_BESSEL_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Bessel (Thomson) filter of order N (1 ≤ N ≤ 8).
//
// The denominator polynomial is the Nth-order Bessel polynomial with
// coefficients c[k] = (2N-k)! / (2^(N-k) * k! * (N-k)!), k = 0..N.
// The prototype is frequency-normalized so that the −3 dB point lies at
// cutoff_hz. This preserves the maximally-flat group-delay property while
// providing a consistent cutoff definition across all orders.
//
// LPF and HPF variants share the same denominator (HPF via s → wc/s).
template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class Bessel : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 1u && N <= 8u, "Bessel order must be 1–8");

  public:
    constexpr Bessel(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method>(compute_continuous_tf(cutoff_hz),
                                     sample_rate_hz)
    {
    }

  private:
    // Normalized −3 dB frequencies of the delay-normalized Bessel prototype
    // (ω₀ = 1) for orders 1–8.  Multiply by wc to get the actual −3 dB edge.
    // Values from standard filter design tables (Zverev, Williams & Taylor).
    static constexpr T bessel_w3db(consteig::Size n)
    {
        // clang-format off
        constexpr double tbl[8] = {
            1.0000000000, 1.3617021277, 1.7557037788, 2.1139203656,
            2.4274334658, 2.7033990536, 2.9517256720, 3.1796268716
        };
        // clang-format on
        return static_cast<T>(tbl[n - 1u]);
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // Bessel polynomial B_N coefficients in descending power order (monic).
    // result[k] is the coefficient of s^(N-k), so result[0]=1, result[N]=B_N(0).
    // Raw (delay-normalized): c[k] = (2N-k)! / (2^(N-k) * k! * (N-k)!)
    // Stored descending: result[k] = c[N-k].
    static constexpr void bessel_poly_coeffs(T (&result)[N + 1u])
    {
        for (consteig::Size k = 0; k <= N; ++k)
        {
            // coefficient of s^(N-k), i.e. c[k] in ascending order
            // c[k] = (2N-k)! / (2^(N-k) * k! * (N-k)!)
            // Compute as a floating-point product to avoid integer overflow.
            T val = static_cast<T>(1);
            // numerator: product (N+1)..(2N-k)  divided by 2^(N-k) / k! / (N-k)!
            // Use the recurrence: c[k] = c[k-1] * (2N-k) / (2*(N-k+1)*(N-k))
            // Start from c[0] = (2N)! / (2^N * N! * N!)  and recurse.
            // Simpler: compute directly via product formula.
            // numerator: prod_{i=k+1}^{2N-k} ... too complex; use iterative.
            // Direct: c[k] = factorial(2N-k) / (2^(N-k) * factorial(k) * factorial(N-k))
            // Build using floating-point to avoid overflow for N<=8.

            // factorial(2N-k)
            T num = static_cast<T>(1);
            for (consteig::Size i = 1u; i <= 2u * N - k; ++i)
                num *= static_cast<T>(i);

            // 2^(N-k)
            T pow2 = static_cast<T>(1);
            for (consteig::Size i = 0u; i < N - k; ++i)
                pow2 *= static_cast<T>(2);

            // factorial(k)
            T factk = static_cast<T>(1);
            for (consteig::Size i = 1u; i <= k; ++i)
                factk *= static_cast<T>(i);

            // factorial(N-k)
            T factnk = static_cast<T>(1);
            for (consteig::Size i = 1u; i <= N - k; ++i)
                factnk *= static_cast<T>(i);

            val = num / (pow2 * factk * factnk);

            // result in descending power: result[j] = c[N-j]
            result[N - k] = val;
        }
    }

    // --- Low-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        // Delay-normalized polynomial coefficients (descending power, monic).
        T p[N + 1u]{};
        bessel_poly_coeffs(p);
        const T w3 = bessel_w3db(N);

        // Frequency-denormalize: substitute s → s * w3 / wc so that −3 dB is at wc.
        // coefficient of s^(N-k) in p[k] → p[k] * (w3/wc)^(N-k)
        // After substitution and scaling to keep monic:
        //   a[k] = p[k] * (w3/wc)^(N-k)  ... but butterworth.hpp uses a[k]*wc^k
        // Equivalent: a[k] = p[k] * w3^(N-k) / wc^(N-k)
        // Then scale entire denominator by (wc/w3)^N to match standard form:
        //   a[k] *= (wc/w3)^k * (wc/w3)^(... ) ... let's keep it simple:
        // Standard analog-filter convention (constfilt): a[k] is the coefficient
        // of s^(N-k) in descending order. After s → s*w3/wc:
        //   coefficient of s^(N-k) = p[k] * (w3/wc)^(N-k)
        // Multiply through by (wc/w3)^N to restore monic s^N coefficient = 1:
        //   a[k] = p[k] * (w3/wc)^(N-k) * (wc/w3)^N = p[k] * w3^k / wc^k * wc^N/w3^N
        //        = p[k] * (w3/wc)^(N-2k) ...
        // This is getting complicated. Use the cleaner approach:
        //   substitute s → s/(wc/w3) = s*w3/wc, then multiply through by (wc/w3)^N.
        // After both steps: a[k] = p[k] * (wc/w3)^(N-k) ... wait that's for s->s*wc/w3.
        //
        // Easiest: a[k] = p[k] * (wc/w3)^k  following butterworth.hpp pattern
        // (substituting s → s/wc and using the w3-normalized prototype).
        // The prototype with w3-normalization has a[k] = p[k]*w3^k in the
        // standard butterworth scaling. Substituting wc for w3:
        //   a[k] = p[k] * (wc/w3)^(N-k) * w3^(N-k) ... nope.
        //
        // Let's do it directly. The raw Bessel poly B_N(s) has:
        //   B_N(j*w3) is the point where |H|^2 = 0.5 (for unit numerator).
        // So the -3dB normalized polynomial is B_N(s * w3) renormalized to monic.
        // We want -3dB at wc, so use B_N(s * w3 / wc), renormalized.
        // B_N(s * w3 / wc) = sum_k c[k] * (w3/wc)^k * s^k  (ascending).
        // In descending, coefficient of s^(N-k) = c[k] * (w3/wc)^k = p[N-k] * (w3/wc)^k.
        // So: a[k] (descending index k) = p[k] * (w3/wc)^(N-k).
        // Then normalize to make monic (divide by a[0] = p[0]*(w3/wc)^N = (w3/wc)^N):
        //   a[k] /= (w3/wc)^N  → a[k] = p[k] * (wc/w3)^N / (wc/w3)^k ... sigh.
        //
        // Actually: p[0] = 1 (monic already). a[0] = 1*(w3/wc)^N which is not 1.
        // To fix: a[k] = p[k] * (w3/wc)^(N-k) / (w3/wc)^N = p[k] * (wc/w3)^k.
        //
        // So: a[k] = p[k] * (wc/w3)^k.  Same pattern as Butterworth with wc→wc/w3!
        const T wc_norm = wc / w3;
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[k] * consteig::pow(wc_norm, static_cast<int>(k));

        // DC gain = 1: b[N] = a[N]  (H(0) = b[N]/a[N] = 1)
        b[N] = a[N];
    }

    // --- High-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        T p[N + 1u]{};
        bessel_poly_coeffs(p);
        const T w3 = bessel_w3db(N);
        const T wc_norm = wc / w3;

        // HPF via LP→HP: a[k] = p[N-k] * (wc/w3)^k  (reversed + scaled)
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[N - k] * consteig::pow(wc_norm, static_cast<int>(k));

        // HF gain = 1
        b[0] = static_cast<T>(1);
    }
};

} // namespace constfilt

#endif // CONSTFILT_BESSEL_HPP
