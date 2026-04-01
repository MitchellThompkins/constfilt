#ifndef CONSTFILT_GAUSSIAN_HPP
#define CONSTFILT_GAUSSIAN_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Gaussian analog-prototype filter of order N (1 ≤ N ≤ 8).
//
// The denominator D_N(s) satisfies
//   D_N(s) * D_N(-s) = E_N(s²)
// where E_N(x) = Σ_{k=0}^{N} x^k/k!  (truncated exp).
// This makes |H(jω)|² = 1/E_N(-ω²) ≈ exp(-ω²) near ω = 0.
//
// The prototype is frequency-scaled so that the −3 dB point is at cutoff_hz.
// Poles are found via companion-matrix eigendecomposition of the polynomial
// Q_N(u) = Σ_{k=0}^{N} (-1)^k u^k/k!  (roots u_k = poles-squared of D_N).
template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class Gaussian : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 1u && N <= 8u, "Gaussian order must be 1–8");

  public:
    constexpr Gaussian(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method>(compute_continuous_tf(cutoff_hz),
                                     sample_rate_hz)
    {
    }

  private:
    using Cx = consteig::Complex<T>;

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // Constexpr factorial (exact for N ≤ 8, so 2N ≤ 16; no overflow for double).
    static constexpr T fact(consteig::Size n)
    {
        T v = static_cast<T>(1);
        for (consteig::Size i = 2u; i <= n; ++i)
            v *= static_cast<T>(i);
        return v;
    }

    // Complex square root: returns the root with non-negative imaginary part;
    // caller negates to select LHP.
    static constexpr Cx cx_sqrt(const Cx &z)
    {
        T r = consteig::abs(z); // |z|
        T re = consteig::sqrt((r + z.real) / static_cast<T>(2));
        T im = consteig::sqrt((r - z.real) / static_cast<T>(2));
        if (z.imag < static_cast<T>(0))
            im = -im;
        return {re, im};
    }

    // Build the Frobenius companion matrix for the monic form of
    //   Q_N(u) = Σ_{k=0}^N (-1)^k u^k / k!
    // Monic form: multiply Q_N by (-1)^N * N!  →  M_N(u) = u^N + Σ_{k<N} b_k u^k
    // with b_k = (-1)^{N+k} * N!/k!.
    static constexpr consteig::Matrix<T, N, N> make_companion()
    {
        consteig::Matrix<T, N, N> C{};
        // Sub-diagonal: C(i+1, i) = 1
        for (consteig::Size i = 0; i + 1u < N; ++i)
            C(i + 1u, i) = static_cast<T>(1);
        // Last column: -b_k where b_k = (-1)^{N+k} * N!/k!
        const T factN = fact(N);
        for (consteig::Size k = 0; k < N; ++k)
        {
            T sign = ((N + k) % 2u == 0u) ? static_cast<T>(1) : static_cast<T>(-1);
            T bk = sign * factN / fact(k);
            C(k, N - 1u) = -bk;
        }
        return C;
    }

    // Evaluate E_N(x) = Σ_{k=0}^N x^k/k! (used for −3 dB bisection).
    static constexpr T eval_EN(T x)
    {
        T sum = static_cast<T>(0);
        T xk = static_cast<T>(1);
        T fk = static_cast<T>(1);
        for (consteig::Size k = 0; k <= N; ++k)
        {
            sum += xk / fk;
            xk *= x;
            fk *= static_cast<T>(k + 1u);
        }
        return sum;
    }

    // Compute the normalized Gaussian denominator polynomial coefficients
    // (descending power, monic, −3 dB at ω = 1).
    static constexpr void gaussian_poly_coeffs(T (&result)[N + 1u])
    {
        // --- 1. Companion eigenvalues → roots u_k of Q_N -----------------------
        auto C = make_companion();
        auto evals = consteig::eigenvalues(C); // Matrix<Cx, N, 1>

        // --- 2. Poles s_k = −sqrt(u_k) (choose LHP root) ----------------------
        Cx poles[N]{};
        for (consteig::Size k = 0; k < N; ++k)
        {
            Cx sq = cx_sqrt(evals(k, 0));
            // We want Re(s_k) < 0: negate if real part is non-negative.
            if (sq.real >= static_cast<T>(0))
                sq = {-sq.real, -sq.imag};
            poles[k] = sq;
        }

        // --- 3. Build poly from poles (ascending, then flip to descending) ------
        // Same pattern as butterworth_poly_coeffs.
        Cx poly[N + 1u]{};
        poly[0] = Cx{static_cast<T>(1), static_cast<T>(0)};
        for (consteig::Size k = 0; k < N; ++k)
        {
            const Cx &p = poles[k];
            for (consteig::Size j = k + 1u; j > 0u; --j)
                poly[j] = poly[j - 1u] - p * poly[j];
            poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - p * poly[0];
        }
        // poly[i] is the coefficient of s^i (ascending). Flip to descending.
        T proto[N + 1u]{};
        for (consteig::Size i = 0; i <= N; ++i)
            proto[i] = poly[N - i].real;

        // --- 4. Find −3 dB frequency via bisection on E_N(ω²) = 2 -------------
        // E_N(x) = 2 where x = ω_3². E_N is strictly increasing for x ≥ 0.
        T lo = static_cast<T>(0);
        T hi = static_cast<T>(20); // generous upper bound (ω² ≤ 20 for N ≤ 8)
        for (int iter = 0; iter < 64; ++iter)
        {
            T mid = (lo + hi) / static_cast<T>(2);
            if (eval_EN(mid) < static_cast<T>(2))
                lo = mid;
            else
                hi = mid;
        }
        const T w3sq = (lo + hi) / static_cast<T>(2);
        const T w3 = consteig::sqrt(w3sq);

        // --- 5. Scale: a_norm[k] = proto[k] * w3^k ----------------------------
        for (consteig::Size k = 0; k <= N; ++k)
            result[k] = proto[k] * consteig::pow(w3, static_cast<int>(k));
    }

    // --- Low-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        T p[N + 1u]{};
        gaussian_poly_coeffs(p);
        // Denormalize: a[k] = p[k] * (wc)^k  (substitute s → s/wc, same pattern
        // as Butterworth but p already accounts for the w3 normalization above).
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[k] * consteig::pow(wc, static_cast<int>(k));
        // DC gain = 1
        b[N] = a[N];
    }

    // --- High-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        T p[N + 1u]{};
        gaussian_poly_coeffs(p);
        // HPF: a[k] = p[N-k] * wc^k
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[N - k] * consteig::pow(wc, static_cast<int>(k));
        // HF gain = 1
        b[0] = static_cast<T>(1);
    }
};

} // namespace constfilt

#endif // CONSTFILT_GAUSSIAN_HPP
