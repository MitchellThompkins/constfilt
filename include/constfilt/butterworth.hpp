#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

struct LowPass
{
};

struct HighPass
{
};

template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class Butterworth : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");

  public:
    // Construct from filter specification; all math is constexpr.
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method>(compute_continuous_tf(cutoff_hz),
                                     sample_rate_hz)
    {
    }

  private:
    // Computes the continuous-time Butterworth transfer function.
    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // --- Low-pass
    // -------------------------------------------------------------
    //
    // Numerator: b[N] = wc^N, all other b[k] = 0  (DC gain = 1)
    //
    // Denominator: a[k] = p[k] * wc^k where p[] are the normalized (wc=1)
    // Butterworth polynomial coefficients in descending order (p[0] = 1).
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        b[N] = consteig::pow(wc, static_cast<int>(N));

        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[k] * consteig::pow(wc, static_cast<int>(k));
    }

    // --- High-pass
    // ------------------------------------------------------------
    //
    // Derived from the LPF via the LP-to-HP frequency transformation s → wc/s.
    //
    // Numerator: b[0] = 1, all other b[k] = 0  (high-frequency gain = 1)
    //
    // Denominator: a[k] = p[N-k] * wc^k - the normalized Butterworth
    // coefficients in reversed order, scaled by wc^k (p[0] = 1 so a[0] = 1,
    // keeping the denominator monic).
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        b[0] = static_cast<T>(1);

        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[N - k] * consteig::pow(wc, static_cast<int>(k));
    }

    // Normalized Butterworth denominator coefficients (wc=1, monic).
    // Fills result in descending power order: [1, p[N-1], ..., p[0]]
    //
    // Computed by multiplying out (s - p_k) for each Butterworth pole:
    //   p_k = cos(theta_k) + j*sin(theta_k),  theta_k = pi*(2k+N-1)/(2N)
    //         k = 1..N
    // Coefficients are real by construction (poles come in conjugate pairs,
    // or are real for odd N).
    static constexpr void butterworth_poly_coeffs(T (&result)[N + 1u])
    {
        using Cx = consteig::Complex<T>;

        // poly[i] holds the coefficient of s^i (ascending order).
        // Start with the constant polynomial 1.
        Cx poly[N + 1u]{};
        poly[0] = Cx{static_cast<T>(1), static_cast<T>(0)};

        // Multiply by (s - p_k) for k = 1..N.
        for (consteig::Size k = 1u; k <= N; ++k)
        {
            T theta = static_cast<T>(CONSTFILT_PI) *
                      static_cast<T>(2u * k + N - 1u) / static_cast<T>(2u * N);
            Cx pole{consteig::cos(theta), consteig::sin(theta)};

            // In-place multiply by (s - pole), traversing backwards
            // to avoid overwriting values still needed this iteration.
            for (consteig::Size j = k; j > 0u; --j)
            {
                poly[j] = poly[j - 1u] - pole * poly[j];
            }
            poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - pole * poly[0];
        }

        // Convert ascending-order complex to descending-order real.
        // The imaginary parts are zero (or near-zero numerical noise).
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            result[i] = poly[N - i].real;
        }
    }
};

// Convenience aliases for first-order RC-equivalent filters.
template <typename T, typename Method = ZOH>
using FirstOrderLowPass = Butterworth<T, 1u, Method, LowPass>;

template <typename T, typename Method = ZOH>
using FirstOrderHighPass = Butterworth<T, 1u, Method, HighPass>;

} // namespace constfilt

#endif // CONSTFILT_BUTTERWORTH_HPP
