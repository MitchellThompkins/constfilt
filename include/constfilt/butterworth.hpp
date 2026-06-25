#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "analog_filter.hpp"
#include "vendor/consteig/consteig.hpp"
#include "vendor/gcem_wrapper.hpp"

namespace constfilt
{

struct LowPass
{
};

struct HighPass
{
};

// Template parameters:
//   T          - floating-point scalar type
//   N          - filter order (>= 1)
//   Method     - TustinPW (default), TustinNW, ZOH, or MatchedZ
//   FilterType - LowPass (default) or HighPass
template <typename T, consteig::Size N, typename Method = TustinPW,
          typename FilterType = LowPass>
class Butterworth
    : public AnalogFilter<T, N, typename bind_method<T, Method>::type>
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");

    using BoundMethod = typename bind_method<T, Method>::type;

  public:
    // Construct from filter specification; all math is constexpr.
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, BoundMethod>(
              compute_continuous_tf(cutoff_hz),
              compute_factored_tf(cutoff_hz, FilterType{}), sample_rate_hz,
              make_tustin_tag(cutoff_hz, BoundMethod{}))
    {
    }

  private:
    // Computes the continuous-time Butterworth transfer function.
    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // Low-pass
    //
    // Numerator: b[N] = wc^N, all other b[k] = 0  (DC gain = 1)
    //
    // Denominator: a[k] = p[k] * wc^k where p[] are the normalized (wc=1)
    // Butterworth polynomial coefficients in descending order (p[0] = 1).
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        b[N] = gcem::pow(wc, static_cast<int>(N));

        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
        {
            a[k] = p[k] * gcem::pow(wc, static_cast<int>(k));
        }
    }

    // High-pass
    //
    // Derived from the LPF via the LP-to-HP frequency transformation s -> wc/s.
    //
    // Numerator: b[0] = 1, all other b[k] = 0  (high-frequency gain = 1)
    //
    // Denominator: a[k] = p[N-k] * wc^k. The normalized Butterworth
    // coefficients appear in reversed order, each scaled by wc^k (p[0] = 1
    // so a[0] = 1, keeping the denominator monic).
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        b[0] = static_cast<T>(1);

        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
        {
            a[k] = p[N - k] * gcem::pow(wc, static_cast<int>(k));
        }
    }

    // LP: poles at wc*exp(j*theta_k), no finite zeros, gain = wc^N.
    static constexpr FactoredTF<T, N> compute_factored_tf(T cutoff_hz, LowPass)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        FactoredTF<T, N> factored_tf{};
        factored_tf.nz = 0;
        factored_tf.gain = gcem::pow(wc, static_cast<int>(N));
        for (consteig::Size k = 1u; k <= N; ++k)
        {
            const T theta = static_cast<T>(GCEM_PI) *
                            static_cast<T>(2u * k + N - 1u) /
                            static_cast<T>(2u * N);
            factored_tf.poles[k - 1u] = {wc * gcem::cos(theta),
                                         wc * gcem::sin(theta)};
        }
        return factored_tf;
    }

    // HP: poles at wc*exp(-j*theta_k) (magnitude wc), N zeros at s=0, gain=1.
    static constexpr FactoredTF<T, N> compute_factored_tf(T cutoff_hz, HighPass)
    {
        using Complex = consteig::Complex<T>;
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        FactoredTF<T, N> factored_tf{};
        factored_tf.nz = N;
        factored_tf.gain = static_cast<T>(1);
        for (consteig::Size k = 1u; k <= N; ++k)
        {
            const T theta = static_cast<T>(GCEM_PI) *
                            static_cast<T>(2u * k + N - 1u) /
                            static_cast<T>(2u * N);
            factored_tf.poles[k - 1u] = {wc * gcem::cos(theta),
                                         -wc * gcem::sin(theta)};
            factored_tf.zeros[k - 1u] =
                Complex{static_cast<T>(0), static_cast<T>(0)};
        }
        return factored_tf;
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
        using Complex = consteig::Complex<T>;

        // poly[i] holds the coefficient of s^i (ascending order).
        // Start with the constant polynomial 1.
        Complex poly[N + 1u]{};
        poly[0] = Complex{static_cast<T>(1), static_cast<T>(0)};

        // Multiply by (s - p_k) for k = 1..N.
        for (consteig::Size k = 1u; k <= N; ++k)
        {
            const T theta = static_cast<T>(GCEM_PI) *
                            static_cast<T>(2u * k + N - 1u) /
                            static_cast<T>(2u * N);
            const Complex pole{gcem::cos(theta), gcem::sin(theta)};

            // In-place multiply by (s - pole), traversing backwards
            // to avoid overwriting values still needed this iteration.
            for (consteig::Size j = k; j > 0u; --j)
            {
                poly[j] = poly[j - 1u] - pole * poly[j];
            }
            poly[0] =
                Complex{static_cast<T>(0), static_cast<T>(0)} - pole * poly[0];
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
template <typename T, typename Method = TustinPW>
using FirstOrderLowPass = Butterworth<T, 1u, Method, LowPass>;

template <typename T, typename Method = TustinPW>
using FirstOrderHighPass = Butterworth<T, 1u, Method, HighPass>;

} // namespace constfilt

#endif // CONSTFILT_BUTTERWORTH_HPP
