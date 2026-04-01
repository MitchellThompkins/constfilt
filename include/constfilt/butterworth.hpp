#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Butterworth poles are always strictly in the left half-plane, so no
// stability check is needed (CheckStab = false).
template <typename T, consteig::Size N, typename Method = ZOH>
class Butterworth : public AnalogFilter<T, N, Method, false>
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");

  public:
    // Construct from filter specification; all math is constexpr.
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method, false>(
              compute_tf(cutoff_hz, sample_rate_hz))
    {
    }

  private:
    // Full pipeline: specs -> continuous TF -> continuous SS -> discrete TF.
    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_tf(
        T cutoff_hz, T sample_rate_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        T b_c[N + 1u]{};
        T a_c[N + 1u]{};
        continuous_tf(wc, b_c, a_c);
        return AnalogFilter<T, N, Method, false>::discretize(b_c, a_c,
                                                             sample_rate_hz);
    }

    // Continuous-time lowpass Butterworth TF coefficients (descending power
    // order).
    //
    // Numerator: b[N] = wc^N, all other b[k] = 0  (DC gain = 1)
    //
    // Denominator: a[k] = p[k] * wc^k where p[] are the normalized (wc=1)
    // Butterworth polynomial coefficients in descending order (p[0] = 1).
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u])
    {
        b[N] = consteig::pow(wc, static_cast<int>(N));

        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[k] * consteig::pow(wc, static_cast<int>(k));
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

} // namespace constfilt

#endif // CONSTFILT_BUTTERWORTH_HPP
