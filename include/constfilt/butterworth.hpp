#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "constfilt_options.hpp"
#include "discretize.hpp"
#include "filter.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

template <typename T, consteig::Size N>
class Butterworth : public Filter<T, N + 1u, N + 1u>
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");

  public:
    // Construct from filter specification; all math is constexpr.
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
        : Butterworth(compute_ba(cutoff_hz, sample_rate_hz))
    {
    }

  private:
    // Delegating constructor: receives pre-computed transfer function.
    constexpr explicit Butterworth(TransferFunction<T, N + 1u, N + 1u> tf)
        : Filter<T, N + 1u, N + 1u>(tf.b, tf.a)
    {
    }

    // Full pipeline: specs -> continuous SS -> ZOH discrete SS -> TF.
    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_ba(
        T cutoff_hz, T sample_rate_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        const T Ts = static_cast<T>(1) / sample_rate_hz;
        auto sys_c = build_continuous_ss(wc);
        auto sys_d = zoh_discretize(sys_c, Ts, ZOH{});
        return ss_to_tf(sys_d);
    }

    // Build controllable-canonical-form state-space for an N-th order
    // Butterworth LP with cutoff wc (rad/s).
    //
    // Denormalized denominator (monic):
    //   d(s) = s^N + p[N-1]*wc * s^{N-1} + ... + p[0]*wc^N
    // where p[] are the normalized (wc=1) Butterworth coefficients.
    //
    // Companion matrix A:
    //   A[i][i+1] = 1   for i = 0..N-2
    //   A[N-1][k] = -(p[k] * wc^{N-k})   for k = 0..N-1
    //
    // B = [0, ..., 0, wc^N]^T,  C = [1, 0, ..., 0],  D = 0
    static constexpr StateSpace<T, N> build_continuous_ss(T wc)
    {
        StateSpace<T, N> sys{};

        // Fill super-diagonal with 1
        for (consteig::Size i = 0; i < N - 1u; ++i)
        {
            sys.A(i, i + 1u) = static_cast<T>(1);
        }

        // Last row: -p[k] * wc^{N-k}
        T p[N + 1u]{};
        butterworth_poly_coeffs(p);
        for (consteig::Size k = 0; k < N; ++k)
        {
            sys.A(N - 1u, k) =
                -p[k] * consteig::pow(wc, static_cast<int>(N - k));
        }

        // B: last entry = wc^N
        sys.B(N - 1u, 0) = consteig::pow(wc, static_cast<int>(N));

        // C: first entry = 1
        sys.C(0, 0) = static_cast<T>(1);

        // D = 0 (already zero-initialized)
        return sys;
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
