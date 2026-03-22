#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "../../consteig/consteig.hpp"
#include "../constfilt_options.hpp"
#include "discretize.hpp"
#include "filter.hpp"

namespace constfilt
{

template <typename T, consteig::Size N>
class Butterworth : public Filter<T, N + 1u, N + 1u>
{
    static_assert(N >= 1u && N <= 8u,
                  "Butterworth order must be between 1 and 8");

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

    // Full pipeline: specs → continuous SS → ZOH discrete SS → TF.
    static constexpr TransferFunction<T, N + 1u, N + 1u>
    compute_ba(T cutoff_hz, T sample_rate_hz)
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
        auto p = butterworth_poly_coeffs();
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
    // Returns Array<T, N+1>: [1, p[N-1], ..., p[1], p[0]]
    // i.e. index 0 = leading coefficient (1), index N = constant term (1).
    static constexpr consteig::Array<T, N + 1u> butterworth_poly_coeffs()
    {
        if constexpr (N == 1u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 2u)
        {
            // sqrt(2) ≈ 1.41421356237309504
            return {static_cast<T>(1),
                    consteig::sqrt(static_cast<T>(2)),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 3u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(2),
                    static_cast<T>(2),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 4u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(2.6131259297527580),
                    static_cast<T>(3.4142135623730950),
                    static_cast<T>(2.6131259297527580),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 5u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(3.2360679774997900),
                    static_cast<T>(5.2360679774997900),
                    static_cast<T>(5.2360679774997900),
                    static_cast<T>(3.2360679774997900),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 6u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(3.8637033051562740),
                    static_cast<T>(7.4641016151377540),
                    static_cast<T>(9.1416201726388200),
                    static_cast<T>(7.4641016151377540),
                    static_cast<T>(3.8637033051562740),
                    static_cast<T>(1)};
        }
        else if constexpr (N == 7u)
        {
            return {static_cast<T>(1),
                    static_cast<T>(4.4939592074349700),
                    static_cast<T>(10.0978347107381400),
                    static_cast<T>(14.5917938970012100),
                    static_cast<T>(14.5917938970012100),
                    static_cast<T>(10.0978347107381400),
                    static_cast<T>(4.4939592074349700),
                    static_cast<T>(1)};
        }
        else // N == 8
        {
            return {static_cast<T>(1),
                    static_cast<T>(5.1258374988234580),
                    static_cast<T>(13.1370711902630900),
                    static_cast<T>(21.8461499978232800),
                    static_cast<T>(25.6884421691498800),
                    static_cast<T>(21.8461499978232800),
                    static_cast<T>(13.1370711902630900),
                    static_cast<T>(5.1258374988234580),
                    static_cast<T>(1)};
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_BUTTERWORTH_HPP
