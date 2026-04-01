#ifndef CONSTFILT_LINKWITZ_RILEY_HPP
#define CONSTFILT_LINKWITZ_RILEY_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Linkwitz–Riley filter of total order N (must be even).
//
// H_LR(s) = [H_BW_{N/2}(s)]^2  — the square of an (N/2)th-order Butterworth.
//
// LPF: b[N] = wc^N, all other b = 0  (DC gain = 1)
// HPF: b[0] = 1,    all other b = 0  (HF gain = 1)
//
// The denominator is the self-convolution of the (N/2)th-order Butterworth
// denominator polynomial (descending power order, monic).
template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class LinkwitzRiley : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 2u && N % 2u == 0u,
                  "Linkwitz-Riley order must be even and >= 2");
    static constexpr consteig::Size M = N / 2u; // half-order

  public:
    constexpr LinkwitzRiley(T crossover_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method>(compute_continuous_tf(crossover_hz),
                                     sample_rate_hz)
    {
    }

  private:
    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T crossover_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * crossover_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // Normalized Butterworth denominator for order M (descending power, wc=1).
    // Same algorithm as Butterworth::butterworth_poly_coeffs.
    static constexpr void bw_poly(T (&result)[M + 1u])
    {
        using Cx = consteig::Complex<T>;
        Cx poly[M + 1u]{};
        poly[0] = Cx{static_cast<T>(1), static_cast<T>(0)};
        for (consteig::Size k = 1u; k <= M; ++k)
        {
            T theta = static_cast<T>(CONSTFILT_PI) *
                      static_cast<T>(2u * k + M - 1u) / static_cast<T>(2u * M);
            Cx pole{consteig::cos(theta), consteig::sin(theta)};
            for (consteig::Size j = k; j > 0u; --j)
                poly[j] = poly[j - 1u] - pole * poly[j];
            poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - pole * poly[0];
        }
        for (consteig::Size i = 0u; i <= M; ++i)
            result[i] = poly[M - i].real;
    }

    // Self-convolve a degree-M polynomial (descending) → degree-N result.
    static constexpr void poly_self_convolve(const T (&p)[M + 1u],
                                             T (&result)[N + 1u])
    {
        for (consteig::Size k = 0; k <= N; ++k)
        {
            T sum = static_cast<T>(0);
            for (consteig::Size i = 0; i <= M; ++i)
            {
                const consteig::Size j = k - i;
                if (j <= M)
                    sum += p[i] * p[j];
            }
            result[k] = sum;
        }
    }

    // --- Low-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        T p[M + 1u]{};
        bw_poly(p);

        // Scale: substitute s -> s/wc in the half-order poly.
        // Descending coeffs p[k] correspond to s^(M-k), so scaled coeff = p[k]*wc^k...
        // wait: p[0]=1 is coeff of s^M, p[M] is coeff of s^0.
        // After substituting s->s/wc: p[k]*(s/wc)^(M-k) -> p[k]/wc^(M-k) * s^(M-k).
        // Multiply through by wc^M to make monic: p[k]*wc^k ... but we want monic.
        // Actually: for the LR filter we square the full scaled Butterworth poly.
        // The Butterworth LPF at wc: a_bw[k] = p[k] * wc^(M-k) ... descending.
        // p[0]=1 -> a_bw[0]=wc^0... hmm let me re-examine.
        //
        // butterworth.hpp: a[k] = p[k] * wc^k  where k=0..N and p is in
        // descending order. So p[0]*wc^0, p[1]*wc^1, ..., p[N]*wc^N.
        T p_scaled[M + 1u]{};
        for (consteig::Size k = 0; k <= M; ++k)
            p_scaled[k] = p[k] * consteig::pow(wc, static_cast<int>(k));

        poly_self_convolve(p_scaled, a);
        b[N] = consteig::pow(wc, static_cast<int>(N)); // DC gain = 1
    }

    // --- High-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        T p[M + 1u]{};
        bw_poly(p);

        // HPF half-order: a_hp[k] = p[M-k] * wc^k  (same as Butterworth HPF)
        T p_hp[M + 1u]{};
        for (consteig::Size k = 0; k <= M; ++k)
            p_hp[k] = p[M - k] * consteig::pow(wc, static_cast<int>(k));

        poly_self_convolve(p_hp, a);
        b[0] = static_cast<T>(1); // HF gain = 1
    }
};

} // namespace constfilt

#endif // CONSTFILT_LINKWITZ_RILEY_HPP
