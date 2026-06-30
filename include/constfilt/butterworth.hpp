#ifndef CONSTFILT_BUTTERWORTH_HPP
#define CONSTFILT_BUTTERWORTH_HPP

#include "analog_filter.hpp"
#include "vendor/consteig/consteig.hpp"
#include "vendor/gcem_wrapper.hpp"

namespace constfilt
{

template <typename T, consteig::Size N, typename Method = TustinPW,
          typename FilterType = LowPass>
class ButterworthImpl
    : public AnalogFilter<T, N, typename bind_method<T, Method>::type>
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");

    using BoundMethod = typename bind_method<T, Method>::type;

  public:
    constexpr ButterworthImpl(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, BoundMethod>(
              compute_continuous_tf(cutoff_hz),
              compute_factored_tf(cutoff_hz, FilterType{}), sample_rate_hz,
              make_tustin_tag(cutoff_hz, BoundMethod{}))
    {
    }

    // Construct with uniform damping ratio zeta across all complex pole pairs.
    // Routes through generic eigendecomposition (no FactoredTF) to avoid the
    // Vandermonde singularity that arises when N>=4 produces identical pole
    // pairs.
    constexpr ButterworthImpl(T cutoff_hz, T sample_rate_hz, T zeta)
        : AnalogFilter<T, N, BoundMethod>(
              compute_continuous_tf_zeta(cutoff_hz, zeta), sample_rate_hz,
              make_tustin_tag(cutoff_hz, BoundMethod{}))
    {
        // Uniform zeta produces repeated complex eigenvalue pairs for N>=4.
        // consteig's QR iteration cannot split a defective eigenvalue block,
        // so matrix_exp returns wrong results. MatchedZ and Tustin are fine.
        static_assert(
            !is_zoh_tag<Method>::value || N <= 3u,
            "Butterworth zeta constructor: ZOH is unreliable for N>=4 due to "
            "repeated eigenvalues. Use MatchedZ or TustinPW instead. "
            "See https://github.com/MitchellThompkins/constfilt/issues/56");
    }

  private:
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

    // Normalized Butterworth denominator coefficients for uniform damping
    // ratio. Each complex pair contributes a quadratic factor (s^2 + 2*zeta*s +
    // 1); an odd-order real pole contributes (s + 1). Pure real arithmetic, no
    // complex types. Fills result in descending power order: [1, ..., 1].
    static constexpr void butterworth_poly_coeffs_zeta(T (&result)[N + 1u],
                                                       T zeta)
    {
        // poly[i] holds the coefficient of s^i (ascending order).
        T poly[N + 1u]{};
        poly[0] = static_cast<T>(1);
        for (consteig::Size pair = 0u; pair < N / 2u; ++pair)
        {
            const consteig::Size cur = 2u * pair;
            for (consteig::Size j = cur + 2u; j > 1u; --j)
            {
                poly[j] +=
                    static_cast<T>(2) * zeta * poly[j - 1u] + poly[j - 2u];
            }
            poly[1] += static_cast<T>(2) * zeta * poly[0];
        }
        if (N % 2u == 1u)
        {
            for (consteig::Size j = N; j > 0u; --j)
            {
                poly[j] += poly[j - 1u];
            }
        }
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            result[i] = poly[N - i];
        }
    }

    static constexpr void continuous_tf_zeta(T wc, T zeta, T (&b)[N + 1u],
                                             T (&a)[N + 1u], LowPass)
    {
        b[N] = gcem::pow(wc, static_cast<int>(N));
        T p[N + 1u]{};
        butterworth_poly_coeffs_zeta(p, zeta);
        for (consteig::Size k = 0; k <= N; ++k)
        {
            a[k] = p[k] * gcem::pow(wc, static_cast<int>(k));
        }
    }

    static constexpr void continuous_tf_zeta(T wc, T zeta, T (&b)[N + 1u],
                                             T (&a)[N + 1u], HighPass)
    {
        b[0] = static_cast<T>(1);
        T p[N + 1u]{};
        butterworth_poly_coeffs_zeta(p, zeta);
        for (consteig::Size k = 0; k <= N; ++k)
        {
            a[k] = p[N - k] * gcem::pow(wc, static_cast<int>(k));
        }
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u>
    compute_continuous_tf_zeta(T cutoff_hz, T zeta)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf_zeta(wc, zeta, tf.b, tf.a, FilterType{});
        return tf;
    }
};

// Primary template: SOS cascade (default, numerically superior).
//
// Stores ceil(N/2) second-order sections (biquads), each independently
// discretized from one conjugate pole pair of the N-th order Butterworth
// design. Odd-N filters have one padded first-order section.
//
// Template parameters:
//   T          - floating-point scalar type
//   N          - filter order (>= 1)
//   Method     - TustinPW (default), TustinNW, ZOH, or MatchedZ
//   FilterType - LowPass (default) or HighPass
//   SOS        - true (default): SOS cascade; false: direct form
template <typename T, consteig::Size N, typename Method = TustinPW,
          typename FilterType = LowPass, bool SOS = true>
class Butterworth
{
    static_assert(N >= 1u, "Butterworth order must be at least 1");
    static_assert(
        !is_zoh_tag<Method>::value,
        "Butterworth SOS + ZOH: ZOH is not separable over cascaded "
        "sections and produces a different filter than full-order ZOH. "
        "Use SOS=false for ZOH, or switch to TustinPW or MatchedZ.");

    using BoundMethod = typename bind_method<T, Method>::type;

    static constexpr consteig::Size kSections = (N + 1u) / 2u;
    Filter<T, 3u, 3u> _sections[kSections]{};

  public:
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        const BoundMethod method_tag =
            make_tustin_tag(cutoff_hz, BoundMethod{});

        const consteig::Size n_pairs = N / 2u;
        for (consteig::Size i = 0u; i < n_pairs; ++i)
        {
            _sections[i] = make_complex_section(wc, i, sample_rate_hz,
                                                method_tag, FilterType{});
        }

        if constexpr (N % 2u == 1u)
        {
            _sections[kSections - 1u] =
                make_real_section(wc, sample_rate_hz, method_tag, FilterType{});
        }
    }

    // Construct with uniform damping ratio zeta across all complex pole pairs.
    // Each 2x2 section has distinct eigenvalues regardless of zeta, so all
    // discretization methods including ZOH work for any N.
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz, T zeta)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        const BoundMethod method_tag =
            make_tustin_tag(cutoff_hz, BoundMethod{});
        const T re = -zeta * wc;
        const T im = wc * gcem::sqrt(static_cast<T>(1) - zeta * zeta);

        const consteig::Size n_pairs = N / 2u;
        for (consteig::Size i = 0u; i < n_pairs; ++i)
        {
            _sections[i] = make_zeta_section(re, im, wc, sample_rate_hz,
                                             method_tag, FilterType{});
        }

        if constexpr (N % 2u == 1u)
        {
            _sections[kSections - 1u] =
                make_real_section(wc, sample_rate_hz, method_tag, FilterType{});
        }
    }

    T operator()(T x) const
    {
        for (consteig::Size i = 0u; i < kSections; ++i)
            x = _sections[i](x);
        return x;
    }

    template <consteig::Size Len>
    constexpr void operator()(const T (&input)[Len], T (&output)[Len]) const
    {
        _sections[0](input, output);
        for (consteig::Size i = 1u; i < kSections; ++i)
        {
            T buf[Len]{};
            _sections[i](output, buf);
            for (consteig::Size n = 0u; n < Len; ++n)
                output[n] = buf[n];
        }
    }

  private:
    // Build a second-order section for the i-th conjugate pole pair.
    //
    // Pole index mapping (1-indexed k = i+1):
    //   theta_k = pi*(2k+N-1)/(2N)
    // The upper-half-plane pole (positive imaginary part) is used; its
    // conjugate is implicit. All Butterworth poles have magnitude wc.

    static constexpr Filter<T, 3u, 3u> make_complex_section(
        T wc, consteig::Size pair_idx, T sample_rate_hz, BoundMethod method_tag,
        LowPass)
    {
        const T theta = static_cast<T>(GCEM_PI) *
                        static_cast<T>(2u * (pair_idx + 1u) + N - 1u) /
                        static_cast<T>(2u * N);
        const T re = wc * gcem::cos(theta);
        const T im = wc * gcem::sin(theta);
        const T mag_sq = re * re + im * im;

        // LP: b(s) = |p|^2 (unity DC gain), a(s) = (s-p)(s-conj(p))
        T b_c[3]{static_cast<T>(0), static_cast<T>(0), mag_sq};
        T a_c[3]{static_cast<T>(1), -static_cast<T>(2) * re, mag_sq};

        FactoredTF<T, 2u> factored{};
        factored.poles[0] = {re, im};
        factored.poles[1] = {re, -im};
        factored.nz = 0;
        factored.gain = mag_sq;

        TransferFunction<T, 3u, 3u> ctf{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            ctf.b[j] = b_c[j];
            ctf.a[j] = a_c[j];
        }

        const AnalogFilter<T, 2u, Method> sec{ctf, factored, sample_rate_hz,
                                              method_tag};
        T b_d[3]{};
        T a_d[3]{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            b_d[j] = sec.coeffs_b()[j];
            a_d[j] = sec.coeffs_a()[j];
        }
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    static constexpr Filter<T, 3u, 3u> make_complex_section(
        T wc, consteig::Size pair_idx, T sample_rate_hz, BoundMethod method_tag,
        HighPass)
    {
        const T theta = static_cast<T>(GCEM_PI) *
                        static_cast<T>(2u * (pair_idx + 1u) + N - 1u) /
                        static_cast<T>(2u * N);
        const T re = wc * gcem::cos(theta);
        const T im = wc * gcem::sin(theta);
        const T mag_sq = re * re + im * im;

        // HP: zeros at s=0 (N zeros total), unity high-freq gain
        T b_c[3]{static_cast<T>(1), static_cast<T>(0), static_cast<T>(0)};
        T a_c[3]{static_cast<T>(1), -static_cast<T>(2) * re, mag_sq};

        FactoredTF<T, 2u> factored{};
        factored.poles[0] = {re, im};
        factored.poles[1] = {re, -im};
        factored.nz = 2;
        factored.zeros[0] = {static_cast<T>(0), static_cast<T>(0)};
        factored.zeros[1] = {static_cast<T>(0), static_cast<T>(0)};
        factored.gain = static_cast<T>(1);

        TransferFunction<T, 3u, 3u> ctf{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            ctf.b[j] = b_c[j];
            ctf.a[j] = a_c[j];
        }

        const AnalogFilter<T, 2u, Method> sec{ctf, factored, sample_rate_hz,
                                              method_tag};
        T b_d[3]{};
        T a_d[3]{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            b_d[j] = sec.coeffs_b()[j];
            a_d[j] = sec.coeffs_a()[j];
        }
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    // Build a second-order section from an explicitly supplied pole (re, im).
    // Used by the uniform-zeta constructor where all sections share the same
    // pole position. mag_sq = wc^2 for any zeta since |p|^2 = wc^2.
    static constexpr Filter<T, 3u, 3u> make_zeta_section(T re, T im, T wc,
                                                         T sample_rate_hz,
                                                         BoundMethod method_tag,
                                                         LowPass)
    {
        const T mag_sq = wc * wc;

        T b_c[3]{static_cast<T>(0), static_cast<T>(0), mag_sq};
        T a_c[3]{static_cast<T>(1), -static_cast<T>(2) * re, mag_sq};

        FactoredTF<T, 2u> factored{};
        factored.poles[0] = {re, im};
        factored.poles[1] = {re, -im};
        factored.nz = 0;
        factored.gain = mag_sq;

        TransferFunction<T, 3u, 3u> ctf{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            ctf.b[j] = b_c[j];
            ctf.a[j] = a_c[j];
        }

        const AnalogFilter<T, 2u, Method> sec{ctf, factored, sample_rate_hz,
                                              method_tag};
        T b_d[3]{};
        T a_d[3]{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            b_d[j] = sec.coeffs_b()[j];
            a_d[j] = sec.coeffs_a()[j];
        }
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    static constexpr Filter<T, 3u, 3u> make_zeta_section(T re, T im, T wc,
                                                         T sample_rate_hz,
                                                         BoundMethod method_tag,
                                                         HighPass)
    {
        const T mag_sq = wc * wc;

        T b_c[3]{static_cast<T>(1), static_cast<T>(0), static_cast<T>(0)};
        T a_c[3]{static_cast<T>(1), -static_cast<T>(2) * re, mag_sq};

        FactoredTF<T, 2u> factored{};
        factored.poles[0] = {re, im};
        factored.poles[1] = {re, -im};
        factored.nz = 2;
        factored.zeros[0] = {static_cast<T>(0), static_cast<T>(0)};
        factored.zeros[1] = {static_cast<T>(0), static_cast<T>(0)};
        factored.gain = static_cast<T>(1);

        TransferFunction<T, 3u, 3u> ctf{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            ctf.b[j] = b_c[j];
            ctf.a[j] = a_c[j];
        }

        const AnalogFilter<T, 2u, Method> sec{ctf, factored, sample_rate_hz,
                                              method_tag};
        T b_d[3]{};
        T a_d[3]{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            b_d[j] = sec.coeffs_b()[j];
            a_d[j] = sec.coeffs_a()[j];
        }
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    // First-order real pole section for odd N, padded to biquad storage.
    // LP real pole is at s = -wc (left half plane), giving unity DC gain.
    static constexpr Filter<T, 3u, 3u> make_real_section(T wc, T sample_rate_hz,
                                                         BoundMethod method_tag,
                                                         LowPass)
    {
        TransferFunction<T, 2u, 2u> ctf{};
        ctf.b[0] = static_cast<T>(0);
        ctf.b[1] = wc;
        ctf.a[0] = static_cast<T>(1);
        ctf.a[1] = wc;

        const AnalogFilter<T, 1u, Method> sec{ctf, sample_rate_hz, method_tag};
        T b_d[3]{sec.coeffs_b()[0], sec.coeffs_b()[1], static_cast<T>(0)};
        T a_d[3]{sec.coeffs_a()[0], sec.coeffs_a()[1], static_cast<T>(0)};
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    // HP real pole at s = -wc, zero at s = 0 (LP-to-HP transform), unity
    // high-freq gain.
    static constexpr Filter<T, 3u, 3u> make_real_section(T wc, T sample_rate_hz,
                                                         BoundMethod method_tag,
                                                         HighPass)
    {
        TransferFunction<T, 2u, 2u> ctf{};
        ctf.b[0] = static_cast<T>(1);
        ctf.b[1] = static_cast<T>(0);
        ctf.a[0] = static_cast<T>(1);
        ctf.a[1] = wc;

        const AnalogFilter<T, 1u, Method> sec{ctf, sample_rate_hz, method_tag};
        T b_d[3]{sec.coeffs_b()[0], sec.coeffs_b()[1], static_cast<T>(0)};
        T a_d[3]{sec.coeffs_a()[0], sec.coeffs_a()[1], static_cast<T>(0)};
        return Filter<T, 3u, 3u>{b_d, a_d};
    }
};

// Partial specialization: direct-form realization (opt-in via SOS=false).
// Inherits from ButterworthImpl which extends AnalogFilter<T, N, Method>.
// Preserves the existing single-filter interface (coeffs_b, coeffs_a, etc.).
template <typename T, consteig::Size N, typename Method, typename FilterType>
class Butterworth<T, N, Method, FilterType, false>
    : public ButterworthImpl<T, N, Method, FilterType>
{
  public:
    constexpr Butterworth(T cutoff_hz, T sample_rate_hz)
        : ButterworthImpl<T, N, Method, FilterType>(cutoff_hz, sample_rate_hz)
    {
    }

    constexpr Butterworth(T cutoff_hz, T sample_rate_hz, T zeta)
        : ButterworthImpl<T, N, Method, FilterType>(cutoff_hz, sample_rate_hz,
                                                    zeta)
    {
    }
};

// Convenience aliases for first-order RC-equivalent filters.
// N=1 has a single real pole; SOS offers no advantage, so these use
// direct form.
template <typename T, typename Method = TustinPW>
using FirstOrderLowPass = Butterworth<T, 1u, Method, LowPass, false>;

template <typename T, typename Method = TustinPW>
using FirstOrderHighPass = Butterworth<T, 1u, Method, HighPass, false>;

} // namespace constfilt

#endif // CONSTFILT_BUTTERWORTH_HPP
