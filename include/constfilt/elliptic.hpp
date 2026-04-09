#ifndef CONSTFILT_ELLIPTIC_HPP
#define CONSTFILT_ELLIPTIC_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>
#include <gcem.hpp>

namespace constfilt
{

// Elliptic (Cauer) IIR filter of order N.
//
// Equiripple in both passband and stopband; minimum-order for a given
// Rp / Rs specification.
//
// Template parameters:
//   T          - floating-point scalar type
//   N          - filter order (>= 1)
//   Method     - discretization method (ZOH default)
//   FilterType - LowPass (default) or HighPass
//
// Constructor parameters:
//   cutoff_hz      - passband edge (-Rp dB point)
//   ripple_db      - passband ripple Rp in dB  (e.g. 0.5)
//   attenuation_db - stopband attenuation Rs in dB (e.g. 40)
//   sample_rate_hz - sample rate in Hz
//
// Algorithm follows Octave's ncauer (theta-function / q-series path).
// All coefficient math is constexpr.
template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class Elliptic : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 1u, "Elliptic order must be at least 1");

  public:
    constexpr Elliptic(T cutoff_hz, T ripple_db, T attenuation_db,
                       T sample_rate_hz)
        : AnalogFilter<T, N, Method>(
              compute_continuous_tf(cutoff_hz, ripple_db, attenuation_db),
              sample_rate_hz)
    {
    }

  private:
    using Cx = consteig::Complex<T>;

    // Number of complex-conjugate pole/zero pairs: floor(N/2).
    static constexpr consteig::Size M = N / 2u;

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz, T ripple_db, T attenuation_db)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        elliptic_tf(wc, ripple_db, attenuation_db, tf.b, tf.a, FilterType{});
        return tf;
    }

    // =========================================================================
    // Math helpers
    // =========================================================================

    // Complete elliptic integral of the first kind K(k) via AGM.
    //   K(k) = pi / (2 * AGM(1, sqrt(1-k^2)))
    static constexpr T elliptic_K(T k)
    {
        T a = static_cast<T>(1);
        T b = consteig::sqrt(static_cast<T>(1) - k * k);
        for (int i = 0; i < 64; ++i)
        {
            T a2 = (a + b) / static_cast<T>(2);
            T b2 = consteig::sqrt(a * b);
            a = a2;
            b = b2;
        }
        return static_cast<T>(CONSTFILT_PI) / (static_cast<T>(2) * a);
    }

    // Convert a power ratio in dB to linear: 10^(x/10) = exp(x * ln10/10).
    static constexpr T from_db10(T x)
    {
        return consteig::exp(x * static_cast<T>(0.23025850929940457));
    }

    // Nome q from modulus k (ncauer q-series approximation).
    //   q0 = 0.5 * (1 - sqrt(k')) / (1 + sqrt(k'))
    //   q  = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13
    static constexpr T compute_nome(T k)
    {
        const T kp = consteig::sqrt(static_cast<T>(1) - k * k);
        const T sqrt_kp = consteig::sqrt(kp);
        const T q0 = static_cast<T>(0.5) * (static_cast<T>(1) - sqrt_kp) /
                     (static_cast<T>(1) + sqrt_kp);
        const T q0_2 = q0 * q0;
        const T q0_4 = q0_2 * q0_2;
        const T q0_5 = q0_4 * q0;
        const T q0_9 = q0_5 * q0_4;
        const T q0_13 = q0_9 * q0_4;
        return q0 + static_cast<T>(2) * q0_5 + static_cast<T>(15) * q0_9 +
               static_cast<T>(150) * q0_13;
    }

    // Recover modulus k from nome q via Jacobi theta functions.
    //   theta2(q) = 2*q^(1/4) * sum_{n=0}^{inf} q^{n(n+1)}
    //   theta3(q) = 1 + 2*sum_{n=1}^{inf} q^{n^2}
    //   k = (theta2/theta3)^2
    static constexpr T modulus_from_nome(T q)
    {
        const T q14 = consteig::sqrt(consteig::sqrt(q));
        const T q2 = q * q;

        T theta2 = static_cast<T>(0);
        T qpow = static_cast<T>(1); // q^(n*(n+1))
        T q_2n = static_cast<T>(1);
        for (int n = 0; n <= 30; ++n)
        {
            if (n > 0)
            {
                q_2n *= q2;
                qpow *= q_2n;
            }
            theta2 += qpow;
        }
        theta2 *= static_cast<T>(2) * q14;

        T theta3 = static_cast<T>(1);
        T qpow3 = q; // q^(n^2), starting at q^1
        T q_2n1 = q; // q^(2n-1), starting at q^1
        for (int n = 1; n <= 30; ++n)
        {
            if (n > 1)
            {
                q_2n1 *= q2;
                qpow3 *= q_2n1;
            }
            theta3 += static_cast<T>(2) * qpow3;
        }

        const T ratio = theta2 / theta3;
        return ratio * ratio;
    }

    // =========================================================================
    // Pole-shift parameter sig0 via theta-function series (ncauer algorithm).
    //
    //   l     = (1/(2N)) * log((10^(0.05*Rp) + 1) / (10^(0.05*Rp) - 1))
    //   sig01 = sum_{m=0..30} (-1)^m * q^(m(m+1)) * sinh((2m+1)*l)
    //   sig02 = sum_{m=1..30} (-1)^m * q^(m^2)    * cosh(2*m*l)
    //   sig0  = abs(2 * q^(1/4) * sig01 / (1 + 2*sig02))
    // =========================================================================
    static constexpr T compute_sig0(T ripple_db, T q)
    {
        const T g = from_db10(ripple_db / static_cast<T>(2));
        const T l =
            gcem::log((g + static_cast<T>(1)) / (g - static_cast<T>(1))) /
            (static_cast<T>(2) * static_cast<T>(N));

        const T q2 = q * q;

        // sig01: q^(m*(m+1)) incremental via ratio q^(2m).
        T sig01 = static_cast<T>(0);
        T qpow1 = static_cast<T>(1); // q^(m*(m+1)), starts at q^0
        T q_2m = static_cast<T>(1);  // q^(2m), updated before use
        for (int m = 0; m <= 30; ++m)
        {
            if (m > 0)
            {
                q_2m *= q2;    // q^2, q^4, q^6, ...
                qpow1 *= q_2m; // q^2, q^6, q^12, ...
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T x = static_cast<T>(2 * m + 1) * l;
            const T ex = consteig::exp(x);
            const T sh = (ex - static_cast<T>(1) / ex) / static_cast<T>(2);
            sig01 += sign * qpow1 * sh;
        }

        // sig02: q^(m^2) incremental via ratio q^(2m-1).
        T sig02 = static_cast<T>(0);
        T qpow2 = q; // q^(1^2) = q
        T q_2m1 = q; // q^(2*1-1) = q
        for (int m = 1; m <= 30; ++m)
        {
            if (m > 1)
            {
                q_2m1 *= q2;    // q^3, q^5, q^7, ...
                qpow2 *= q_2m1; // q^4, q^9, q^16, ...
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T x = static_cast<T>(2 * m) * l;
            const T ex = consteig::exp(x);
            const T ch = (ex + static_cast<T>(1) / ex) / static_cast<T>(2);
            sig02 += sign * qpow2 * ch;
        }

        const T q14 = consteig::sqrt(consteig::sqrt(q));
        T sig0 = static_cast<T>(2) * q14 * sig01 /
                 (static_cast<T>(1) + static_cast<T>(2) * sig02);
        return (sig0 < static_cast<T>(0)) ? -sig0 : sig0;
    }

    // =========================================================================
    // Compute zero position wi via theta-function series (ncauer algorithm).
    //
    //   mu = ii (odd N), mu = ii - 0.5 (even N)
    //   soma1 = sum_{m=0..30} 2*q^(1/4)*(-1)^m*q^(m(m+1))*sin((2m+1)*pi*mu/N)
    //   soma2 = sum_{m=1..30} 2*(-1)^m*q^(m^2)*cos(2*m*pi*mu/N)
    //   wi    = soma1 / (1 + soma2)
    // =========================================================================
    static constexpr T compute_wi(consteig::Size ii, T q)
    {
        const T mu = (N % 2u == 1u) ? static_cast<T>(ii)
                                    : static_cast<T>(ii) - static_cast<T>(0.5);
        const T q14 = consteig::sqrt(consteig::sqrt(q));
        const T q2 = q * q;
        const T pi_mu_n = static_cast<T>(CONSTFILT_PI) * mu / static_cast<T>(N);

        // soma1: q^(m*(m+1)) incremental.
        T soma1 = static_cast<T>(0);
        T qpow1 = static_cast<T>(1);
        T q_2m = static_cast<T>(1);
        for (int m = 0; m <= 30; ++m)
        {
            if (m > 0)
            {
                q_2m *= q2;
                qpow1 *= q_2m;
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T arg = static_cast<T>(2 * m + 1) * pi_mu_n;
            soma1 += sign * qpow1 * consteig::sin(arg);
        }
        soma1 *= static_cast<T>(2) * q14;

        // soma2: q^(m^2) incremental.
        T soma2 = static_cast<T>(0);
        T qpow2 = q;
        T q_2m1 = q;
        for (int m = 1; m <= 30; ++m)
        {
            if (m > 1)
            {
                q_2m1 *= q2;
                qpow2 *= q_2m1;
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T arg = static_cast<T>(2 * m) * pi_mu_n;
            soma2 += sign * qpow2 * consteig::cos(arg);
        }
        soma2 *= static_cast<T>(2);

        return soma1 / (static_cast<T>(1) + soma2);
    }

    // =========================================================================
    // In-place multiply polynomial poly (ascending order, currently degree
    // deg-1) by (s - root), bringing it to degree deg.
    // =========================================================================
    static constexpr void poly_mul_root(Cx (&poly)[N + 1u], consteig::Size deg,
                                        Cx root)
    {
        for (consteig::Size j = deg; j > 0u; --j)
            poly[j] = poly[j - 1u] - root * poly[j];
        poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - root * poly[0];
    }

    // =========================================================================
    // Low-pass elliptic transfer function (ncauer theta-function algorithm).
    //
    // Steps:
    //   1. Compute nome q from k1 via modular identity: q = q1^(1/N).
    //   2. Recover design modulus k from q via theta functions.
    //   3. Pole-shift sig0 via theta series.
    //   4. Zero positions wi via theta series.
    //   5. Build s-domain polynomials from poles/zeros, scale by sqrt(ws).
    //   6. Gain normalization: H(0)=1 (odd N), H(0)=Gp (even N).
    //   7. Scale for passband cutoff wc.
    // =========================================================================
    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], LowPass)
    {
        // Ripple factors.
        const T ep = consteig::sqrt(from_db10(ripple_db) - static_cast<T>(1));
        const T es =
            consteig::sqrt(from_db10(attenuation_db) - static_cast<T>(1));

        const T k1 = ep / es;

        // Nome q via modular identity: q = q1^(1/N) where q1 = nome(k1).
        // This bypasses the degree equation solver entirely, giving full
        // machine precision for q regardless of how k is found.
        const T q1 = compute_nome(k1);
        const T q = consteig::exp(gcem::log(q1) / static_cast<T>(N));

        // Recover design modulus k from q via theta functions.
        const T k = modulus_from_nome(q);
        const T ws = static_cast<T>(1) / k;
        const T sqrt_ws = consteig::sqrt(ws);

        // Pole-shift parameter sig0.
        const T sig0 = compute_sig0(ripple_db, q);

        // Derived quantity w.
        const T w = consteig::sqrt((static_cast<T>(1) + k * sig0 * sig0) *
                                   (static_cast<T>(1) + sig0 * sig0 / k));

        // Ascending-order polynomials: poly[i] = coefficient of s^i.
        Cx poly_a[N + 1u]{};
        Cx poly_b[N + 1u]{};
        poly_a[0] = Cx{static_cast<T>(1), static_cast<T>(0)};
        poly_b[0] = Cx{static_cast<T>(1), static_cast<T>(0)};

        consteig::Size deg_a = 0u;
        consteig::Size deg_b = 0u;

        for (consteig::Size ii = 1u; ii <= M; ++ii)
        {
            const T wi = compute_wi(ii, q);
            const T Vi = consteig::sqrt((static_cast<T>(1) - k * wi * wi) *
                                        (static_cast<T>(1) - wi * wi / k));

            // Zeros (scaled): +/-j * sqrt(ws) / wi
            const T omega_z = sqrt_ws / wi;
            ++deg_b;
            poly_mul_root(poly_b, deg_b, Cx{static_cast<T>(0), omega_z});
            ++deg_b;
            poly_mul_root(poly_b, deg_b, Cx{static_cast<T>(0), -omega_z});

            // Poles (scaled): sqrt(ws)*(-sig0*Vi +/- j*wi*w)/(1+sig0^2*wi^2)
            const T denom = static_cast<T>(1) + sig0 * sig0 * wi * wi;
            const T p_re = sqrt_ws * (-sig0 * Vi) / denom;
            const T p_im = sqrt_ws * (wi * w) / denom;

            ++deg_a;
            poly_mul_root(poly_a, deg_a, Cx{p_re, p_im});
            ++deg_a;
            poly_mul_root(poly_a, deg_a, Cx{p_re, -p_im});
        }

        // Real pole for odd N: -sig0 * sqrt(ws).
        if (N % 2u == 1u)
        {
            ++deg_a;
            poly_mul_root(poly_a, deg_a,
                          Cx{-sig0 * sqrt_ws, static_cast<T>(0)});
        }

        // Convert ascending -> descending; take real parts (imag ~ 0).
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            a[i] = poly_a[N - i].real;
            b[i] = poly_b[N - i].real;
        }

        // Gain normalization to match Octave's ellipap:
        //   odd  N -> H(0) = 1
        //   even N -> H(0) = Gp = 1/sqrt(1+ep^2)
        const T Gp =
            static_cast<T>(1) / consteig::sqrt(static_cast<T>(1) + ep * ep);
        const T H0 = (N % 2u == 1u) ? static_cast<T>(1) : Gp;
        const T gain = H0 * a[N] / b[N];
        for (consteig::Size i = 0u; i <= N; ++i)
            b[i] *= gain;

        // Scale for passband cutoff wc: coefficient of s^i picks up wc^i.
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            const T sc = consteig::pow(wc, static_cast<int>(i));
            a[i] *= sc;
            b[i] *= sc;
        }
    }

    // =========================================================================
    // High-pass via LP-to-HP transform.
    //
    // Computes the normalized LP prototype (wc = 1 rad/s), then:
    //   a_hp[j] = a_lp[N-j] * wc^j
    //   b_hp[j] = b_lp[N-j] * wc^j
    // =========================================================================
    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], HighPass)
    {
        T b_lp[N + 1u]{};
        T a_lp[N + 1u]{};
        elliptic_tf(static_cast<T>(1), ripple_db, attenuation_db, b_lp, a_lp,
                    LowPass{});

        for (consteig::Size j = 0u; j <= N; ++j)
        {
            const T sc = consteig::pow(wc, static_cast<int>(j));
            a[j] = a_lp[N - j] * sc;
            b[j] = b_lp[N - j] * sc;
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_ELLIPTIC_HPP
