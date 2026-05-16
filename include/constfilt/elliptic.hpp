#ifndef CONSTFILT_ELLIPTIC_HPP
#define CONSTFILT_ELLIPTIC_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include "vendor/consteig/consteig.hpp"
#include "vendor/gcem_wrapper.hpp"

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
// Algorithm follows Octave's ncauer (theta-function / q-series path), with one
// deviation: Step 2 uses the modular identity q = q1^(1/N) instead of
// ncauer's iterative degree-equation solver. Steps 1 and 3 onward are
// identical. All coefficient math is constexpr.
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
    using Complex = consteig::Complex<T>;

    // Number of complex-conjugate pole/zero pairs: floor(N/2).
    static constexpr consteig::Size M{N / 2u};

    // ln(10)/10: converts dB power ratio to natural-log exponent (10^(x/10) =
    // exp(x*ln10/10)).
    static constexpr double LN10_OVER_10{0.23025850929940457};
    // AGM iteration count for elliptic_K; 64 rounds gives full double
    // precision.
    static constexpr int AGM_ITERATIONS{64};
    // Truncation depth for theta-function and nome q-series. Since q < 1,
    // terms decay as q^(n^2) and are below double precision well before n=30.
    // Any value >= ~15 would give the same result; 30 is a conservative margin.
    static constexpr int SERIES_TERMS{30};
    // Coefficients of the nome q-series: q = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13.
    static constexpr int NOME_COEFF_Q0_5{2};
    static constexpr int NOME_COEFF_Q0_9{15};
    static constexpr int NOME_COEFF_Q0_13{150};

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

    // Complete elliptic integral of the first kind K(k) via AGM
    // (https://dlmf.nist.gov/19.2#ii).
    //   K(k) = pi / (2 * AGM(1, sqrt(1-k^2)))
    static constexpr T elliptic_K(T k)
    {
        T a = static_cast<T>(1);
        T b = gcem::sqrt(static_cast<T>(1) - k * k);
        for (int i = 0; i < AGM_ITERATIONS; ++i)
        {
            T a2 = (a + b) / static_cast<T>(2);
            T b2 = gcem::sqrt(a * b);
            a = a2;
            b = b2;
        }
        return static_cast<T>(CONSTFILT_PI) / (static_cast<T>(2) * a);
    }

    // Convert a power ratio in dB to linear: 10^(x/10) = exp(x * ln10/10).
    static constexpr T from_db10(T x)
    {
        return gcem::exp(x * static_cast<T>(LN10_OVER_10));
    }

    // Nome q from modulus k (ncauer q-series approximation).
    //   q0 = 0.5 * (1 - sqrt(k')) / (1 + sqrt(k'))
    //   q  = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13
    static constexpr T compute_nome(T k)
    {
        const T kp = gcem::sqrt(static_cast<T>(1) - k * k);
        const T sqrt_kp = gcem::sqrt(kp);
        const T q0 = static_cast<T>(0.5) * (static_cast<T>(1) - sqrt_kp) /
                     (static_cast<T>(1) + sqrt_kp);
        const T q0_2 = q0 * q0;
        const T q0_4 = q0_2 * q0_2;
        const T q0_5 = q0_4 * q0;
        const T q0_9 = q0_5 * q0_4;
        const T q0_13 = q0_9 * q0_4;
        return q0 + static_cast<T>(NOME_COEFF_Q0_5) * q0_5 +
               static_cast<T>(NOME_COEFF_Q0_9) * q0_9 +
               static_cast<T>(NOME_COEFF_Q0_13) * q0_13;
    }

    // Recover modulus k from nome q (invert q = exp(-pi*K(k')/K(k)),
    // where K(k) and K(k') are the complete elliptic integrals
    // https://dlmf.nist.gov/19.2#ii).
    // No closed form exists for this inversion, so Jacobi theta functions
    // (https://dlmf.nist.gov/20.2#i) are used -- power series in q whose
    // ratio gives k exactly via the identity k = theta2^2/theta3^2
    // (https://dlmf.nist.gov/22.2, Whittaker & Watson ch. 22, Zverev s4.3):
    //   theta2(q) = 2*q^(1/4) * sum_{n=0}^{inf} q^{n(n+1)}
    //   theta3(q) = 1 + 2*sum_{n=1}^{inf} q^{n^2}
    //   k = (theta2/theta3)^2
    static constexpr T modulus_from_nome(T q)
    {
        const T q14 = gcem::sqrt(gcem::sqrt(q)); // q^(1/4)
        const T q2 = q * q;

        // theta2(0,q) power series: 2*q^(1/4) * sum_{n=0}^{inf} q^{n(n+1)}
        // Derived from the DLMF form 2*sum q^{(n+1/2)^2} by factoring out
        // q^(1/4) since (n+1/2)^2 = n(n+1) + 1/4.
        T theta2 = static_cast<T>(0);
        T qpow = static_cast<T>(1); // q^(n*(n+1)), starts at q^0 when n=0
        T q_2n = static_cast<T>(1);
        for (int n = 0; n <= SERIES_TERMS; ++n)
        {
            if (n > 0)
            {
                q_2n *= q2;
                qpow *= q_2n;
            }
            theta2 += qpow;
        }
        theta2 *= static_cast<T>(2) * q14;

        // theta3(0,q) power series: 1 + 2*sum_{n=1}^{inf} q^{n^2}
        T theta3 = static_cast<T>(1);
        T qpow3 = q; // q^(n^2), starting at q^1
        T q_2n1 = q; // q^(2n-1), starting at q^1
        for (int n = 1; n <= SERIES_TERMS; ++n)
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
        const T gain = from_db10(ripple_db / static_cast<T>(2)); // from_db10
                                                                 // divides by
                                                                 // 10; so to
                                                                 // divide by 20
                                                                 // we only need
                                                                 // to divide by
                                                                 // 2 here

        // hyperbolic angle encoding ripple spec, spread across N poles it is
        // the imaginary argument at which we evaluate the Jacobi theta
        // functions used from step 2.(https://dlmf.nist.gov/20.2#E1)
        const T l =
            gcem::log((gain + static_cast<T>(1)) / (gain - static_cast<T>(1))) /
            (static_cast<T>(2) * static_cast<T>(N));

        const T q2 = q * q; // q^2

        // sig01: q^(m*(m+1)) incremental via ratio q^(2m).
        T sig01 = static_cast<T>(0);
        T qpow1 = static_cast<T>(1); // q^(m*(m+1)), starts at q^0
        T q_2m = static_cast<T>(1);  // q^(2m), updated before use
        for (int m = 0; m <= SERIES_TERMS; ++m)
        {
            if (m > 0)
            {
                q_2m *= q2;    // q^2, q^4, q^6, ...
                qpow1 *= q_2m; // q^2, q^6, q^12, ...
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T x = static_cast<T>(2 * m + 1) * l;
            sig01 += sign * qpow1 * gcem::sinh(x);
        }

        // sig02: q^(m^2) incremental via ratio q^(2m-1).
        T sig02 = static_cast<T>(0);
        T qpow2 = q; // q^(1^2) = q
        T q_2m1 = q; // q^(2*1-1) = q
        for (int m = 1; m <= SERIES_TERMS; ++m)
        {
            if (m > 1)
            {
                q_2m1 *= q2;    // q^3, q^5, q^7, ...
                qpow2 *= q_2m1; // q^4, q^9, q^16, ...
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T x = static_cast<T>(2 * m) * l;
            sig02 += sign * qpow2 * gcem::cosh(x);
        }

        const T q14 = gcem::sqrt(gcem::sqrt(q)); // q^(1/4)

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
        const T q14 = gcem::sqrt(gcem::sqrt(q)); // q^(1/4)
        const T q2 = q * q;
        const T pi_mu_n = static_cast<T>(CONSTFILT_PI) * mu / static_cast<T>(N);

        // soma1: q^(m*(m+1)) incremental.
        T soma1 = static_cast<T>(0);
        T qpow1 = static_cast<T>(1);
        T q_2m = static_cast<T>(1);
        for (int m = 0; m <= SERIES_TERMS; ++m)
        {
            if (m > 0)
            {
                q_2m *= q2;
                qpow1 *= q_2m;
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T arg = static_cast<T>(2 * m + 1) * pi_mu_n;
            soma1 += sign * qpow1 * gcem::sin(arg);
        }
        soma1 *= static_cast<T>(2) * q14;

        // soma2: q^(m^2) incremental.
        T soma2 = static_cast<T>(0);
        T qpow2 = q;
        T q_2m1 = q;
        for (int m = 1; m <= SERIES_TERMS; ++m)
        {
            if (m > 1)
            {
                q_2m1 *= q2;
                qpow2 *= q_2m1;
            }
            const T sign =
                (m % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);
            const T arg = static_cast<T>(2 * m) * pi_mu_n;
            soma2 += sign * qpow2 * gcem::cos(arg);
        }
        soma2 *= static_cast<T>(2);

        return soma1 / (static_cast<T>(1) + soma2);
    }

    // =========================================================================
    // In-place multiply polynomial poly (ascending order, currently degree
    // deg-1) by (s - root), bringing it to degree deg. Builds up the full
    // polynomial one root at a time: start with (s - r1), call with r2 to get
    // (s - r1)(s - r2), call again with r3 to get (s - r1)(s - r2)(s - r3).
    // =========================================================================
    static constexpr void poly_mul_root(Complex (&poly)[N + 1u],
                                        consteig::Size deg, Complex root)
    {
        for (consteig::Size j = deg; j > 0u; --j)
        {
            poly[j] = poly[j - 1u] - root * poly[j];
        }
        poly[0] =
            Complex{static_cast<T>(0), static_cast<T>(0)} - root * poly[0];
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
    //
    // The filter is fully determined after step 4. Steps 2-4 are the elliptic
    // function machinery that places poles and zeros for equiripple in both
    // bands. Step 1 is unit conversion; steps 5-7 are extraction and scaling.
    // =========================================================================
    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], LowPass)
    {
        // Step 1: selectivity ratio k1 = ep/es.
        //   ep = sqrt(10^(Rp/10) - 1),  es = sqrt(10^(Rs/10) - 1)
        const T ep = gcem::sqrt(from_db10(ripple_db) - static_cast<T>(1));
        const T es = gcem::sqrt(from_db10(attenuation_db) - static_cast<T>(1));
        const T k1 = ep / es;

        // Step 2: find design modulus k.
        //   q1 = exp(-pi * K(k1') / K(k1))       nome of k1
        //   q  = q1^(1/N)                         modular equation (Zverev
        //   s4.3) k  = (theta2(q) / theta3(q))^2        recover modulus from
        //   nome
        const T q1 = compute_nome(k1);
        const T q = gcem::exp(gcem::log(q1) / static_cast<T>(N));
        const T k = modulus_from_nome(q);

        // Step 3: pole-shift parameter sig0.
        //   Controls how far poles sit in the LHP, setting the equiripple
        //   level. Moving sigma further left makes the filter more dampled,
        //   which is what changes the ripple level.
        //   Derived from sn(j*K(k')*l,
        //   k) via theta series, where l = acoth(10^(Rp/20)) / N.
        const T sig0 = compute_sig0(ripple_db, q);

        // Step 4: zero and pole positions for each conjugate pair ii=1..M.
        //   ws  = 1/k                            normalized stopband edge (>1)
        //   wi  = sn(mu*K(k)/N, k)  via theta series  (zero/pole spacing in
        //                                               elliptic frequency
        //                                               space)
        //   Vi  = cn(mu*K(k)/N, k) * dn(mu*K(k)/N, k)
        //   w   = sqrt((1 + k*sig0^2)(1 + sig0^2/k))
        //
        //   zeros: +/-j * sqrt(ws) / wi          (on the imaginary axis)
        //   poles: sqrt(ws) * (-sig0*Vi +/- j*wi*w) / (1 + sig0^2*wi^2)
        //   real pole (odd N only): -sig0 * sqrt(ws)
        const T ws = static_cast<T>(1) / k;
        const T sqrt_ws = gcem::sqrt(ws);
        const T w = gcem::sqrt((static_cast<T>(1) + k * sig0 * sig0) *
                               (static_cast<T>(1) + sig0 * sig0 / k));

        Complex poly_a[N + 1u]{};
        Complex poly_b[N + 1u]{};
        poly_a[0] = Complex{static_cast<T>(1), static_cast<T>(0)};
        poly_b[0] = Complex{static_cast<T>(1), static_cast<T>(0)};

        consteig::Size deg_a = 0u;
        consteig::Size deg_b = 0u;

        for (consteig::Size ii = 1u; ii <= M; ++ii)
        {
            const T wi = compute_wi(ii, q);
            const T Vi = gcem::sqrt((static_cast<T>(1) - k * wi * wi) *
                                    (static_cast<T>(1) - wi * wi / k));

            const T omega_z = sqrt_ws / wi;
            ++deg_b;
            poly_mul_root(poly_b, deg_b, Complex{static_cast<T>(0), omega_z});
            ++deg_b;
            poly_mul_root(poly_b, deg_b, Complex{static_cast<T>(0), -omega_z});

            const T denom = static_cast<T>(1) + sig0 * sig0 * wi * wi;
            const T p_re = sqrt_ws * (-sig0 * Vi) / denom;
            const T p_im = sqrt_ws * (wi * w) / denom;

            ++deg_a;
            poly_mul_root(poly_a, deg_a, Complex{p_re, p_im});
            ++deg_a;
            poly_mul_root(poly_a, deg_a, Complex{p_re, -p_im});
        }

        if (N % 2u == 1u)
        {
            ++deg_a;
            poly_mul_root(poly_a, deg_a,
                          Complex{-sig0 * sqrt_ws, static_cast<T>(0)});
        }

        // Step 5: convert ascending -> descending; take real parts (imag ~ 0).
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            a[i] = poly_a[N - i].real;
            b[i] = poly_b[N - i].real;
        }

        // Step 6: gain normalization.
        //   odd  N -> H(0) = 1      (DC is a ripple peak)
        //   even N -> H(0) = Gp = 1/sqrt(1+ep^2)  (DC is a ripple trough)
        const T Gp =
            static_cast<T>(1) / gcem::sqrt(static_cast<T>(1) + ep * ep);
        const T H0 = (N % 2u == 1u) ? static_cast<T>(1) : Gp;
        const T gain = H0 * a[N] / b[N];
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            b[i] *= gain;
        }

        // Step 7: scale for passband cutoff wc.
        //   Substituting s -> s/wc multiplies the coefficient of s^(N-i) by
        //   wc^i.
        for (consteig::Size i = 0u; i <= N; ++i)
        {
            const T sc = gcem::pow(wc, static_cast<int>(i));
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
            const T sc = gcem::pow(wc, static_cast<int>(j));
            a[j] = a_lp[N - j] * sc;
            b[j] = b_lp[N - j] * sc;
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_ELLIPTIC_HPP
