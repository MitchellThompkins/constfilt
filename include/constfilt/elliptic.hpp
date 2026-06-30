#ifndef CONSTFILT_ELLIPTIC_HPP
#define CONSTFILT_ELLIPTIC_HPP

#include "analog_filter.hpp"
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
//   Method     - TustinPW (default), TustinNW, ZOH, or MatchedZ
//   FilterType - LowPass (default) or HighPass
//   SOS        - true (default): SOS cascade; false: direct form
//
// Constructor parameters:
//   cutoff_hz      - passband edge (-Rp dB point)
//   ripple_db      - passband ripple Rp in dB  (e.g. 0.5)
//   attenuation_db - stopband attenuation Rs in dB (e.g. 40)
//   sample_rate_hz - sample rate in Hz
//
// The implementation follows Octave's ncauer (theta-function / q-series path),
// with one deviation: Step 2 uses the modular identity q = q1^(1/N) instead of
// ncauer's iterative degree-equation solver. Steps 1 and 3 onward are
// identical. All coefficient math is constexpr.
template <typename T, consteig::Size N, typename Method = TustinPW,
          typename FilterType = LowPass>
class EllipticImpl
    : public AnalogFilter<T, N, typename bind_method<T, Method>::type>
{
    static_assert(N >= 1u, "Elliptic order must be at least 1");

    using BoundMethod = typename bind_method<T, Method>::type;

  public:
    constexpr EllipticImpl(T cutoff_hz, T ripple_db, T attenuation_db,
                           T sample_rate_hz)
        : AnalogFilter<T, N, BoundMethod>(
              compute_continuous_tf(cutoff_hz, ripple_db, attenuation_db),
              compute_factored_tf(cutoff_hz, ripple_db, attenuation_db,
                                  FilterType{}),
              sample_rate_hz, make_tustin_tag(cutoff_hz, BoundMethod{}))
    {
    }

    static constexpr FactoredTF<T, N> compute_factored_tf(T cutoff_hz,
                                                          T ripple_db,
                                                          T attenuation_db,
                                                          LowPass);
    static constexpr FactoredTF<T, N> compute_factored_tf(T cutoff_hz,
                                                          T ripple_db,
                                                          T attenuation_db,
                                                          HighPass);

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
    // Coefficients of the nome q-series: q = q0 + 2*q0^5 + 15*q0^9 +
    // 150*q0^13.
    static constexpr int NOME_COEFF_Q0_5{2};
    static constexpr int NOME_COEFF_Q0_9{15};
    static constexpr int NOME_COEFF_Q0_13{150};

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz, T ripple_db, T attenuation_db)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        elliptic_tf(wc, ripple_db, attenuation_db, tf.b, tf.a, FilterType{});
        return tf;
    }

    // Math helpers

    // Complete elliptic integral of the first kind K(k) via AGM
    // (https://dlmf.nist.gov/19.2#ii).
    //   K(k) = pi / (2 * AGM(1, sqrt(1-k^2)))
    static constexpr T elliptic_K(T k)
    {
        T a = static_cast<T>(1);
        T b = gcem::sqrt(static_cast<T>(1) - k * k);
        for (int i = 0; i < AGM_ITERATIONS; ++i)
        {
            const T a2 = (a + b) / static_cast<T>(2);
            const T b2 = gcem::sqrt(a * b);
            a = a2;
            b = b2;
        }
        return static_cast<T>(GCEM_PI) / (static_cast<T>(2) * a);
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

    // Pole-shift parameter sig0 via theta-function series (ncauer algorithm).
    //
    //   l     = (1/(2N)) * log((10^(0.05*Rp) + 1) / (10^(0.05*Rp) - 1))
    //   sig01 = sum_{m=0..30} (-1)^m * q^(m(m+1)) * sinh((2m+1)*l)
    //   sig02 = sum_{m=1..30} (-1)^m * q^(m^2)    * cosh(2*m*l)
    //   sig0  = abs(2 * q^(1/4) * sig01 / (1 + 2*sig02))
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

        const T sig0 = static_cast<T>(2) * q14 * sig01 /
                       (static_cast<T>(1) + static_cast<T>(2) * sig02);

        return (sig0 < static_cast<T>(0)) ? -sig0 : sig0;
    }

    // Compute zero position wi via theta-function series (ncauer algorithm).
    //
    //   mu = ii (odd N), mu = ii - 0.5 (even N)
    //   soma1 = sum_{m=0..30} 2*q^(1/4)*(-1)^m*q^(m(m+1))*sin((2m+1)*pi*mu/N)
    //   soma2 = sum_{m=1..30} 2*(-1)^m*q^(m^2)*cos(2*m*pi*mu/N)
    //   wi    = soma1 / (1 + soma2)
    static constexpr T compute_wi(consteig::Size ii, T q)
    {
        const T mu = (N % 2u == 1u) ? static_cast<T>(ii)
                                    : static_cast<T>(ii) - static_cast<T>(0.5);
        const T q14 = gcem::sqrt(gcem::sqrt(q)); // q^(1/4)
        const T q2 = q * q;
        const T pi_mu_n = static_cast<T>(GCEM_PI) * mu / static_cast<T>(N);

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

    // In-place multiply polynomial poly (ascending order, currently degree
    // deg-1) by (s - root), bringing it to degree deg. Builds up the full
    // polynomial one root at a time: start with (s - r1), call with r2 to get
    // (s - r1)(s - r2), call again with r3 to get (s - r1)(s - r2)(s - r3).
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

    // Fills poles[] and zeros[] with the normalized (wc=1) prototype poles and
    // zeros derived from the elliptic machinery (steps 3-4 of elliptic_tf).
    // Poles: M conjugate pairs + optional real pole for odd N.
    // Zeros: M conjugate pairs on the imaginary axis (+/-j*omega_z).
    static constexpr void compute_prototype_poles_zeros(
        T q, T sig0, T k, Complex (&poles)[N], consteig::Size &pole_cnt,
        Complex (&zeros)[N], consteig::Size &zero_cnt)
    {
        const T ws = static_cast<T>(1) / k;
        const T sqrt_ws = gcem::sqrt(ws);
        const T w = gcem::sqrt((static_cast<T>(1) + k * sig0 * sig0) *
                               (static_cast<T>(1) + sig0 * sig0 / k));

        pole_cnt = 0u;
        zero_cnt = 0u;

        for (consteig::Size ii = 1u; ii <= M; ++ii)
        {
            const T wi = compute_wi(ii, q);
            const T Vi = gcem::sqrt((static_cast<T>(1) - k * wi * wi) *
                                    (static_cast<T>(1) - wi * wi / k));

            const T omega_z = sqrt_ws / wi;
            zeros[zero_cnt++] = Complex{static_cast<T>(0), omega_z};
            zeros[zero_cnt++] = Complex{static_cast<T>(0), -omega_z};

            const T denom = static_cast<T>(1) + sig0 * sig0 * wi * wi;
            const T p_re = sqrt_ws * (-sig0 * Vi) / denom;
            const T p_im = sqrt_ws * (wi * w) / denom;
            poles[pole_cnt++] = Complex{p_re, p_im};
            poles[pole_cnt++] = Complex{p_re, -p_im};
        }

        if (N % 2u == 1u)
        {
            poles[pole_cnt++] = Complex{-sig0 * sqrt_ws, static_cast<T>(0)};
        }
    }

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

    // High-pass via LP-to-HP transform.
    //
    // Computes the normalized LP prototype (wc = 1 rad/s), then:
    //   a_hp[j] = a_lp[N-j] * wc^j
    //   b_hp[j] = b_lp[N-j] * wc^j
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

    // LP FactoredTF: prototype poles/zeros scaled by wc; gain from polynomial.
    // Poles: M conjugate pairs at indices [0..2M-1]; real pole at [2M] for odd
    // N. Zeros: M conjugate pairs at indices [0..2M-1], pure imaginary.
    static constexpr FactoredTF<T, N> compute_factored_tf_impl(T cutoff_hz,
                                                               T ripple_db,
                                                               T attenuation_db,
                                                               LowPass)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;
        const T ep = gcem::sqrt(from_db10(ripple_db) - static_cast<T>(1));
        const T es = gcem::sqrt(from_db10(attenuation_db) - static_cast<T>(1));
        const T k1 = ep / es;
        const T q1 = compute_nome(k1);
        const T q = gcem::exp(gcem::log(q1) / static_cast<T>(N));
        const T k = modulus_from_nome(q);
        const T sig0 = compute_sig0(ripple_db, q);

        Complex poles_proto[N]{};
        Complex zeros_proto[N]{};
        consteig::Size pole_cnt = 0u;
        consteig::Size zero_cnt = 0u;
        compute_prototype_poles_zeros(q, sig0, k, poles_proto, pole_cnt,
                                      zeros_proto, zero_cnt);

        FactoredTF<T, N> factored_tf{};
        factored_tf.nz = zero_cnt;
        for (consteig::Size i = 0u; i < pole_cnt; ++i)
        {
            factored_tf.poles[i] =
                Complex{wc * poles_proto[i].real, wc * poles_proto[i].imag};
        }
        for (consteig::Size i = 0u; i < zero_cnt; ++i)
        {
            factored_tf.zeros[i] =
                Complex{wc * zeros_proto[i].real, wc * zeros_proto[i].imag};
        }

        // Gain from the polynomial TF (captures the normalization from steps
        // 6-7).
        T b_tmp[N + 1u]{};
        T a_tmp[N + 1u]{};
        elliptic_tf(wc, ripple_db, attenuation_db, b_tmp, a_tmp, LowPass{});
        consteig::Size d_b = 0u;
        while (d_b <= N && b_tmp[d_b] == static_cast<T>(0))
        {
            ++d_b;
        }
        factored_tf.gain =
            (d_b > N) ? static_cast<T>(0) : b_tmp[d_b] / a_tmp[0];

        return factored_tf;
    }

    // HP FactoredTF: derive from LP prototype via LP-to-HP transform (s ->
    // wc/s). LP pole p_lp -> HP pole wc/p_lp; LP zero j*omega_z -> HP zero
    // -j*wc/omega_z. For odd N: one extra zero at s=0 (LP strictly proper -> HP
    // has zero at origin).
    // Poles: M conjugate pairs at [0..2M-1]; real pole at [2M] for odd N.
    // Zeros: M conjugate pairs at [0..2M-1]; zero at origin at [2M] for odd N.
    static constexpr FactoredTF<T, N> compute_factored_tf_impl(T cutoff_hz,
                                                               T ripple_db,
                                                               T attenuation_db,
                                                               HighPass)
    {
        const T wc = static_cast<T>(2) * static_cast<T>(GCEM_PI) * cutoff_hz;

        // Normalized LP prototype at wc=1 (cutoff_hz = 1/(2*pi)).
        const T norm_cutoff =
            static_cast<T>(1) / (static_cast<T>(2) * static_cast<T>(GCEM_PI));
        const FactoredTF<T, N> lp = compute_factored_tf_impl(
            norm_cutoff, ripple_db, attenuation_db, LowPass{});

        FactoredTF<T, N> factored_tf{};

        // HP poles: wc / lp_pole  (complex division)
        for (consteig::Size i = 0u; i < N; ++i)
        {
            const Complex &p = lp.poles[i];
            const T denom_sq = p.real * p.real + p.imag * p.imag;
            factored_tf.poles[i] =
                Complex{wc * p.real / denom_sq, -wc * p.imag / denom_sq};
        }

        // HP zeros: wc / lp_zero  (LP zeros are pure imaginary: {0,
        // +/-omega_z})
        consteig::Size hp_nz = 0u;
        for (consteig::Size i = 0u; i < lp.nz; ++i)
        {
            const Complex &z = lp.zeros[i];
            const T denom_sq = z.real * z.real + z.imag * z.imag;
            factored_tf.zeros[hp_nz++] =
                Complex{wc * z.real / denom_sq, -wc * z.imag / denom_sq};
        }
        // For odd N: LP->HP adds a zero at s=0 (from the strictly-proper LP).
        if (N % 2u == 1u)
        {
            factored_tf.zeros[hp_nz++] =
                Complex{static_cast<T>(0), static_cast<T>(0)};
        }
        factored_tf.nz = hp_nz;

        // Gain from the HP polynomial TF.
        T b_tmp[N + 1u]{};
        T a_tmp[N + 1u]{};
        elliptic_tf(wc, ripple_db, attenuation_db, b_tmp, a_tmp, HighPass{});
        consteig::Size d_b = 0u;
        while (d_b <= N && b_tmp[d_b] == static_cast<T>(0))
        {
            ++d_b;
        }
        factored_tf.gain =
            (d_b > N) ? static_cast<T>(0) : b_tmp[d_b] / a_tmp[0];

        return factored_tf;
    }
};

template <typename T, consteig::Size N, typename Method, typename FilterType>
constexpr FactoredTF<T, N> EllipticImpl<
    T, N, Method, FilterType>::compute_factored_tf(T cutoff_hz, T ripple_db,
                                                   T attenuation_db,
                                                   LowPass tag)
{
    return compute_factored_tf_impl(cutoff_hz, ripple_db, attenuation_db, tag);
}

template <typename T, consteig::Size N, typename Method, typename FilterType>
constexpr FactoredTF<T, N> EllipticImpl<
    T, N, Method, FilterType>::compute_factored_tf(T cutoff_hz, T ripple_db,
                                                   T attenuation_db,
                                                   HighPass tag)
{
    return compute_factored_tf_impl(cutoff_hz, ripple_db, attenuation_db, tag);
}

// Primary template: SOS cascade (default, numerically superior).
//
// Stores ceil(N/2) second-order sections, each independently discretized
// from one conjugate pole/zero pair of the N-th order Elliptic design.
// Odd-N filters have one padded first-order section.
template <typename T, consteig::Size N, typename Method = TustinPW,
          typename FilterType = LowPass, bool SOS = true>
class Elliptic
{
    static_assert(N >= 1u, "Elliptic order must be at least 1");
    static_assert(
        !is_zoh_tag<Method>::value,
        "Elliptic SOS + ZOH: ZOH is not separable over cascaded "
        "sections and produces a different filter than full-order ZOH. "
        "Use SOS=false for ZOH, or switch to TustinPW or MatchedZ.");

    using BoundMethod = typename bind_method<T, Method>::type;
    using Impl = EllipticImpl<T, N, Method, FilterType>;

    static constexpr consteig::Size kSections = (N + 1u) / 2u;
    Filter<T, 3u, 3u> _sections[kSections]{};

  public:
    constexpr Elliptic(T cutoff_hz, T ripple_db, T attenuation_db,
                       T sample_rate_hz)
    {
        const FactoredTF<T, N> factored = Impl::compute_factored_tf(
            cutoff_hz, ripple_db, attenuation_db, FilterType{});
        const consteig::Size n_pairs = N / 2u;

        if constexpr (is_matchedz_tag<BoundMethod>::value)
        {
            const T Ts = static_cast<T>(1) / sample_rate_hz;
            const T w_c_global = compute_mz_test_freq<T, N>(
                factored.poles, factored.zeros, factored.nz, Ts);
            for (consteig::Size i = 0u; i < n_pairs; ++i)
            {
                _sections[i] =
                    make_complex_section_mz(factored, i, Ts, w_c_global);
            }
            if constexpr (N % 2u == 1u)
            {
                const BoundMethod method_tag =
                    make_tustin_tag(cutoff_hz, BoundMethod{});
                _sections[kSections - 1u] = make_real_section(
                    factored, sample_rate_hz, method_tag, FilterType{});
            }
        }
        else
        {
            const BoundMethod method_tag =
                make_tustin_tag(cutoff_hz, BoundMethod{});
            for (consteig::Size i = 0u; i < n_pairs; ++i)
            {
                _sections[i] = make_complex_section(factored, i, sample_rate_hz,
                                                    method_tag);
            }
            if constexpr (N % 2u == 1u)
            {
                _sections[kSections - 1u] = make_real_section(
                    factored, sample_rate_hz, method_tag, FilterType{});
            }
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
    // Pure Matched-Z second-order section: n_extra=0, caller-supplied w_c.
    // Used when Method=MatchedZ so that all sections share the global test
    // frequency derived from the full-filter FactoredTF, avoiding
    // inconsistent per-section gain matching (especially for HP odd N, where
    // the full filter has a zero at s=0 that shifts the test frequency).
    static constexpr Filter<T, 3u, 3u> make_complex_section_mz(
        const FactoredTF<T, N> &factored, consteig::Size pair_idx, T Ts,
        T w_c_global)
    {
        using Complex = consteig::Complex<T>;
        const T section_gain =
            (pair_idx == 0u) ? factored.gain : static_cast<T>(1);
        Complex poles[2u]{factored.poles[2u * pair_idx],
                          factored.poles[2u * pair_idx + 1u]};
        Complex zeros[2u]{factored.zeros[2u * pair_idx],
                          factored.zeros[2u * pair_idx + 1u]};
        const auto dtf = matched_z_assemble_sos<T, 2u>(
            poles, zeros, 2u, section_gain, Ts, w_c_global, 0u);
        T b_d[3]{dtf.b[0], dtf.b[1], dtf.b[2]};
        T a_d[3]{dtf.a[0], dtf.a[1], dtf.a[2]};
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    // Build a second-order section for conjugate pole/zero pair i.
    //
    // Poles and zeros are stored as conjugate pairs in the FactoredTF:
    //   poles[2*i], poles[2*i+1]  -- conjugate pair
    //   zeros[2*i], zeros[2*i+1]  -- conjugate pair (pure imaginary for both
    //                                 LP and HP elliptic)
    //
    // The overall filter gain is placed on section 0 so the cascade product
    // equals the true filter response. Remaining sections use unit gain in the
    // factored sense.
    static constexpr Filter<T, 3u, 3u> make_complex_section(
        const FactoredTF<T, N> &factored, consteig::Size pair_idx,
        T sample_rate_hz, BoundMethod method_tag)
    {
        const T re = factored.poles[2u * pair_idx].real;
        const T im = factored.poles[2u * pair_idx].imag;
        const T mag_sq_p = re * re + im * im;

        // zeros are pure imaginary; use the imaginary part of the first zero
        const T omega_z = factored.zeros[2u * pair_idx].imag;
        const T mag_sq_z = omega_z * omega_z;

        const T section_gain =
            (pair_idx == 0u) ? factored.gain : static_cast<T>(1);

        // 2nd-order section TF: section_gain*(s^2+omega_z^2)/(s^2-2*re*s+|p|^2)
        T b_c[3]{section_gain, static_cast<T>(0), section_gain * mag_sq_z};
        T a_c[3]{static_cast<T>(1), -static_cast<T>(2) * re, mag_sq_p};

        FactoredTF<T, 2u> sec_factored{};
        sec_factored.poles[0] = factored.poles[2u * pair_idx];
        sec_factored.poles[1] = factored.poles[2u * pair_idx + 1u];
        sec_factored.zeros[0] = factored.zeros[2u * pair_idx];
        sec_factored.zeros[1] = factored.zeros[2u * pair_idx + 1u];
        sec_factored.nz = 2u;
        sec_factored.gain = section_gain;

        TransferFunction<T, 3u, 3u> ctf{};
        for (consteig::Size j = 0u; j < 3u; ++j)
        {
            ctf.b[j] = b_c[j];
            ctf.a[j] = a_c[j];
        }

        const AnalogFilter<T, 2u, Method> sec{ctf, sec_factored, sample_rate_hz,
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

    // LP real pole at factored.poles[N-1]; no zeros.
    // TF: 1/(s - r)  where r < 0. Unit contribution in cascade; overall gain
    // is on section 0.
    static constexpr Filter<T, 3u, 3u> make_real_section(
        const FactoredTF<T, N> &factored, T sample_rate_hz,
        BoundMethod method_tag, LowPass)
    {
        const T r = factored.poles[N - 1u].real; // negative
        TransferFunction<T, 2u, 2u> ctf{};
        ctf.b[0] = static_cast<T>(0);
        ctf.b[1] = static_cast<T>(1);
        ctf.a[0] = static_cast<T>(1);
        ctf.a[1] = -r; // |r| > 0

        const AnalogFilter<T, 1u, Method> sec{ctf, sample_rate_hz, method_tag};
        T b_d[3]{sec.coeffs_b()[0], sec.coeffs_b()[1], static_cast<T>(0)};
        T a_d[3]{sec.coeffs_a()[0], sec.coeffs_a()[1], static_cast<T>(0)};
        return Filter<T, 3u, 3u>{b_d, a_d};
    }

    // HP real pole at factored.poles[N-1]; zero at origin
    // (factored.zeros[N-1]). TF: s/(s - r)  where r < 0. Unity high-frequency
    // gain.
    static constexpr Filter<T, 3u, 3u> make_real_section(
        const FactoredTF<T, N> &factored, T sample_rate_hz,
        BoundMethod method_tag, HighPass)
    {
        const T r = factored.poles[N - 1u].real; // negative
        TransferFunction<T, 2u, 2u> ctf{};
        ctf.b[0] = static_cast<T>(1);
        ctf.b[1] = static_cast<T>(0);
        ctf.a[0] = static_cast<T>(1);
        ctf.a[1] = -r; // |r| > 0

        const AnalogFilter<T, 1u, Method> sec{ctf, sample_rate_hz, method_tag};
        T b_d[3]{sec.coeffs_b()[0], sec.coeffs_b()[1], static_cast<T>(0)};
        T a_d[3]{sec.coeffs_a()[0], sec.coeffs_a()[1], static_cast<T>(0)};
        return Filter<T, 3u, 3u>{b_d, a_d};
    }
};

// Partial specialization: direct-form realization (opt-in via SOS=false).
template <typename T, consteig::Size N, typename Method, typename FilterType>
class Elliptic<T, N, Method, FilterType, false>
    : public EllipticImpl<T, N, Method, FilterType>
{
  public:
    constexpr Elliptic(T cutoff_hz, T ripple_db, T attenuation_db,
                       T sample_rate_hz)
        : EllipticImpl<T, N, Method, FilterType>(cutoff_hz, ripple_db,
                                                 attenuation_db, sample_rate_hz)
    {
    }
};

} // namespace constfilt

#endif // CONSTFILT_ELLIPTIC_HPP
