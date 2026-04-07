#ifndef CONSTFILT_ELLIPTIC_HPP
#define CONSTFILT_ELLIPTIC_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

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
// Algorithm follows Orfanidis, "Lecture Notes on Elliptic Filter Design."
// All coefficient math is constexpr and uses only consteig primitives.
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

    static constexpr TransferFunction<T, N + 1u, N + 1u>
    compute_continuous_tf(T cutoff_hz, T ripple_db, T attenuation_db)
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
            a    = a2;
            b    = b2;
        }
        return static_cast<T>(CONSTFILT_PI) / (static_cast<T>(2) * a);
    }

    // Convert a power ratio in dB to linear: 10^(x/10) = exp(x * ln10/10).
    // ln(10)/10 = 0.23025850929940457
    static constexpr T from_db10(T x)
    {
        return consteig::exp(x * static_cast<T>(0.23025850929940457));
    }

    // =========================================================================
    // Jacobi elliptic functions sn, cn, dn via descending Landen transform.
    //
    // u   - absolute argument
    // k   - modulus in (0, 1)
    //
    // Descending Landen step:
    //   k_{n+1} = (1 - sqrt(1-k_n^2)) / (1 + sqrt(1-k_n^2))
    //
    // Scale accumulation:
    //   scale = product(1 + k_{n+1}) = K(k) / K(k_final)
    //
    // Bottom (k_final ~ 0): sn = sin(u/scale), cn = cos(u/scale), dn = 1.
    //
    // Ascending formulas (from level n+1 back to level n):
    //   denom  = 1 + k_{n+1} * sn^2
    //   sn_new = (1 + k_{n+1}) * sn  / denom
    //   cn_new = cn * dn              / denom
    //   dn_new = (1 - k_{n+1} * sn^2) / denom
    // =========================================================================
    static constexpr void jac_scd(T u, T k, T &sn_r, T &cn_r, T &dn_r)
    {
        if (k < static_cast<T>(1e-12))
        {
            sn_r = consteig::sin(u);
            cn_r = consteig::cos(u);
            dn_r = static_cast<T>(1);
            return;
        }

        T k_seq[32]{};
        int nsteps = 0;
        T   scale  = static_cast<T>(1);
        T   kn     = k;

        while (kn > static_cast<T>(1e-12) && nsteps < 32)
        {
            T knp           = consteig::sqrt(static_cast<T>(1) - kn * kn);
            T kn1           = (static_cast<T>(1) - knp) /
                              (static_cast<T>(1) + knp);
            k_seq[nsteps++] = kn1;
            scale *= (static_cast<T>(1) + kn1);
            kn = kn1;
        }

        T phi = u / scale;
        T sn  = consteig::sin(phi);
        T cn  = consteig::cos(phi);
        T dn  = static_cast<T>(1);

        for (int i = nsteps - 1; i >= 0; --i)
        {
            T ki1   = k_seq[i];
            T denom = static_cast<T>(1) + ki1 * sn * sn;
            T sn_n  = (static_cast<T>(1) + ki1) * sn / denom;
            T cn_n  = cn * dn / denom;
            T dn_n  = (static_cast<T>(1) - ki1 * sn * sn) / denom;
            sn = sn_n;
            cn = cn_n;
            dn = dn_n;
        }
        sn_r = sn;
        cn_r = cn;
        dn_r = dn;
    }

    // =========================================================================
    // Compute v0 via bisection on sn(u, k1p) = 1/sqrt(1+ep^2).
    //
    //   v0 = asn(1/sqrt(1+ep^2), k1p) / (N * K1)
    //
    // ep  - passband ripple factor
    // k1p - complement of selectivity k1 = ep/es
    // K1p - K(k1p), upper bound for bisection
    // K1  - K(k1), used in final normalization
    // =========================================================================
    static constexpr T compute_v0(T ep, T k1p, T K1p, T K1)
    {
        const T y_target =
            static_cast<T>(1) /
            consteig::sqrt(static_cast<T>(1) + ep * ep);

        T u_lo = static_cast<T>(0);
        T u_hi = K1p;

        for (int iter = 0; iter < 64; ++iter)
        {
            T um = (u_lo + u_hi) / static_cast<T>(2);
            T sn_m{}, cn_m{}, dn_m{};
            jac_scd(um, k1p, sn_m, cn_m, dn_m);
            if (sn_m < y_target)
                u_lo = um;
            else
                u_hi = um;
        }

        return (u_lo + u_hi) /
               (static_cast<T>(2) * static_cast<T>(N) * K1);
    }

    // =========================================================================
    // In-place multiply polynomial poly (ascending order, currently degree
    // deg-1) by (s - root), bringing it to degree deg.
    //
    // Traverse j from deg downto 1 to avoid overwriting values still needed,
    // then handle j=0 separately.
    // =========================================================================
    static constexpr void poly_mul_root(Cx (&poly)[N + 1u],
                                        consteig::Size   deg,
                                        Cx               root)
    {
        for (consteig::Size j = deg; j > 0u; --j)
            poly[j] = poly[j - 1u] - root * poly[j];
        poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - root * poly[0];
    }

    // =========================================================================
    // Low-pass elliptic transfer function.
    //
    // Produces b[] and a[] in descending power order for passband edge wc.
    //
    // Poles and zeros of the normalized prototype (passband edge = 1 rad/s):
    //
    //   Zero pairs (imaginary axis): +-j*Omega_z_m
    //     Omega_z_m = 1 / (k * cd(u_m * K, k)),  u_m = (2m-1)/N
    //
    //   Pole pairs (LHP): p_m = j * cd((u_m - j*v0) * K, k)
    //     Evaluated via Jacobi addition theorem (all real arithmetic):
    //       Let sn_u,cn_u,dn_u = Jacobi(u_m*K, k)
    //           sn_v,cn_v,dn_v = Jacobi(v0*K,  k')
    //       A = cn_u*cn_v, B = sn_u*dn_u*sn_v*dn_v,
    //       C = dn_u*dn_v*cn_v, E = k^2*sn_u*cn_u*sn_v
    //       cd_re = (A*C + B*E)/(C^2+E^2),  cd_im = (B*C - A*E)/(C^2+E^2)
    //       p_m = -cd_im + j*cd_re
    //
    //   Real pole (odd N): p_0 = -sn(v0*K, k') / cn(v0*K, k')
    //
    // Gain: H(0) = 1 for odd N;  H(0) = Gp = 1/sqrt(1+ep^2) for even N.
    // =========================================================================
    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], LowPass)
    {
        // Ripple factors.
        const T ep =
            consteig::sqrt(from_db10(ripple_db) - static_cast<T>(1));
        const T es =
            consteig::sqrt(from_db10(attenuation_db) - static_cast<T>(1));

        const T k1  = ep / es;
        const T k1p = consteig::sqrt(static_cast<T>(1) - k1 * k1);
        const T K1  = elliptic_K(k1);
        const T K1p = elliptic_K(k1p);

        // Degree equation: bisect k in (0,1) until K(k')/K(k) = K1p/(N*K1).
        // Derivation: N = K(k)*K1p / (K(kp)*K1)  =>  K(kp)/K(k) = K1p/(N*K1).
        // K'/K is monotone decreasing, so ratio > target => k_lo = km.
        const T target = K1p / (static_cast<T>(N) * K1);
        T       k_lo   = static_cast<T>(1e-6);
        T       k_hi   = static_cast<T>(1) - static_cast<T>(1e-6);
        for (int iter = 0; iter < 64; ++iter)
        {
            T km  = (k_lo + k_hi) / static_cast<T>(2);
            T kmp = consteig::sqrt(static_cast<T>(1) - km * km);
            if (elliptic_K(kmp) / elliptic_K(km) > target)
                k_lo = km;
            else
                k_hi = km;
        }
        const T k    = (k_lo + k_hi) / static_cast<T>(2);
        const T kp   = consteig::sqrt(static_cast<T>(1) - k * k);
        const T Kval = elliptic_K(k);

        // Pole-shift parameter v0.
        const T v0 = compute_v0(ep, k1p, K1p, K1);

        // Jacobi at (v0*Kval, kp): shared across all pole computations.
        T sn_v{}, cn_v{}, dn_v{};
        jac_scd(v0 * Kval, kp, sn_v, cn_v, dn_v);

        // Ascending-order polynomials: poly[i] = coefficient of s^i.
        Cx poly_a[N + 1u]{};
        Cx poly_b[N + 1u]{};
        poly_a[0] = Cx{static_cast<T>(1), static_cast<T>(0)};
        poly_b[0] = Cx{static_cast<T>(1), static_cast<T>(0)};

        consteig::Size deg_a = 0u;
        consteig::Size deg_b = 0u;

        for (consteig::Size m = 1u; m <= M; ++m)
        {
            const T u_m = static_cast<T>(2u * m - 1u) /
                          static_cast<T>(N);

            T sn_u{}, cn_u{}, dn_u{};
            jac_scd(u_m * Kval, k, sn_u, cn_u, dn_u);

            // Zero pair at +-j*omega_z.
            const T omega_z =
                static_cast<T>(1) / (k * (cn_u / dn_u));

            ++deg_b;
            poly_mul_root(poly_b, deg_b,
                          Cx{static_cast<T>(0), omega_z});
            ++deg_b;
            poly_mul_root(poly_b, deg_b,
                          Cx{static_cast<T>(0), -omega_z});

            // Complex pole via addition theorem: p_m = j * cd((u_m - j*v0)*K, k).
            const T A  = cn_u * cn_v;
            const T B  = sn_u * dn_u * sn_v * dn_v;
            const T C  = dn_u * dn_v * cn_v;
            const T E  = k * k * sn_u * cn_u * sn_v;
            const T D2 = C * C + E * E;

            const T cd_re = (A * C + B * E) / D2;
            const T cd_im = (B * C - A * E) / D2;

            // p_m = j*(cd_re + j*cd_im) = -cd_im + j*cd_re
            ++deg_a;
            poly_mul_root(poly_a, deg_a,
                          Cx{-cd_im, cd_re});
            ++deg_a;
            poly_mul_root(poly_a, deg_a,
                          Cx{-cd_im, -cd_re}); // conjugate pole

        }

        // Real pole for odd N.
        if (N % 2u == 1u)
        {
            ++deg_a;
            poly_mul_root(poly_a, deg_a,
                          Cx{-sn_v / cn_v, static_cast<T>(0)});
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
        const T Gp   = static_cast<T>(1) /
                       consteig::sqrt(static_cast<T>(1) + ep * ep);
        const T H0   = (N % 2u == 1u) ? static_cast<T>(1) : Gp;
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
        elliptic_tf(static_cast<T>(1), ripple_db, attenuation_db,
                    b_lp, a_lp, LowPass{});

        for (consteig::Size j = 0u; j <= N; ++j)
        {
            const T sc = consteig::pow(wc, static_cast<int>(j));
            a[j]       = a_lp[N - j] * sc;
            b[j]       = b_lp[N - j] * sc;
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_ELLIPTIC_HPP
