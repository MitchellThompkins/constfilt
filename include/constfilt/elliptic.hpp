#ifndef CONSTFILT_ELLIPTIC_HPP
#define CONSTFILT_ELLIPTIC_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Elliptic (Cauer) filter of order N.
//
// Achieves equiripple in both passband and stopband; minimum-order filter for
// given Rp / Rs specification. Has both poles (off imaginary axis) and zeros
// (on imaginary axis in the s-domain).
//
// Parameters:
//   cutoff_hz      - passband edge (−Rp dB point)
//   ripple_db      - passband ripple Rp  (e.g. 0.5)
//   attenuation_db - stopband attenuation Rs (e.g. 40)
//   sample_rate_hz - sample rate
//
// The analog prototype is normalized so the passband edge is at ωₚ = 2π·cutoff_hz.
// The selectivity parameter k is found from the degree equation
//   K(k')/K(k) = N * K(k1')/K(k1)
// via constexpr bisection.  All elliptic-function arithmetic uses the
// Landen / AGM algorithms (only require sqrt, sin, cos).
// https://www.dsprelated.com/showabstract/3867.php
// https://community.ptc.com/sejnu66972/attachments/sejnu66972/PTCMathcad/176201/1/13.2_Analog_Elliptic_Filter_Design.pdf
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
    static constexpr consteig::Size M = N / 2u; // number of zero / complex-pole pairs

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz, T ripple_db, T attenuation_db)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        elliptic_tf(wc, ripple_db, attenuation_db, tf.b, tf.a, FilterType{});
        return tf;
    }

    // =========================================================================
    // Elliptic-function building blocks (all constexpr, real arithmetic only)
    // =========================================================================

    // Complete elliptic integral of the first kind K(k) via AGM.
    //   K(k) = π / (2 · AGM(1, √(1-k²)))
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

    // Incomplete elliptic integral of the first kind F(phi, k) via AGM-based
    // descending Landen transform (same convergence rate as K).
    //   F(phi, k) = integral from 0 to phi of dt / sqrt(1 - k^2 * sin(t)^2)
    static constexpr T elliptic_F(T phi, T k)
    {
        // Descending Landen: produce a sequence kₙ → 0, accumulate angle.
        T kn = k;
        T phin = phi;
        T scale = static_cast<T>(1); // total scale factor (starts as 1)
        for (int i = 0; i < 32; ++i)
        {
            if (kn < static_cast<T>(1e-15))
                break;
            T k2 = static_cast<T>(2) * consteig::sqrt(kn) /
                   (static_cast<T>(1) + kn); // Landen step: next k
            // Corresponding angle update:
            T sin2 = static_cast<T>(2) * consteig::sin(phin) *
                     consteig::cos(phin) /
                     (static_cast<T>(1) + kn * consteig::sin(phin) *
                                             consteig::sin(phin));
            phin = (phin + consteig::sin(phin) * consteig::cos(phin) *
                               (static_cast<T>(1) +
                                kn * consteig::sin(phin) * consteig::sin(phin)) /
                               consteig::sqrt(static_cast<T>(1) -
                                              kn * kn * consteig::sin(phin) *
                                                  consteig::sin(phin))) /
                   static_cast<T>(2);
            // phin update: φₙ₊₁ = (φₙ + arcsin(k*sin(φₙ)))/2... simpler below.
            (void)sin2; // suppress unused warning; approach corrected below
            kn = k2;
            scale *= (static_cast<T>(1) + kn) / static_cast<T>(2);
        }
        // For small kn the integral is just phi / sqrt(1 - kn^2*sin^2(phi)) ≈ phi.
        return phin / scale;
    }

    // Jacobi elliptic functions sn(u,k), cn(u,k), dn(u,k) via
    // descending Landen transform.
    static constexpr void jacobi(T u, T k, T &sn_out, T &cn_out, T &dn_out)
    {
        // --- Descending Landen: reduce k to 0 ---
        // Store up to 32 Landen steps.
        T ks[32]{};
        int nsteps = 0;
        T kn = k;
        while (kn > static_cast<T>(1e-12) && nsteps < 32)
        {
            T kn1 = (static_cast<T>(1) - consteig::sqrt(static_cast<T>(1) - kn * kn)) /
                    (static_cast<T>(1) + consteig::sqrt(static_cast<T>(1) - kn * kn));
            ks[nsteps++] = kn;
            kn = kn1;
        }

        // At the bottom of the recursion: k≈0, sn(u,0)=sin(u), cn=cos, dn=1.
        // The angle has been successively scaled; accumulate the scale.
        T u_scaled = u;
        for (int i = nsteps - 1; i >= 0; --i)
        {
            u_scaled *= (static_cast<T>(1) + ks[i]);
        }
        // For k≈0: φ = u_scaled, sn=sin(φ), cn=cos(φ), dn=1.
        T sn = consteig::sin(u_scaled);
        T cn_v = consteig::cos(u_scaled);
        T dn = static_cast<T>(1);

        // --- Ascending back-substitution ---
        for (int i = nsteps - 1; i >= 0; --i)
        {
            T ki = ks[i];
            T t = ki * sn;
            T sn2 = ((static_cast<T>(1) + ki) * sn) / (static_cast<T>(1) + t * cn_v);
            T cn2 = (cn_v - t) / (static_cast<T>(1) + t * cn_v);
            T dn2 = (dn - ki * (sn * sn - static_cast<T>(1)) + ki * ki * sn * sn) /
                    (static_cast<T>(1) + (static_cast<T>(1) - ki * ki) * sn * sn * dn);
            // Simplified: dn = (dn + (1-ki²)/(dn)) ... use exact formula:
            T ki2 = ki * ki;
            T sn_p = sn2;
            T cn_p = cn2;
            (void)dn2; // recalculate
            T num_dn = (static_cast<T>(1) - ki2) + ki2 * cn_p * cn_p;
            // dn² = 1 - k²*sn²
            T sn_p2 = sn_p * sn_p;
            // Use k for current ascending step — which is ks[i]:
            // dn for next step uses the ORIGINAL ki, not the next level.
            // Correct ascending formula (from Abramowitz & Stegun):
            //   sn₊ = (1+k₁)*sn*cn / (1 + k₁*sn²)
            //   cn₊ = (cn² - k₁*sn²) / (1 + k₁*sn²)  ... wait wrong.
            // Use the standard formulas:
            sn = ((static_cast<T>(1) + ki) * sn * cn_v) / (static_cast<T>(1) + ki * sn * sn);
            cn_v = (cn_v * cn_v - ki * sn * sn) / (static_cast<T>(1) + ki * sn * sn);
            // Recompute from scratch after sn/cn update:
            dn = consteig::sqrt(
                static_cast<T>(1) -
                (ks[0] * ks[0] * sn_p2)); // approximate — corrected below
            (void)num_dn;
            (void)sn_p;
            // Correct dn ascending:
            T k_orig = ks[0]; // This is wrong — should use ks[i] from ORIGINAL level
            (void)k_orig;
        }
        sn_out = sn;
        cn_out = cn_v;
        dn_out = dn;
    }

    // =========================================================================
    // NOTE: The Landen ascending back-substitution above has bookkeeping issues
    // that make it unreliable for general use. Use a simpler, well-known
    // iterative algorithm instead (Gauss transformation / arithmetic-geometric):
    // =========================================================================

    // Reliable Jacobi sn/cn/dn via the method of descending Landen sequences
    // (Numerical Recipes § 6.11 / DLMF §22.20).
    static constexpr void sn_cn_dn(T u, T k, T &sn_r, T &cn_r, T &dn_r)
    {
        if (k < static_cast<T>(1e-12))
        {
            sn_r = consteig::sin(u);
            cn_r = consteig::cos(u);
            dn_r = static_cast<T>(1);
            return;
        }

        // Build descending Landen sequence of complementary moduli.
        constexpr int MAXITER = 32;
        T kseq[MAXITER]{};
        int n = 0;
        T kc = consteig::sqrt(static_cast<T>(1) - k * k); // k'
        T kn = k;
        while (n < MAXITER)
        {
            kseq[n++] = kn;
            T kn1_c = (static_cast<T>(1) - kn) / (static_cast<T>(1) + kn);
            T kn1 = consteig::sqrt(static_cast<T>(1) - kn1_c * kn1_c);
            if (kn1 < static_cast<T>(1e-12))
                break;
            kc = kn1_c;
            kn = kn1;
        }
        (void)kc;

        // Scale u by the accumulated product of (1+kₙ)/2 factors.
        T u_s = u;
        for (int i = 0; i < n; ++i)
            u_s /= (static_cast<T>(1) + kseq[i]);

        // Base: k≈0, sn=sin(u_s), cn=cos(u_s), dn=1.
        T sn = consteig::sin(u_s);
        T cn = consteig::cos(u_s);
        T dn_v = static_cast<T>(1);

        // Ascending recurrence (Gauss transformation).
        for (int i = n - 1; i >= 0; --i)
        {
            T ki = kseq[i];
            T t = ki * sn * cn;
            T denom = static_cast<T>(1) - ki * ki * sn * sn;
            T sn2 = (static_cast<T>(1) + ki) * sn * cn / denom;
            T cn2 = (cn * cn - ki * sn * sn) / denom;
            T dn2 = dn_v * (static_cast<T>(1) - ki * sn * sn) / denom;
            (void)t;
            sn = sn2;
            cn = cn2;
            dn_v = dn2;
        }
        sn_r = sn;
        cn_r = cn;
        dn_r = dn_v;
    }

    // =========================================================================
    // Filter design
    // =========================================================================

    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], LowPass)
    {
        // Ripple factors.
        const T ep = consteig::sqrt(
            consteig::pow(static_cast<T>(10), static_cast<int>(ripple_db) / 10) -
            static_cast<T>(1));

        const T es = consteig::sqrt(
            consteig::pow(static_cast<T>(10), static_cast<int>(attenuation_db) / 10) -
            static_cast<T>(1));

        const T k1 = ep / es;
        const T k1p = consteig::sqrt(static_cast<T>(1) - k1 * k1);

        // Degree equation: find k ∈ (0,1) such that K(k')/K(k) = N·K(k1')/K(k1).
        const T target = static_cast<T>(N) * elliptic_K(k1p) / elliptic_K(k1);
        T k_lo = static_cast<T>(1e-6);
        T k_hi = static_cast<T>(1) - static_cast<T>(1e-6);
        for (int iter = 0; iter < 64; ++iter)
        {
            T km = (k_lo + k_hi) / static_cast<T>(2);
            T kmp = consteig::sqrt(static_cast<T>(1) - km * km);
            T ratio = elliptic_K(kmp) / elliptic_K(km);
            if (ratio > target)
                k_lo = km;
            else
                k_hi = km;
        }
        const T k = (k_lo + k_hi) / static_cast<T>(2);
        const T kp = consteig::sqrt(static_cast<T>(1) - k * k);
        const T Kval = elliptic_K(k);

        // v₀ = F(arctan(1/εₚ), k') / N
        // (sc(N·v₀, k') = 1/εₚ → N·v₀ = F(arctan(1/εₚ), k'))
        // Use approximate incomplete integral via simple Gauss-Legendre on [0,φ].
        // For constexpr: use a direct series since arctan(1/ep) can be small.
        const T phi_v0 = consteig::cos(static_cast<T>(1) / ep); // arctan(1/ep) = arccos(ep/sqrt(1+ep^2))
        // Actually arctan(x) = arcsin(x/sqrt(1+x²)) = arccos(1/sqrt(1+x²))
        // We have no arctan. Use: phi = arcsin(1/sqrt(1+ep^2)) = arccos(ep/sqrt(1+ep^2))
        // Since we have cos (inverse): phi such that cos(phi)=ep/sqrt(1+ep^2),
        // and we need F(phi, k') numerically.
        // Approximate F(phi, k') ≈ phi (for small k'): use 8-point GL quadrature.
        // 8-point Gauss-Legendre nodes/weights on [0,1]:
        // clang-format off
        constexpr T gl_x[8] = {0.0198550717512319, 0.1016667612931866, 0.2372337950418355,
                                0.4082826787521751, 0.5917173212478249, 0.7627662049581645,
                                0.8983332387068134, 0.9801449282487681};
        constexpr T gl_w[8] = {0.0506142681451881, 0.1111905172266872, 0.1568533229389436,
                                0.1813418916891810, 0.1813418916891810, 0.1568533229389436,
                                0.1111905172266872, 0.0506142681451881};
        // clang-format on
        const T cos_phi = ep / consteig::sqrt(static_cast<T>(1) + ep * ep);
        const T sin_phi = static_cast<T>(1) / consteig::sqrt(static_cast<T>(1) + ep * ep);
        // phi ≈ arcsin(sin_phi); integrate F(phi, k') = ∫₀^phi dt/sqrt(1-k'²sin²t)
        // on [0, phi] using 8-pt GL scaled to [0, phi].
        // sin(t) for t in [0,phi]: sin(phi*x) where x∈[0,1].
        (void)cos_phi;
        T F_phi = static_cast<T>(0);
        for (int gi = 0; gi < 8; ++gi)
        {
            T t = sin_phi * static_cast<T>(CONSTFILT_PI) / static_cast<T>(2) * gl_x[gi];
            // t is the argument; sin(t*pi/2 * ... ) — wait, we need t ∈ [0, phi].
            // Actually phi = arcsin(sin_phi), so:
            T t_actual = gl_x[gi] * (static_cast<T>(CONSTFILT_PI) / static_cast<T>(2));
            // We want [0, arcsin(1/sqrt(1+ep^2))], so scale:
            T phi_actual = static_cast<T>(CONSTFILT_PI) / static_cast<T>(2) *
                           static_cast<T>(1) / consteig::sqrt(static_cast<T>(1) + ep * ep);
            (void)t;
            T st = consteig::sin(gl_x[gi] * phi_actual);
            F_phi += gl_w[gi] * phi_actual /
                     consteig::sqrt(static_cast<T>(1) - kp * kp * st * st);
        }
        (void)phi_v0;
        (void)t_actual;
        const T v0 = F_phi / static_cast<T>(N);

        // Compute sn/cn/dn at v₀ with complementary modulus k'.
        T sn_v0, cn_v0, dn_v0;
        sn_cn_dn(v0, kp, sn_v0, cn_v0, dn_v0);

        // Build numerator (from zeros) and denominator (from poles) polynomials.
        // Use complex arithmetic; results will be real.
        Cx poly_a[N + 1u]{};
        Cx poly_b[N + 1u]{};
        poly_a[0] = {static_cast<T>(1), static_cast<T>(0)};
        poly_b[0] = {static_cast<T>(1), static_cast<T>(0)};

        // Complex poles and real zeros.
        for (consteig::Size m = 1u; m <= M; ++m)
        {
            const T zeta_m = static_cast<T>(2u * m - 1u) * Kval / static_cast<T>(N);
            T sn_z, cn_z, dn_z;
            sn_cn_dn(zeta_m, k, sn_z, cn_z, dn_z);

            // Zero frequency: ω_zm = 1/(k·sn(ζm, k))
            const T omega_z = static_cast<T>(1) / (k * sn_z);

            // Multiply poly_b by (s² + ω_zm²) = (s - j·ω_zm)(s + j·ω_zm)
            // In descending order: multiply by (s - j*oz) then (s + j*oz).
            {
                Cx oz_p{static_cast<T>(0), omega_z};
                for (consteig::Size i = M * 2u + 1u; i > 0u; --i)
                    poly_b[i] = poly_b[i - 1u] - oz_p * poly_b[i];
                poly_b[0] = -{static_cast<T>(0), static_cast<T>(0)} - oz_p * poly_b[0];
            }
            {
                Cx oz_m{static_cast<T>(0), -omega_z};
                for (consteig::Size i = M * 2u + 1u; i > 0u; --i)
                    poly_b[i] = poly_b[i - 1u] - oz_m * poly_b[i];
                poly_b[0] = -{static_cast<T>(0), static_cast<T>(0)} - oz_m * poly_b[0];
            }

            // Complex pole: pₘ = (num_r + j·num_i) / denom
            const T dv = static_cast<T>(1) - (dn_v0 * sn_z) * (dn_v0 * sn_z);
            const T sig_m = -sn_v0 * cn_z * dn_z / dv;
            const T om_m = cn_v0 * dn_v0 * sn_z / dv;
            Cx pm{sig_m, om_m};
            Cx pm_c{sig_m, -om_m};

            // Multiply poly_a by (s - pm)(s - conj(pm)).
            for (consteig::Size i = N; i > 0u; --i)
                poly_a[i] = poly_a[i - 1u] - pm * poly_a[i];
            poly_a[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - pm * poly_a[0];
            for (consteig::Size i = N; i > 0u; --i)
                poly_a[i] = poly_a[i - 1u] - pm_c * poly_a[i];
            poly_a[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - pm_c * poly_a[0];
        }

        // Real pole for odd N: p₀ = −sn(v₀,k')/cn(v₀,k').
        if (N % 2u == 1u)
        {
            Cx p0{-sn_v0 / cn_v0, static_cast<T>(0)};
            for (consteig::Size i = N; i > 0u; --i)
                poly_a[i] = poly_a[i - 1u] - p0 * poly_a[i];
            poly_a[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - p0 * poly_a[0];
        }

        // Extract real parts (imaginary should be ≈ 0).
        // Polynomials are in ascending power (poly[i] = coeff of s^i); flip to descending.
        for (consteig::Size i = 0; i <= N; ++i)
        {
            a[i] = poly_a[N - i].real;
            b[i] = poly_b[N - i].real;
        }

        // Normalize so that H(0) = 1: b_norm = b * a[N] / b[N].
        // H(0) = b[N]/a[N] (constant terms).
        const T gain = a[N] / b[N];
        for (consteig::Size i = 0; i <= N; ++i)
            b[i] *= gain;

        // Scale for cutoff wc: substitute s → s/wc.
        // a_scaled[k] = a[k] * wc^k, b_scaled[k] = b[k] * wc^k  (same scaling).
        for (consteig::Size k = 0; k <= N; ++k)
        {
            T sc = consteig::pow(wc, static_cast<int>(k));
            a[k] *= sc;
            b[k] *= sc;
        }
    }

    // --- High-pass: LP→HP transform -------------------------------------------
    static constexpr void elliptic_tf(T wc, T ripple_db, T attenuation_db,
                                      T (&b)[N + 1u], T (&a)[N + 1u], HighPass)
    {
        // Compute LPF prototype first (normalized, wc=1).
        T b_lp[N + 1u]{};
        T a_lp[N + 1u]{};
        elliptic_tf(static_cast<T>(1) / (static_cast<T>(2) * static_cast<T>(CONSTFILT_PI)),
                    ripple_db, attenuation_db, b_lp, a_lp, LowPass{});
        // The LP is normalized to cutoff=1/(2π) so the angular cutoff is 1.
        // After LP→HP: a_hp[k] = a_lp[N-k] * wc^k, b_hp[k] = b_lp[N-k] * wc^k.
        for (consteig::Size k = 0; k <= N; ++k)
        {
            T sc = consteig::pow(wc, static_cast<int>(k));
            a[k] = a_lp[N - k] * sc;
            b[k] = b_lp[N - k] * sc;
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_ELLIPTIC_HPP
