#ifndef CONSTFILT_LEGENDRE_PAPOULIS_HPP
#define CONSTFILT_LEGENDRE_PAPOULIS_HPP

#include "analog_filter.hpp"
#include "constfilt_options.hpp"
#include <consteig/consteig.hpp>

namespace constfilt
{

// Legendre–Papoulis (optimal monotone) filter of order N.
//
// Maximizes the transition-band slope at ω = 1 while keeping the magnitude
// response strictly monotone — no passband or stopband ripple.
//
// The squared magnitude is  |H(jω)|² = 1 / (1 + f(ω²))  where
//   m = ⌊(N+1)/2⌋
//   c = 2m − 1
//   g(x)  = ∫₀ˣ [P_{m-1}(2t−1)]² dt      (degree 2m−1 polynomial, g(0)=0)
//   f(x)  = c · g(x)                        (odd  N = 2m−1, so f(1) = 1)
//   f(x)  = x · c · g(x)                   (even N = 2m,   so f(1) = 1)
//
// The denominator D_N(s) is the degree-N polynomial whose roots are the
// left-half-plane roots of D_N(s)·D_N(−s) = 1 + f(−s²).
// These roots are found via companion-matrix eigendecomposition (consteig).
template <typename T, consteig::Size N, typename Method = ZOH,
          typename FilterType = LowPass>
class LegendrePapoulis : public AnalogFilter<T, N, Method>
{
    static_assert(N >= 1u, "Legendre-Papoulis order must be at least 1");

  public:
    constexpr LegendrePapoulis(T cutoff_hz, T sample_rate_hz)
        : AnalogFilter<T, N, Method>(compute_continuous_tf(cutoff_hz),
                                     sample_rate_hz)
    {
    }

  private:
    using Cx = consteig::Complex<T>;
    // m = ⌊(N+1)/2⌋
    static constexpr consteig::Size M_IDX = (N + 1u) / 2u;

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_continuous_tf(
        T cutoff_hz)
    {
        const T wc = static_cast<T>(2.0 * CONSTFILT_PI) * cutoff_hz;
        TransferFunction<T, N + 1u, N + 1u> tf{};
        continuous_tf(wc, tf.b, tf.a, FilterType{});
        return tf;
    }

    // =========================================================================
    // Polynomial helpers (all ascending-power order: coeff[k] = coefficient of x^k)
    // =========================================================================

    // Compute coefficients of P_{m-1}(2x-1) of degree m-1 into p[0..m-1].
    // Uses Bonnet's recurrence on P_k(y) with y = 2x-1, converted to poly in x.
    static constexpr void legendre_shifted(T (&p)[M_IDX])
    {
        // P_0(y) = 1 → poly [1].
        // P_1(y) = y = 2x-1 → poly [-1, 2].
        // P_k(y) = ((2k-1)*y*P_{k-1}(y) - (k-1)*P_{k-2}(y)) / k
        //        = ((2k-1)*(2x-1)*P_{k-1} - (k-1)*P_{k-2}) / k   [poly multiplication]

        if (M_IDX == 0u)
            return;

        // Previous two polynomials stored as full M_IDX-length arrays.
        T prev2[M_IDX]{}; // P_{k-2}
        T prev1[M_IDX]{}; // P_{k-1}
        T cur[M_IDX]{};

        // P_0 = 1
        prev1[0] = static_cast<T>(1);

        if (M_IDX == 1u)
        {
            p[0] = prev1[0];
            return;
        }

        // P_1(2x-1) = 2x-1 → [-1, 2]
        prev2[0] = prev1[0];
        prev1[0] = static_cast<T>(-1);
        if (M_IDX > 1u)
            prev1[1] = static_cast<T>(2);

        if (M_IDX == 2u)
        {
            for (consteig::Size i = 0; i < M_IDX; ++i)
                p[i] = prev1[i];
            return;
        }

        // Recurrence for k = 2..M_IDX-1
        for (consteig::Size k = 2u; k < M_IDX; ++k)
        {
            for (consteig::Size i = 0; i < M_IDX; ++i)
                cur[i] = static_cast<T>(0);

            // (2k-1)*(2x-1)*P_{k-1}:
            //   First multiply P_{k-1} by (2x-1):
            //     (2x-1)*p(x) = 2*x*p(x) - p(x)
            //   coeff[j] of result = 2*prev1[j-1] - prev1[j]
            T tmp[M_IDX]{};
            for (consteig::Size j = 0; j < M_IDX; ++j)
            {
                tmp[j] = -prev1[j];
                if (j > 0u)
                    tmp[j] += static_cast<T>(2) * prev1[j - 1u];
            }
            const T c1 = static_cast<T>(2u * k - 1u);
            for (consteig::Size j = 0; j < M_IDX; ++j)
                cur[j] += c1 * tmp[j];

            // - (k-1)*P_{k-2}
            const T c2 = -static_cast<T>(k - 1u);
            for (consteig::Size j = 0; j < M_IDX; ++j)
                cur[j] += c2 * prev2[j];

            // Divide by k
            const T inv_k = static_cast<T>(1) / static_cast<T>(k);
            for (consteig::Size j = 0; j < M_IDX; ++j)
                cur[j] *= inv_k;

            // Advance
            for (consteig::Size j = 0; j < M_IDX; ++j)
            {
                prev2[j] = prev1[j];
                prev1[j] = cur[j];
            }
        }

        for (consteig::Size i = 0; i < M_IDX; ++i)
            p[i] = prev1[i];
    }

    // Square polynomial p (degree m-1) → result of degree 2m-2.
    // result has 2*M_IDX-1 coefficients.
    static constexpr void poly_square(const T (&p)[M_IDX],
                                      T (&result)[2u * M_IDX - 1u])
    {
        for (consteig::Size i = 0; i < 2u * M_IDX - 1u; ++i)
            result[i] = static_cast<T>(0);
        for (consteig::Size i = 0; i < M_IDX; ++i)
            for (consteig::Size j = 0; j < M_IDX; ++j)
                result[i + j] += p[i] * p[j];
    }

    // Integrate p2 (degree 2m-2, 2m-1 coefficients) from 0 to x:
    // g[k] = p2[k-1]/k for k=1..2m-1, g[0]=0.  Result: 2m coefficients.
    static constexpr void poly_integrate(const T (&p2)[2u * M_IDX - 1u],
                                         T (&g)[2u * M_IDX])
    {
        g[0] = static_cast<T>(0);
        for (consteig::Size k = 1u; k <= 2u * M_IDX - 1u; ++k)
            g[k] = p2[k - 1u] / static_cast<T>(k);
    }

    // Build the G polynomial (ascending, N+1 coefficients) from g.
    // G(u) = 1 + f(-u²) where:
    //   odd  N: f(x) = c*g(x)      → G(u) = 1 + c*g(-u)
    //   even N: f(x) = x*c*g(x)    → G(u) = 1 + (-u)*c*g(-u) = 1 - u*c*g(-u)
    static constexpr void build_G(const T (&g)[2u * M_IDX], T (&G)[N + 1u])
    {
        const T c = static_cast<T>(2u * M_IDX - 1u); // 2m-1
        for (consteig::Size i = 0; i <= N; ++i)
            G[i] = static_cast<T>(0);
        G[0] = static_cast<T>(1);

        if (N % 2u == 1u)
        {
            // Odd N = 2m-1: G(u) = 1 + c*g(-u)
            // g(-u) = Σ g[k]*(-u)^k = Σ g[k]*(-1)^k*u^k
            for (consteig::Size k = 1u; k <= N; ++k)
            {
                T sign = (k % 2u == 0u) ? static_cast<T>(1) : static_cast<T>(-1);
                G[k] += c * g[k] * sign;
            }
        }
        else
        {
            // Even N = 2m: G(u) = 1 - u*c*g(-u)
            // -u*g(-u): shift g(-u) up by 1 power and negate.
            // G[k] = -c * (-1)^(k-1) * g[k-1] for k=1..N
            for (consteig::Size k = 1u; k <= N; ++k)
            {
                T sign_gm1 = ((k - 1u) % 2u == 0u) ? static_cast<T>(1)
                                                     : static_cast<T>(-1);
                G[k] += -c * g[k - 1u] * sign_gm1;
            }
        }
    }

    // Complex sqrt (same as Gaussian helper).
    static constexpr Cx cx_sqrt(const Cx &z)
    {
        T r = consteig::abs(z);
        T re = consteig::sqrt((r + z.real) / static_cast<T>(2));
        T im = consteig::sqrt((r - z.real) / static_cast<T>(2));
        if (z.imag < static_cast<T>(0))
            im = -im;
        return {re, im};
    }

    // Compute the normalized Legendre-Papoulis denominator polynomial
    // (descending power, monic, −3 dB at ω = 1).
    static constexpr void lp_poly_coeffs(T (&result)[N + 1u])
    {
        // --- 1. Build Legendre poly P_{m-1}(2x-1) ----------------------------
        T p[M_IDX]{};
        legendre_shifted(p);

        // --- 2. Square it -----------------------------------------------------
        T p2[2u * M_IDX - 1u]{};
        poly_square(p, p2);

        // --- 3. Integrate from 0 to x → g(x) ---------------------------------
        T g[2u * M_IDX]{};
        poly_integrate(p2, g);

        // --- 4. Build G(u) = 1 + f(-u) (degree N) ----------------------------
        T G[N + 1u]{};
        build_G(g, G);

        // --- 5. Companion matrix for G, eigenvalues → roots u_k ---------------
        // Monic form: divide by leading coefficient G[N].
        consteig::Matrix<T, N, N> C{};
        const T gN = G[N];
        for (consteig::Size i = 0; i + 1u < N; ++i)
            C(i + 1u, i) = static_cast<T>(1);
        for (consteig::Size k = 0; k < N; ++k)
            C(k, N - 1u) = -G[k] / gN;

        auto evals = consteig::eigenvalues(C);

        // --- 6. Poles s_k = −√u_k  (LHP selection) ---------------------------
        Cx poles[N]{};
        for (consteig::Size k = 0; k < N; ++k)
        {
            Cx sq = cx_sqrt(evals(k, 0));
            if (sq.real >= static_cast<T>(0))
                sq = {-sq.real, -sq.imag};
            poles[k] = sq;
        }

        // --- 7. Build polynomial from poles (same as Butterworth) -------------
        Cx poly[N + 1u]{};
        poly[0] = Cx{static_cast<T>(1), static_cast<T>(0)};
        for (consteig::Size k = 0; k < N; ++k)
        {
            const Cx &p_k = poles[k];
            for (consteig::Size j = k + 1u; j > 0u; --j)
                poly[j] = poly[j - 1u] - p_k * poly[j];
            poly[0] = Cx{static_cast<T>(0), static_cast<T>(0)} - p_k * poly[0];
        }
        // Flip ascending → descending.
        for (consteig::Size i = 0; i <= N; ++i)
            result[i] = poly[N - i].real;
    }

    // --- Low-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        LowPass)
    {
        T p[N + 1u]{};
        lp_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[k] * consteig::pow(wc, static_cast<int>(k));
        b[N] = a[N]; // DC gain = 1
    }

    // --- High-pass
    static constexpr void continuous_tf(T wc, T (&b)[N + 1u], T (&a)[N + 1u],
                                        HighPass)
    {
        T p[N + 1u]{};
        lp_poly_coeffs(p);
        for (consteig::Size k = 0; k <= N; ++k)
            a[k] = p[N - k] * consteig::pow(wc, static_cast<int>(k));
        b[0] = static_cast<T>(1); // HF gain = 1
    }
};

} // namespace constfilt

#endif // CONSTFILT_LEGENDRE_PAPOULIS_HPP
