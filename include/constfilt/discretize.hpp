#ifndef CONSTFILT_DISCRETIZE_HPP
#define CONSTFILT_DISCRETIZE_HPP

#include "constfilt_options.hpp"
#include "vendor/consteig/consteig.hpp"

namespace constfilt
{

// --------------------------- Tag types ---------------------------------------

struct ZOH
{
};

struct MatchedZ
{
};

// --------------------------- Data structures ---------------------------------

template <typename T, consteig::Size N> struct StateSpace
{
    consteig::Matrix<T, N, N> A{};
    consteig::Matrix<T, N, 1> B{};
    consteig::Matrix<T, 1, N> C{};
    T D{};
};

template <typename T, consteig::Size NB, consteig::Size NA>
struct TransferFunction
{
    T b[NB]{};
    T a[NA]{};
};

// --------------------------- Matrix exponential ------------------------------

// matrix_exp(A) via eigendecomposition:
//   Ad = V * diag(exp(lam_i)) * V^{-1}
// where V = eigenvectors, lam_i = eigenvalues (complex).
// Real part extracted at the end (imaginary parts cancel for real A).
template <typename T, consteig::Size N>
constexpr consteig::Matrix<T, N, N> matrix_exp(
    const consteig::Matrix<T, N, N> &A)
{
    using Complex = consteig::Complex<T>;
    using ComplexMat_NN = consteig::Matrix<Complex, N, N>;
    using ComplexMat_N1 = consteig::Matrix<Complex, N, 1>;

    // 1. Eigenvalues and eigenvectors
    auto evals = consteig::eigenvalues(A);     // Matrix<Complex, N, 1>
    auto V = consteig::eigenvectors(A, evals); // Matrix<Complex, N, N>

    // 2. Invert V column-by-column via LU
    auto lu_V = consteig::lu(V);

    ComplexMat_NN V_inv{};
    for (consteig::Size col = 0; col < N; ++col)
    {
        ComplexMat_N1 e_col{};
        e_col(col, 0) = Complex{static_cast<T>(1), static_cast<T>(0)};
        auto col_vec = consteig::lu_solve(lu_V, e_col);
        for (consteig::Size row = 0; row < N; ++row)
        {
            V_inv(row, col) = col_vec(row, 0);
        }
    }

    // 3. Accumulate: matrix_exp = sum_i exp(lam_i) * v_i * w_i^T
    //    where v_i = column i of V, w_i^T = row i of V_inv
    ComplexMat_NN result_c{};
    for (consteig::Size i = 0; i < N; ++i)
    {
        Complex exp_lambda = consteig::exp(evals(i, 0));
        for (consteig::Size r = 0; r < N; ++r)
        {
            for (consteig::Size c = 0; c < N; ++c)
            {
                result_c(r, c) =
                    result_c(r, c) + exp_lambda * V(r, i) * V_inv(i, c);
            }
        }
    }

    // 4. Extract real part
    consteig::Matrix<T, N, N> result{};
    for (consteig::Size r = 0; r < N; ++r)
    {
        for (consteig::Size c = 0; c < N; ++c)
        {
            result(r, c) = result_c(r, c).real;
        }
    }
    return result;
}

// --------------------------- ZOH discretization ------------------------------

// ZOH: Ad = matrix_exp(Ac*Ts),  Bd = Ac^{-1} * (Ad - I) * Bc
// Solve  Ac * Bd = (Ad - I) * Bc  via LU.
// Cc and Dc are unchanged.
template <typename T, consteig::Size N>
constexpr StateSpace<T, N> zoh_discretize(const StateSpace<T, N> &sys_c, T Ts,
                                          ZOH /*tag*/)
{
    const auto &Ac = sys_c.A;
    const auto &Bc = sys_c.B;

    // Ac * Ts
    consteig::Matrix<T, N, N> AcTs{};
    for (consteig::Size r = 0; r < N; ++r)
    {
        for (consteig::Size c = 0; c < N; ++c)
        {
            AcTs(r, c) = Ac(r, c) * Ts;
        }
    }

    consteig::Matrix<T, N, N> Ad = matrix_exp(AcTs);

    // (Ad - I)
    consteig::Matrix<T, N, N> AdmI = Ad - consteig::eye<T, N>();

    // (Ad - I) * Bc
    consteig::Matrix<T, N, 1> rhs{};
    for (consteig::Size r = 0; r < N; ++r)
    {
        T sum = static_cast<T>(0);
        for (consteig::Size k = 0; k < N; ++k)
        {
            sum += AdmI(r, k) * Bc(k, 0);
        }
        rhs(r, 0) = sum;
    }

    // Bd = Ac^{-1} * rhs  ->  solve Ac * Bd = rhs
    auto lu_Ac = consteig::lu(Ac);
    auto Bd = consteig::lu_solve(lu_Ac, rhs);

    StateSpace<T, N> sys_d{};
    sys_d.A = Ad;
    sys_d.B = Bd;
    sys_d.C = sys_c.C;
    sys_d.D = sys_c.D;
    return sys_d;
}

// --------------------------- Characteristic polynomial -----------------------

// Fills monic characteristic polynomial of Ad:
//   [1, c_1, c_2, ..., c_N]   (N+1 coefficients)
// Uses consteig::eigenvalues to obtain lam_1..lam_N, then builds
//   (z - lam_1)(z - lam_2)...(z - lam_N) in complex arithmetic.
// Real parts are extracted at the end (imaginary parts cancel for real A).
// Note: Faddeev-LeVerrier is an acceptable alternative if consteig were
// unavailable; it computes the same coefficients via matrix traces without
// eigenvalue decomposition.
template <typename T, consteig::Size N>
constexpr void char_poly(const consteig::Matrix<T, N, N> &Ad,
                         T (&coeffs)[N + 1u])
{
    using Complex = consteig::Complex<T>;

    auto evals = consteig::eigenvalues(Ad); // Matrix<Complex, N, 1>

    // p[0..k] holds the monic polynomial of degree k after k iterations.
    Complex p[N + 1u]{};
    p[0] = Complex{static_cast<T>(1), static_cast<T>(0)};

    for (consteig::Size k = 0; k < N; ++k)
    {
        const Complex lam = evals(k, 0);
        // Multiply degree-k poly by (z - lam), working high-to-low in-place.
        p[k + 1u] = Complex{static_cast<T>(0), static_cast<T>(0)} - lam * p[k];
        for (consteig::Size i = k; i > 0u; --i)
        {
            p[i] = p[i] - lam * p[i - 1u];
        }
        // p[0] is unchanged (stays 1)
    }

    for (consteig::Size i = 0; i <= N; ++i)
    {
        coeffs[i] = p[i].real;
    }
}

// --------------------------- Markov numerator --------------------------------

// Computes numerator polynomial b from Markov parameters and denominator a.
//   h[0] = D
//   h[k] = C * A^{k-1} * B   for k = 1..N
//   b[k] = sum_{j=0}^{k} a[j] * h[k-j]
template <typename T, consteig::Size N>
constexpr void markov_numerator(const StateSpace<T, N> &sys_d,
                                const T (&a_coeffs)[N + 1u], T (&b)[N + 1u])
{
    const auto &Ad = sys_d.A;
    const auto &Bd = sys_d.B;
    const auto &Cd = sys_d.C;

    // Compute Markov parameters h[0..N]
    T h[N + 1u]{};
    h[0] = sys_d.D;

    // A_pow = A^{k-1}, starting with A^0 = I
    consteig::Matrix<T, N, N> A_pow = consteig::eye<T, N>();

    for (consteig::Size k = 1u; k <= N; ++k)
    {
        // h[k] = C * A_pow * B  (scalar)
        T val = static_cast<T>(0);
        for (consteig::Size c = 0; c < N; ++c)
        {
            // (A_pow * B)[c] = sum_col A_pow(c,col)*B(col,0)
            T AB_c = static_cast<T>(0);
            for (consteig::Size col = 0; col < N; ++col)
            {
                AB_c += A_pow(c, col) * Bd(col, 0);
            }
            val += Cd(0, c) * AB_c;
        }
        h[k] = val;

        // Advance: A_pow = A_pow * Ad
        consteig::Matrix<T, N, N> next{};
        for (consteig::Size r = 0; r < N; ++r)
        {
            for (consteig::Size cc = 0; cc < N; ++cc)
            {
                T sum = static_cast<T>(0);
                for (consteig::Size col = 0; col < N; ++col)
                {
                    sum += A_pow(r, col) * Ad(col, cc);
                }
                next(r, cc) = sum;
            }
        }
        A_pow = next;
    }

    // Convolve: b[k] = sum_{j=0}^{k} a[j] * h[k-j]
    for (consteig::Size k = 0; k <= N; ++k)
    {
        T sum = static_cast<T>(0);
        for (consteig::Size jj = 0; jj <= k; ++jj)
        {
            sum += a_coeffs[jj] * h[k - jj];
        }
        b[k] = sum;
    }
}

// --------------------------- ss_to_tf ----------------------------------------

// Full pipeline: discrete state-space -> (b, a) transfer function.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> ss_to_tf(
    const StateSpace<T, N> &sys_d)
{
    TransferFunction<T, N + 1u, N + 1u> tf{};
    char_poly(sys_d.A, tf.a);
    markov_numerator(sys_d, tf.a, tf.b);
    return tf;
}

// --------------------------- tf_to_ss ----------------------------------------

// Build controllable canonical form continuous-time state-space from
// analog transfer function coefficients in descending power order:
//
//   H(s) = (b[0]*s^N + b[1]*s^{N-1} + ... + b[N])
//          / (a[0]*s^N + a[1]*s^{N-1} + ... + a[N])
//
// Handles proper (deg b = deg a) and strictly proper (b[0]=0) cases.
// a[0] must be non-zero; the denominator is normalized to monic form
// internally.
//
// Resulting state-space (controllable canonical form):
//   A: super-diagonal = 1; last row = [-a_N, -a_{N-1}, ..., -a_1] / a_0
//   B: [0, ..., 0, 1]^T
//   C: [e_N, e_{N-1}, ..., e_1]  where e_k = b_k/a_0 - D*(a_k/a_0)
//   D: b[0] / a[0]
template <typename T, consteig::Size N>
constexpr StateSpace<T, N> tf_to_ss(const T (&b)[N + 1u], const T (&a)[N + 1u])
{
    StateSpace<T, N> sys{};

    const T inv_a0 = static_cast<T>(1) / a[0];

    sys.D = b[0] * inv_a0;

    // Numerator residual: e[k] = b[k]/a[0] - D*(a[k]/a[0])  for k=1..N
    T e[N + 1u]{};
    for (consteig::Size k = 1u; k <= N; ++k)
    {
        e[k] = b[k] * inv_a0 - sys.D * (a[k] * inv_a0);
    }

    // A: super-diagonal
    for (consteig::Size row = 0; row < N - 1u; ++row)
    {
        sys.A(row, row + 1u) = static_cast<T>(1);
    }

    // A: last row = -a[N-k]/a[0]  for k=0..N-1
    for (consteig::Size k = 0; k < N; ++k)
    {
        sys.A(N - 1u, k) = -(a[N - k] * inv_a0);
    }

    // B: last entry = 1
    sys.B(N - 1u, 0) = static_cast<T>(1);

    // C: C[0][k] = e[N-k]
    for (consteig::Size k = 0; k < N; ++k)
    {
        sys.C(0, k) = e[N - k];
    }

    return sys;
}

// --- Matched-Z discretization (TF entry point) ------------------------------

// Full matched-Z from analog transfer function coefficients.
// Maps each finite analog zero via z = exp(s*Ts), pads with zeros at z = -1
// for strictly proper systems, and matches gain at a test frequency w_c.
// Reference: Octave control pkg @tf/__c2d__.m, lines 32-66.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> matched_z_discretize_tf(
    const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T Ts, MatchedZ /*tag*/)
{
    using Complex = consteig::Complex<T>;

    // Step 1: count leading zeros in b_c; derive nz (finite analog zero count).
    consteig::Size d_b = 0;
    while (d_b <= N && b_c[d_b] == static_cast<T>(0))
    {
        ++d_b;
    }
    const consteig::Size nz = (d_b > N) ? 0u : N - d_b;

    // Step 2: continuous leading-coefficient gain.
    const T k_c = (d_b > N) ? static_cast<T>(0) : b_c[d_b] / a_c[0];

    // Step 3: find analog poles (roots of a_c) via companion matrix.
    // consteig exposes eigenvalues, not polynomial roots; a companion matrix
    // has eigenvalues equal to the polynomial roots by construction.
    consteig::Matrix<T, N, N> A_pole{};
    for (consteig::Size i = 0; i < N - 1u; ++i)
        A_pole(i + 1u, i) = static_cast<T>(1);
    for (consteig::Size i = 0; i < N; ++i)
        A_pole(i, N - 1u) = -(a_c[N - i] / a_c[0]);
    const auto p_c_evals = consteig::eigenvalues(A_pole);

    // Step 4: map poles to z-domain; build monic denominator polynomial.
    Complex p_d_vals[N]{};
    Complex pole_poly[N + 1u]{};
    pole_poly[0] = Complex{static_cast<T>(1), static_cast<T>(0)};
    for (consteig::Size k = 0; k < N; ++k)
    {
        p_d_vals[k] =
            consteig::exp(p_c_evals(k, 0) * Complex{Ts, static_cast<T>(0)});
        const Complex zk = p_d_vals[k];
        pole_poly[k + 1u] =
            Complex{static_cast<T>(0), static_cast<T>(0)} - zk * pole_poly[k];
        for (consteig::Size i = k; i > 0u; --i)
        {
            pole_poly[i] = pole_poly[i] - zk * pole_poly[i - 1u];
        }
    }

    // Steps 5-7: find analog zeros (roots of b_c) via companion matrix,
    // same trick as above. Embedded in NxN so consteig::eigenvalues can be
    // called with a fixed size; the top-left nz x nz block holds the actual
    // companion, the remaining entries are 0, producing d_b spurious
    // eigenvalues at the origin that are discarded afterward.
    Complex z_c_finite[N]{};
    Complex z_d_finite[N]{};
    Complex zero_poly[N + 1u]{};
    zero_poly[0] = Complex{static_cast<T>(1), static_cast<T>(0)};

    if (nz > 0u)
    {
        consteig::Matrix<T, N, N> A_zero{};
        for (consteig::Size i = 0; i + 1u < nz; ++i)
            A_zero(i + 1u, i) = static_cast<T>(1);
        for (consteig::Size i = 0; i < nz; ++i)
            A_zero(i, nz - 1u) = -(b_c[d_b + (nz - i)] / b_c[d_b]);
        const auto z_c_evals = consteig::eigenvalues(A_zero);

        // Mark d_b spurious eigenvalues (smallest magnitude = from zero block).
        bool spurious[N]{};
        for (consteig::Size m = 0; m < d_b; ++m)
        {
            consteig::Size min_idx = 0;
            T min_mag_sq = static_cast<T>(-1);
            for (consteig::Size i = 0; i < N; ++i)
            {
                if (!spurious[i])
                {
                    const T mag_sq =
                        z_c_evals(i, 0).real * z_c_evals(i, 0).real +
                        z_c_evals(i, 0).imag * z_c_evals(i, 0).imag;
                    if (min_mag_sq < static_cast<T>(0) || mag_sq < min_mag_sq)
                    {
                        min_mag_sq = mag_sq;
                        min_idx = i;
                    }
                }
            }
            spurious[min_idx] = true;
        }

        // Collect finite zeros and build zero polynomial prod(z - z_d_k).
        consteig::Size nz_cnt = 0u;
        for (consteig::Size k = 0; k < N; ++k)
        {
            if (!spurious[k])
            {
                z_c_finite[nz_cnt] = z_c_evals(k, 0);
                z_d_finite[nz_cnt] = consteig::exp(
                    z_c_evals(k, 0) * Complex{Ts, static_cast<T>(0)});
                const Complex zk = z_d_finite[nz_cnt];
                zero_poly[nz_cnt + 1u] =
                    Complex{static_cast<T>(0), static_cast<T>(0)} -
                    zk * zero_poly[nz_cnt];
                for (consteig::Size i = nz_cnt; i > 0u; --i)
                {
                    zero_poly[i] = zero_poly[i] - zk * zero_poly[i - 1u];
                }
                ++nz_cnt;
            }
        }
    }

    // Step 8: pad with zeros at z = -1 to reach numerator degree N-1.
    const consteig::Size n_extra = (nz + 1u < N) ? (N - nz - 1u) : 0u;
    for (consteig::Size e = 0; e < n_extra; ++e)
    {
        const consteig::Size cur_deg = nz + e;
        zero_poly[cur_deg + 1u] = zero_poly[cur_deg + 1u] + zero_poly[cur_deg];
        for (consteig::Size i = cur_deg; i > 0u; --i)
        {
            zero_poly[i] = zero_poly[i] + zero_poly[i - 1u];
        }
    }
    const consteig::Size num_deg = nz + n_extra;

    // Step 9: find matching frequency w_c (avoid collision with poles/zeros).
    const T tol = static_cast<T>(1.49e-8);
    T w_c = static_cast<T>(0);
    for (consteig::Size attempt = 0; attempt < 1000u; ++attempt)
    {
        bool collision = false;
        for (consteig::Size i = 0; i < N && !collision; ++i)
        {
            const T dr = w_c - p_c_evals(i, 0).real;
            const T di = p_c_evals(i, 0).imag;
            if (dr * dr + di * di < tol * tol)
            {
                collision = true;
            }
        }
        for (consteig::Size i = 0; i < nz && !collision; ++i)
        {
            const T dr = w_c - z_c_finite[i].real;
            const T di = z_c_finite[i].imag;
            if (dr * dr + di * di < tol * tol)
            {
                collision = true;
            }
        }
        if (!collision)
        {
            break;
        }
        w_c += static_cast<T>(0.1) / Ts;
    }

    // Step 10: compute discrete gain k_d matching H_d(w_d) = H_c(w_c).
    const Complex w_c_cx{w_c, static_cast<T>(0)};
    const Complex w_d_cx =
        consteig::exp(w_c_cx * Complex{Ts, static_cast<T>(0)});

    Complex num_c_cx{static_cast<T>(1), static_cast<T>(0)};
    for (consteig::Size i = 0; i < nz; ++i)
    {
        num_c_cx = num_c_cx * (w_c_cx - z_c_finite[i]);
    }

    Complex den_c_cx{static_cast<T>(1), static_cast<T>(0)};
    for (consteig::Size i = 0; i < N; ++i)
    {
        den_c_cx = den_c_cx * (w_c_cx - p_c_evals(i, 0));
    }

    Complex num_d_cx{static_cast<T>(1), static_cast<T>(0)};
    for (consteig::Size i = 0; i < N; ++i)
    {
        num_d_cx = num_d_cx * (w_d_cx - p_d_vals[i]);
    }

    Complex den_d_cx{static_cast<T>(1), static_cast<T>(0)};
    for (consteig::Size i = 0; i < nz; ++i)
    {
        den_d_cx = den_d_cx * (w_d_cx - z_d_finite[i]);
    }
    for (consteig::Size e = 0; e < n_extra; ++e)
    {
        den_d_cx = den_d_cx *
                   (w_d_cx - Complex{static_cast<T>(-1), static_cast<T>(0)});
    }

    const Complex gain_num =
        Complex{k_c, static_cast<T>(0)} * num_c_cx * num_d_cx;
    const Complex gain_den = den_c_cx * den_d_cx;
    const T gain_den_sq =
        gain_den.real * gain_den.real + gain_den.imag * gain_den.imag;
    const T k_d =
        (gain_num.real * gain_den.real + gain_num.imag * gain_den.imag) /
        gain_den_sq;

    // Step 11: assemble output TF.
    TransferFunction<T, N + 1u, N + 1u> tf{};
    for (consteig::Size i = 0; i <= N; ++i)
    {
        tf.a[i] = pole_poly[i].real;
    }
    const consteig::Size pad = N - num_deg;
    for (consteig::Size i = 0; i <= num_deg; ++i)
    {
        tf.b[pad + i] = k_d * zero_poly[i].real;
    }

    return tf;
}

// Backward-compatible wrapper: recovers (b_c, a_c) from SS and delegates.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> matched_z_discretize(
    const StateSpace<T, N> &sys_c, T Ts, MatchedZ /*tag*/)
{
    T a_c[N + 1u]{};
    char_poly(sys_c.A, a_c);
    T b_c[N + 1u]{};
    markov_numerator(sys_c, a_c, b_c);
    return matched_z_discretize_tf<T, N>(b_c, a_c, Ts, MatchedZ{});
}

// --------------------------- analog_to_digital (state-space overloads) -------

template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const StateSpace<T, N> &sys_c, T Ts, ZOH)
{
    return ss_to_tf(zoh_discretize(sys_c, Ts, ZOH{}));
}

template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const StateSpace<T, N> &sys_c, T Ts, MatchedZ)
{
    return matched_z_discretize(sys_c, Ts, MatchedZ{});
}

// --------------------------- analog_to_digital (TF overloads) ----------------

template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T Ts, ZOH)
{
    return ss_to_tf(zoh_discretize(tf_to_ss<T, N>(b_c, a_c), Ts, ZOH{}));
}

template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T Ts, MatchedZ)
{
    return matched_z_discretize_tf<T, N>(b_c, a_c, Ts, MatchedZ{});
}

} // namespace constfilt

#endif // CONSTFILT_DISCRETIZE_HPP
