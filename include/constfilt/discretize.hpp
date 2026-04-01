#ifndef CONSTFILT_DISCRETIZE_HPP
#define CONSTFILT_DISCRETIZE_HPP

#include <consteig/consteig.hpp>

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
    using Cx = consteig::Complex<T>;
    using CxMat_NN = consteig::Matrix<Cx, N, N>;
    using CxMat_N1 = consteig::Matrix<Cx, N, 1>;

    // 1. Eigenvalues and eigenvectors
    auto evals = consteig::eigenvalues(A);     // Matrix<Cx, N, 1>
    auto V = consteig::eigenvectors(A, evals); // Matrix<Cx, N, N>

    // 2. Invert V column-by-column via LU
    auto lu_V = consteig::lu(V);

    CxMat_NN V_inv{};
    for (consteig::Size col = 0; col < N; ++col)
    {
        CxMat_N1 e_col{};
        e_col(col, 0) = Cx{static_cast<T>(1), static_cast<T>(0)};
        auto col_vec = consteig::lu_solve(lu_V, e_col);
        for (consteig::Size row = 0; row < N; ++row)
        {
            V_inv(row, col) = col_vec(row, 0);
        }
    }

    // 3. Accumulate: matrix_exp = sum_i exp(lam_i) * v_i * w_i^T
    //    where v_i = column i of V, w_i^T = row i of V_inv
    CxMat_NN result_c{};
    for (consteig::Size i = 0; i < N; ++i)
    {
        Cx exp_lambda = consteig::exp(evals(i, 0));
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
    using Cx = consteig::Complex<T>;

    auto evals = consteig::eigenvalues(Ad); // Matrix<Cx, N, 1>

    // p[0..k] holds the monic polynomial of degree k after k iterations.
    Cx p[N + 1u]{};
    p[0] = Cx{static_cast<T>(1), static_cast<T>(0)};

    for (consteig::Size k = 0; k < N; ++k)
    {
        const Cx lam = evals(k, 0);
        // Multiply degree-k poly by (z - lam), working high-to-low in-place.
        p[k + 1u] = Cx{static_cast<T>(0), static_cast<T>(0)} - lam * p[k];
        for (consteig::Size i = k; i > 0u; --i)
            p[i] = p[i] - lam * p[i - 1u];
        // p[0] is unchanged (stays 1)
    }

    for (consteig::Size i = 0; i <= N; ++i)
        coeffs[i] = p[i].real;
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
        e[k] = b[k] * inv_a0 - sys.D * (a[k] * inv_a0);

    // A: super-diagonal
    for (consteig::Size row = 0; row < N - 1u; ++row)
        sys.A(row, row + 1u) = static_cast<T>(1);

    // A: last row = -a[N-k]/a[0]  for k=0..N-1
    for (consteig::Size k = 0; k < N; ++k)
        sys.A(N - 1u, k) = -(a[N - k] * inv_a0);

    // B: last entry = 1
    sys.B(N - 1u, 0) = static_cast<T>(1);

    // C: C[0][k] = e[N-k]
    for (consteig::Size k = 0; k < N; ++k)
        sys.C(0, k) = e[N - k];

    return sys;
}

// --------------------------- analog_to_digital
// --------------------------------

// Discretize a continuous-time state-space via ZOH and extract the (b, a) TF.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const StateSpace<T, N> &sys_c, T Ts, ZOH)
{
    return ss_to_tf(zoh_discretize(sys_c, Ts, ZOH{}));
}

// Discretize a continuous-time state-space via matched-Z and extract the (b, a)
// TF.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> analog_to_digital(
    const StateSpace<T, N> &sys_c, T Ts, MatchedZ)
{
    return matched_z_discretize(sys_c, Ts, MatchedZ{});
}

// --- Matched-Z discretization -----------------------------------------------

// Matched-Z: discrete poles at z_k = exp(s_k * Ts).
// For all-pole continuous systems (no finite zeros), places N-1 zeros at
// z = -1 with b[0] = 0, and matches DC gain: H_d(1) = H_c(0).
//
// Steps:
//   1. Ad = matrix_exp(Ac * Ts)                  -- same pole mapping as ZOH
//   2. a  = char_poly(Ad)                   -- discrete denominator
//   3. H_c(0) = D - C * Ac^{-1} * B        -- continuous DC gain
//   4. b[0] = 0
//      b[k] = K * C(N-1, k-1)  for k = 1..N  (coefficients of (z+1)^{N-1})
//      where K = H_c(0) * a(1) / 2^{N-1}
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> matched_z_discretize(
    const StateSpace<T, N> &sys_c, T Ts, MatchedZ /*tag*/)
{
    // 1. Ad = matrix_exp(Ac * Ts)
    consteig::Matrix<T, N, N> AcTs{};
    for (consteig::Size r = 0; r < N; ++r)
        for (consteig::Size c = 0; c < N; ++c)
            AcTs(r, c) = sys_c.A(r, c) * Ts;
    const auto Ad = matrix_exp(AcTs);

    TransferFunction<T, N + 1u, N + 1u> tf{};

    // 2. Discrete denominator polynomial
    char_poly(Ad, tf.a);

    // 3. Continuous DC gain: solve Ac*x = B, then H_c(0) = D - C*x
    auto lu_Ac = consteig::lu(sys_c.A);
    auto x = consteig::lu_solve(lu_Ac, sys_c.B);
    T dc_gain = sys_c.D;
    for (consteig::Size col = 0; col < N; ++col)
        dc_gain -= sys_c.C(0, col) * x(col, 0);

    // 4. a(1) = sum of denominator coefficients
    T a_at_1 = static_cast<T>(0);
    for (consteig::Size k = 0; k <= N; ++k)
        a_at_1 += tf.a[k];

    // 5. 2^{N-1}
    T two_pow_Nm1 = static_cast<T>(1);
    for (consteig::Size k = 0; k < N - 1u; ++k)
        two_pow_Nm1 *= static_cast<T>(2);

    // 6. Scale: K = dc_gain * a(1) / 2^{N-1}
    const T K = dc_gain * a_at_1 / two_pow_Nm1;

    // 7. b[0] = 0; b[k] = K * C(N-1, k-1) for k = 1..N
    tf.b[0] = static_cast<T>(0);
    T binom = static_cast<T>(1);
    for (consteig::Size k = 1u; k <= N; ++k)
    {
        tf.b[k] = K * binom;
        if (k < N)
            binom = binom * static_cast<T>(N - k) / static_cast<T>(k);
    }

    return tf;
}

} // namespace constfilt

#endif // CONSTFILT_DISCRETIZE_HPP
