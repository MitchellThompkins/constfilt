#ifndef CONSTFILT_DISCRETIZE_HPP
#define CONSTFILT_DISCRETIZE_HPP

#include "../../consteig/consteig.hpp"

namespace constfilt
{

// ─────────────────────────── Tag types ───────────────────────────────────────

struct ZOH
{
};

struct MatchedZ
{
};

// ─────────────────────────── Data structures ─────────────────────────────────

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
    consteig::Array<T, NB> b{};
    consteig::Array<T, NA> a{};
};

// ─────────────────────────── Matrix exponential ──────────────────────────────

// expm(A) via eigendecomposition:
//   Ad = V * diag(exp(λ_i)) * V^{-1}
// where V = eigenvectors, λ_i = eigenvalues (complex).
// Real part extracted at the end (imaginary parts cancel for real A).
template <typename T, consteig::Size N>
constexpr consteig::Matrix<T, N, N> expm(const consteig::Matrix<T, N, N> &A)
{
    using Cx = consteig::Complex<T>;
    using CxMat_NN = consteig::Matrix<Cx, N, N>;
    using CxMat_N1 = consteig::Matrix<Cx, N, 1>;

    // 1. Eigenvalues and eigenvectors
    auto evals = consteig::eigenvalues(A);       // Matrix<Cx, N, 1>
    auto V = consteig::eigenvectors(A, evals);   // Matrix<Cx, N, N>

    // 2. Invert V column-by-column via LU
    auto lu_V = consteig::lu(V);

    CxMat_NN V_inv{};
    for (consteig::Size j = 0; j < N; ++j)
    {
        CxMat_N1 e_j{};
        e_j(j, 0) = Cx{static_cast<T>(1), static_cast<T>(0)};
        auto col = consteig::lu_solve(lu_V, e_j);
        for (consteig::Size i = 0; i < N; ++i)
        {
            V_inv(i, j) = col(i, 0);
        }
    }

    // 3. Accumulate: expm = sum_i exp(λ_i) * v_i * w_i^T
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

// ─────────────────────────── ZOH discretization ──────────────────────────────

// ZOH: Ad = expm(Ac*Ts),  Bd = Ac^{-1} * (Ad - I) * Bc
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

    consteig::Matrix<T, N, N> Ad = expm(AcTs);

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

    // Bd = Ac^{-1} * rhs  →  solve Ac * Bd = rhs
    auto lu_Ac = consteig::lu(Ac);
    auto Bd = consteig::lu_solve(lu_Ac, rhs);

    StateSpace<T, N> sys_d{};
    sys_d.A = Ad;
    sys_d.B = Bd;
    sys_d.C = sys_c.C;
    sys_d.D = sys_c.D;
    return sys_d;
}

// ─────────────────────────── Faddeev-LeVerrier ───────────────────────────────

// Returns monic characteristic polynomial of Ad:
//   [1, c_1, c_2, ..., c_N]   (N+1 coefficients)
// Algorithm (Berkowitz / Faddeev-LeVerrier):
//   p_0 = I
//   c_k = -trace(Ad * p_{k-1}) / k
//   p_k = Ad * p_{k-1} + c_k * I
template <typename T, consteig::Size N>
constexpr consteig::Array<T, N + 1u> faddeev_leverrier(
    const consteig::Matrix<T, N, N> &Ad)
{
    consteig::Array<T, N + 1u> coeffs{};
    coeffs[0] = static_cast<T>(1);

    consteig::Matrix<T, N, N> M = consteig::eye<T, N>(); // p_0 = I

    for (consteig::Size k = 1u; k <= N; ++k)
    {
        // Ad * M  (= Ad * p_{k-1})
        consteig::Matrix<T, N, N> AdM{};
        for (consteig::Size r = 0; r < N; ++r)
        {
            for (consteig::Size c = 0; c < N; ++c)
            {
                T sum = static_cast<T>(0);
                for (consteig::Size j = 0; j < N; ++j)
                {
                    sum += Ad(r, j) * M(j, c);
                }
                AdM(r, c) = sum;
            }
        }

        // c_k = -trace(Ad * p_{k-1}) / k
        T tr = static_cast<T>(0);
        for (consteig::Size i = 0; i < N; ++i)
        {
            tr += AdM(i, i);
        }
        coeffs[k] = -tr / static_cast<T>(k);

        // p_k = Ad * p_{k-1} + c_k * I
        M = AdM;
        for (consteig::Size i = 0; i < N; ++i)
        {
            M(i, i) += coeffs[k];
        }
    }

    return coeffs;
}

// ─────────────────────────── Markov numerator ────────────────────────────────

// Computes numerator polynomial b from Markov parameters and denominator a.
//   h[0] = D
//   h[k] = C * A^{k-1} * B   for k = 1..N
//   b[k] = sum_{j=0}^{k} a[j] * h[k-j]
template <typename T, consteig::Size N>
constexpr consteig::Array<T, N + 1u> markov_numerator(
    const StateSpace<T, N> &sys_d, const consteig::Array<T, N + 1u> &a_coeffs)
{
    const auto &Ad = sys_d.A;
    const auto &Bd = sys_d.B;
    const auto &Cd = sys_d.C;

    // Compute Markov parameters h[0..N]
    consteig::Array<T, N + 1u> h{};
    h[0] = sys_d.D;

    // A_pow = A^{k-1}, starting with A^0 = I
    consteig::Matrix<T, N, N> A_pow = consteig::eye<T, N>();

    for (consteig::Size k = 1u; k <= N; ++k)
    {
        // h[k] = C * A_pow * B  (scalar)
        T val = static_cast<T>(0);
        for (consteig::Size c = 0; c < N; ++c)
        {
            // (A_pow * B)[c] = sum_j A_pow(c,j)*B(j,0)
            T AB_c = static_cast<T>(0);
            for (consteig::Size j = 0; j < N; ++j)
            {
                AB_c += A_pow(c, j) * Bd(j, 0);
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
                for (consteig::Size j = 0; j < N; ++j)
                {
                    sum += A_pow(r, j) * Ad(j, cc);
                }
                next(r, cc) = sum;
            }
        }
        A_pow = next;
    }

    // Convolve: b[k] = sum_{j=0}^{k} a[j] * h[k-j]
    consteig::Array<T, N + 1u> b{};
    for (consteig::Size k = 0; k <= N; ++k)
    {
        T sum = static_cast<T>(0);
        for (consteig::Size jj = 0; jj <= k; ++jj)
        {
            sum += a_coeffs[jj] * h[k - jj];
        }
        b[k] = sum;
    }

    return b;
}

// ─────────────────────────── ss_to_tf ────────────────────────────────────────

// Full pipeline: discrete state-space → (b, a) transfer function.
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> ss_to_tf(
    const StateSpace<T, N> &sys_d)
{
    auto a = faddeev_leverrier(sys_d.A);
    auto b = markov_numerator(sys_d, a);
    TransferFunction<T, N + 1u, N + 1u> tf{};
    tf.b = b;
    tf.a = a;
    return tf;
}

// ─── Matched-Z discretization ───────────────────────────────────────────────

// Matched-Z: discrete poles at z_k = exp(s_k * Ts).
// For all-pole continuous systems (no finite zeros), places N zeros at z = -1
// and matches DC gain: H_d(1) = H_c(0).
//
// Steps:
//   1. Ad = expm(Ac * Ts)               -- same pole mapping as ZOH
//   2. a  = faddeev_leverrier(Ad)        -- discrete denominator
//   3. H_c(0) = D - C * Ac^{-1} * B     -- continuous DC gain
//   4. b = K * (z+1)^N  where K = H_c(0) * a(1) / 2^N
template <typename T, consteig::Size N>
constexpr TransferFunction<T, N + 1u, N + 1u> matched_z_discretize(
    const StateSpace<T, N> &sys_c, T Ts, MatchedZ /*tag*/)
{
    // 1. Ad = expm(Ac * Ts)
    consteig::Matrix<T, N, N> AcTs{};
    for (consteig::Size r = 0; r < N; ++r)
        for (consteig::Size c = 0; c < N; ++c)
            AcTs(r, c) = sys_c.A(r, c) * Ts;
    const auto Ad = expm(AcTs);

    // 2. Discrete denominator polynomial
    const auto a = faddeev_leverrier(Ad);

    // 3. Continuous DC gain: solve Ac*x = B, then H_c(0) = D - C*x
    auto lu_Ac = consteig::lu(sys_c.A);
    auto x = consteig::lu_solve(lu_Ac, sys_c.B);
    T dc_gain = sys_c.D;
    for (consteig::Size j = 0; j < N; ++j)
        dc_gain -= sys_c.C(0, j) * x(j, 0);

    // 4. a(1) = sum of denominator coefficients
    T a_at_1 = static_cast<T>(0);
    for (consteig::Size k = 0; k <= N; ++k)
        a_at_1 += a[k];

    // 5. 2^N
    T two_pow_N = static_cast<T>(1);
    for (consteig::Size k = 0; k < N; ++k)
        two_pow_N *= static_cast<T>(2);

    // 6. Scale: K = dc_gain * a(1) / 2^N
    const T K = dc_gain * a_at_1 / two_pow_N;

    // 7. Numerator b[k] = K * C(N, k)  (binomial coefficients of (z+1)^N)
    consteig::Array<T, N + 1u> b{};
    T binom = static_cast<T>(1);
    for (consteig::Size k = 0; k <= N; ++k)
    {
        b[k] = K * binom;
        if (k < N)
            binom = binom * static_cast<T>(N - k) / static_cast<T>(k + 1u);
    }

    TransferFunction<T, N + 1u, N + 1u> tf{};
    tf.b = b;
    tf.a = a;
    return tf;
}

} // namespace constfilt

#endif // CONSTFILT_DISCRETIZE_HPP
