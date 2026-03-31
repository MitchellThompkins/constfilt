#ifndef CONSTFILT_ANALOG_FILTER_HPP
#define CONSTFILT_ANALOG_FILTER_HPP

#include "discretize.hpp"
#include "filter.hpp"
#include "stability.hpp"

namespace constfilt
{

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

// Discretize an analog (continuous-time, s-domain) transfer function into a
// digital Filter.
//
// B and A are references to static-storage-duration constexpr arrays of
// Laplace-domain polynomial coefficients in descending power order:
//   coeff[0]*s^N + coeff[1]*s^{N-1} + ... + coeff[N]
//
// Because B and A are non-type template parameters, stability can be
// static_asserted at instantiation time with no runtime overhead.
//
// Template parameters:
//   T          - numeric type (float, double, …)
//   N          - filter order (degree of denominator)
//   B          - reference to static constexpr s-domain numerator array [N+1]
//   A          - reference to static constexpr s-domain denominator array [N+1]
//   Method     - ZOH (default) or MatchedZ
//   CheckStab  - when true (default), static_asserts that the analog system
//                is not Unstable at instantiation time. Both Stable and
//                MarginallyStable are accepted. Set to false to skip.
//
// Constructor:
//   AnalogFilter(sample_rate_hz)
//     sample_rate_hz - sample rate in Hz; Ts = 1/sample_rate_hz
template <typename T, consteig::Size N, const T (&B)[N + 1u],
          const T (&A)[N + 1u], typename Method = ZOH, bool CheckStab = true>
class AnalogFilter : public Filter<T, N + 1u, N + 1u>
{
    static_assert(N >= 1u, "Filter order must be at least 1");
    static_assert(!CheckStab || check_stability(tf_to_ss<T, N>(B, A)) !=
                                    Stability::Unstable,
                  "constfilt: unstable analog filter");

  public:
    constexpr AnalogFilter(T sample_rate_hz)
        : AnalogFilter(compute_ba(static_cast<T>(1) / sample_rate_hz))
    {
    }

  private:
    constexpr explicit AnalogFilter(TransferFunction<T, N + 1u, N + 1u> tf)
        : Filter<T, N + 1u, N + 1u>(tf.b, tf.a)
    {
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u> compute_ba(T Ts)
    {
        return discretize(tf_to_ss<T, N>(B, A), Ts, Method{});
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u> discretize(
        const StateSpace<T, N> &sys_c, T Ts, ZOH)
    {
        return ss_to_tf(zoh_discretize(sys_c, Ts, ZOH{}));
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u> discretize(
        const StateSpace<T, N> &sys_c, T Ts, MatchedZ)
    {
        return matched_z_discretize(sys_c, Ts, MatchedZ{});
    }
};

} // namespace constfilt

#endif // CONSTFILT_ANALOG_FILTER_HPP
