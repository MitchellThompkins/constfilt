#ifndef CONSTFILT_ANALOG_FILTER_HPP
#define CONSTFILT_ANALOG_FILTER_HPP

#include "discretize.hpp"
#include "filter.hpp"
#include "stability.hpp"

namespace constfilt
{

// Discretize an analog (continuous-time, s-domain) transfer function into a
// digital Filter.
//
// b_c and a_c are Laplace-domain polynomial coefficients in descending power
// order:
//   coeff[0]*s^N + coeff[1]*s^{N-1} + ... + coeff[N]
//
// Template parameters:
//   T         - numeric type (float, double, …)
//   N         - filter order (degree of denominator)
//   Method    - ZOH (default) or MatchedZ
//   CheckStab - when true (default), throws at construction (a compile-time
//               error for constexpr instances) if the analog system is
//               Unstable. Both Stable and MarginallyStable are accepted.
//               Set to false to skip.
//
// Constructor:
//   AnalogFilter(b_c, a_c, sample_rate_hz)
//     b_c            - s-domain numerator coefficients [N+1], descending order
//     a_c            - s-domain denominator coefficients [N+1], descending
//     order sample_rate_hz - sample rate in Hz; Ts = 1/sample_rate_hz
template <typename T, consteig::Size N, typename Method = ZOH,
          bool CheckStab = true>
class AnalogFilter : public Filter<T, N + 1u, N + 1u>
{
    static_assert(N >= 1u, "Filter order must be at least 1");

  public:
    constexpr AnalogFilter(const T (&b_c)[N + 1u], const T (&a_c)[N + 1u],
                           T sample_rate_hz)
        : AnalogFilter(checked_discretize(b_c, a_c, sample_rate_hz))
    {
    }

  protected:
    // For subclasses that supply a pre-computed digital TF (e.g. Butterworth).
    constexpr explicit AnalogFilter(TransferFunction<T, N + 1u, N + 1u> tf)
        : Filter<T, N + 1u, N + 1u>(tf.b, tf.a)
    {
    }

    // Continuous-to-digital discretization; exposed so subclasses can reuse.
    static constexpr TransferFunction<T, N + 1u, N + 1u> discretize(
        const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T sample_rate_hz)
    {
        return analog_to_digital(tf_to_ss<T, N>(b_c, a_c),
                                 static_cast<T>(1) / sample_rate_hz, Method{});
    }

  private:
    static constexpr TransferFunction<T, N + 1u, N + 1u> checked_discretize(
        const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T sample_rate_hz)
    {
        if (CheckStab &&
            check_stability(tf_to_ss<T, N>(b_c, a_c)) == Stability::Unstable)
            throw "constfilt: unstable analog filter";
        return discretize(b_c, a_c, sample_rate_hz);
    }
};

} // namespace constfilt

#endif // CONSTFILT_ANALOG_FILTER_HPP
