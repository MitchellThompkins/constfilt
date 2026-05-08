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
// Coefficients are in descending power order:
//   coeff[0]*s^N + coeff[1]*s^{N-1} + ... + coeff[N]
//
// Template parameters:
//   T         - numeric type (float, double, ...)
//   N         - filter order (degree of denominator)
//   Method    - ZOH (default) or MatchedZ
//   CheckStab - when true (default), throws at construction (a compile-time
//               error for constexpr instances) if the analog system is
//               Unstable. Both Stable and MarginallyStable are accepted.
//               Set to false to skip.
//
// Constructors:
//   AnalogFilter(b_c, a_c, sample_rate_hz)
//     b_c            - s-domain numerator array [N+1], descending order
//     a_c            - s-domain denominator array [N+1], descending order
//     sample_rate_hz - sample rate in Hz
//
//   AnalogFilter(continuous_tf, sample_rate_hz)
//     continuous_tf  - s-domain TransferFunction (e.g. from a subclass)
//     sample_rate_hz - sample rate in Hz
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

    constexpr AnalogFilter(TransferFunction<T, N + 1u, N + 1u> continuous_tf,
                           T sample_rate_hz)
        : AnalogFilter(checked_discretize(continuous_tf.b, continuous_tf.a,
                                          sample_rate_hz))
    {
    }

  private:
    constexpr explicit AnalogFilter(
        TransferFunction<T, N + 1u, N + 1u> digital_tf)
        : Filter<T, N + 1u, N + 1u>(digital_tf.b, digital_tf.a)
    {
    }

    static constexpr TransferFunction<T, N + 1u, N + 1u> checked_discretize(
        const T (&b_c)[N + 1u], const T (&a_c)[N + 1u], T sample_rate_hz)
    {
        if (CheckStab &&
            check_stability(tf_to_ss<T, N>(b_c, a_c)) == Stability::Unstable)
            throw "constfilt: unstable analog filter";
        return analog_to_digital<T, N>(
            b_c, a_c, static_cast<T>(1) / sample_rate_hz, Method{});
    }
};

} // namespace constfilt

#endif // CONSTFILT_ANALOG_FILTER_HPP
