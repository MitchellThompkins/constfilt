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
        : AnalogFilter(analog_to_digital(tf_to_ss<T, N>(B, A),
                                         static_cast<T>(1) / sample_rate_hz,
                                         Method{}))
    {
    }

  private:
    constexpr explicit AnalogFilter(TransferFunction<T, N + 1u, N + 1u> tf)
        : Filter<T, N + 1u, N + 1u>(tf.b, tf.a)
    {
    }
};

} // namespace constfilt

#endif // CONSTFILT_ANALOG_FILTER_HPP
