#ifndef CONSTFILT_STABILITY_HPP
#define CONSTFILT_STABILITY_HPP

#include "constfilt_options.hpp"
#include "discretize.hpp"
#include "vendor/consteig/consteig.hpp"

namespace constfilt
{

enum class Stability
{
    Stable,
    MarginallyStable,
    Unstable
};

// Checks continuous-time stability via eigenvalues of the A matrix.
//
//   Stable:           all poles have Re(λ) < 0
//   MarginallyStable: all poles have Re(λ) ≤ 0, and any imaginary-axis poles
//                     are simple (no repeated imaginary-axis poles)
//   Unstable:         any pole has Re(λ) > 0, or repeated imaginary-axis poles
//
// Tolerance used for "near imaginary axis" and "repeated" comparisons: 1e-8.
template <typename T, consteig::Size N>
constexpr Stability check_stability(const StateSpace<T, N> &sys)
{
    auto evals = consteig::eigenvalues(sys.A);

    constexpr T tol = static_cast<T>(1e-8);

    bool has_imag_axis = false;

    for (consteig::Size i = 0; i < N; ++i)
    {
        const T re = evals(i, 0).real;
        if (re > tol)
        {
            return Stability::Unstable;
        }
        if (re > -tol)
        {
            has_imag_axis = true;
        }
    }

    if (has_imag_axis)
    {
        // Check for repeated imaginary-axis poles (repeated → Unstable).
        for (consteig::Size i = 0; i < N; ++i)
        {
            if (evals(i, 0).real > -tol)
            {
                const T im_i = evals(i, 0).imag;
                for (consteig::Size j = i + 1u; j < N; ++j)
                {
                    if (evals(j, 0).real > -tol)
                    {
                        const T diff = im_i - evals(j, 0).imag;
                        const T absdiff =
                            diff < static_cast<T>(0) ? -diff : diff;
                        if (absdiff < tol)
                            return Stability::Unstable;
                    }
                }
            }
        }
        return Stability::MarginallyStable;
    }

    return Stability::Stable;
}

} // namespace constfilt

#endif // CONSTFILT_STABILITY_HPP
