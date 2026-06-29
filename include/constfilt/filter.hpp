#ifndef CONSTFILT_FILTER_HPP
#define CONSTFILT_FILTER_HPP

#include "vendor/consteig/consteig.hpp"

namespace constfilt
{

// NB = number of numerator coefficients, NA = number of denominator
// coefficients. Convention: a[0] = 1 (monic denominator); stored for
// uniformity. Direct Form II Transposed implementation.
template <typename T, consteig::Size NB, consteig::Size NA> class Filter
{
    static constexpr consteig::Size M = (NA > NB ? NA : NB) - 1u; // STATE_LEN
    // NB == NA == 1 gives H(z) = b[0]/a[0] which is a pure gain. It is
    // mathematically valid but breaks the DF2T implementation (zero-sized
    // state, unsigned underflow in loops).
    static_assert(NB >= 1u,
                  "Filter requires at least one numerator coefficient");
    static_assert(NA >= 1u,
                  "Filter requires at least one denominator coefficient");
    static_assert(M >= 1u, "Filter requires max(NB, NA) - 1 >= 1");

  protected:
    T _b[NB]{};
    T _a[NA]{};
    mutable T _state[M]{}; // DF2T delay line, zero-initialized

  public:
    constexpr Filter() = default;

    constexpr Filter(const T (&b)[NB], const T (&a)[NA])
    {
        for (consteig::Size i = 0; i < NB; ++i)
        {
            _b[i] = b[i];
        }
        for (consteig::Size i = 0; i < NA; ++i)
        {
            _a[i] = a[i];
        }
    }

  private:
    constexpr T b_coeff(consteig::Size i) const
    {
        return i < NB ? _b[i] : T(0);
    }
    constexpr T a_coeff(consteig::Size i) const
    {
        return i < NA ? _a[i] : T(0);
    }

  public:
    // Real-time: process one sample. Declared const so the filter can be
    // static constexpr (coefficients at compile time) while still supporting
    // runtime sample-by-sample processing. _state is mutable for this reason.
    // Not thread-safe on shared instances without external synchronization.
    T operator()(T x) const
    {
        const T y = _b[0] * x + _state[0];

        for (consteig::Size k = 0; k < M - 1u; ++k)
        {
            _state[k] =
                b_coeff(k + 1u) * x - a_coeff(k + 1u) * y + _state[k + 1u];
        }
        _state[M - 1u] = b_coeff(M) * x - a_coeff(M) * y;

        return y;
    }

    // Batch: process an array with a local zero-initialized state.
    // constexpr-capable (does not touch member state).
    template <consteig::Size N>
    constexpr void operator()(const T (&input)[N], T (&output)[N]) const
    {
        T local_state[M]{};

        for (consteig::Size n = 0; n < N; ++n)
        {
            const T x = input[n];
            const T y = _b[0] * x + local_state[0];

            for (consteig::Size k = 0; k < M - 1u; ++k)
            {
                local_state[k] = b_coeff(k + 1u) * x - a_coeff(k + 1u) * y +
                                 local_state[k + 1u];
            }
            local_state[M - 1u] = b_coeff(M) * x - a_coeff(M) * y;

            output[n] = y;
        }
    }

    // Accessors (for testing)
    constexpr const T *coeffs_b() const
    {
        return _b;
    }
    constexpr const T *coeffs_a() const
    {
        return _a;
    }

    // Reset real-time state to zero.
    void reset() const
    {
        for (consteig::Size k = 0; k < M; ++k)
        {
            _state[k] = static_cast<T>(0);
        }
    }
};

} // namespace constfilt

#endif // CONSTFILT_FILTER_HPP
