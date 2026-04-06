#ifndef CONSTFILT_FILTER_HPP
#define CONSTFILT_FILTER_HPP

#include <consteig/consteig.hpp>

namespace constfilt
{

// NB = number of numerator coefficients, NA = number of denominator
// coefficients. Convention: a[0] = 1 (monic denominator); stored for
// uniformity. Direct Form II Transposed implementation.
template <typename T, consteig::Size NB, consteig::Size NA> class Filter
{
    static constexpr consteig::Size M = (NA > NB ? NA : NB) - 1u; // STATE_LEN

  protected:
    T _b[NB]{};
    T _a[NA]{};
    T _state[M]{}; // DF2T delay line, zero-initialized

    constexpr Filter() = default;

    constexpr Filter(const T (&b)[NB], const T (&a)[NA])
    {
        for (consteig::Size i = 0; i < NB; ++i)
            _b[i] = b[i];
        for (consteig::Size i = 0; i < NA; ++i)
            _a[i] = a[i];
    }

  public:
    // Real-time: process one sample, mutate _state.
    // NOT constexpr (writes to member state).
    T operator()(T x)
    {
        T y = _b[0] * x + _state[0];

        for (consteig::Size k = 0; k < M - 1u; ++k)
        {
            _state[k] = _b[k + 1u] * x - _a[k + 1u] * y + _state[k + 1u];
        }
        _state[M - 1u] = _b[M] * x - _a[M] * y;

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
            T x = input[n];
            T y = _b[0] * x + local_state[0];

            for (consteig::Size k = 0; k < M - 1u; ++k)
            {
                local_state[k] =
                    _b[k + 1u] * x - _a[k + 1u] * y + local_state[k + 1u];
            }
            local_state[M - 1u] = _b[M] * x - _a[M] * y;

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
    void reset()
    {
        for (consteig::Size k = 0; k < M; ++k)
        {
            _state[k] = static_cast<T>(0);
        }
    }
};

// Value-semantic array wrapper enabling compile-time capture of batch filter
// output. Use with batch_filter() to evaluate Filter::operator() at compile
// time and verify results with static_assert.
template <typename T, consteig::Size N> struct ConstexprArray
{
    T data[N]{};
};

// Returns the output of filt(input) as a ConstexprArray. Because this
// function returns by value it is suitable for constexpr evaluation:
//
//   static constexpr auto out = constfilt::batch_filter(filt, input);
//   static_assert(out.data[0] == expected, "compile-time check");
template <typename FilterType, typename T, consteig::Size N>
constexpr ConstexprArray<T, N> batch_filter(const FilterType &filt,
                                            const T (&input)[N])
{
    ConstexprArray<T, N> out{};
    filt(input, out.data);
    return out;
}

} // namespace constfilt

#endif // CONSTFILT_FILTER_HPP
