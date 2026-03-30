#ifndef CONSTFILT_FILTER_HPP
#define CONSTFILT_FILTER_HPP

#include <array>
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
    std::array<T, NB> _b{};
    std::array<T, NA> _a{};
    std::array<T, M> _state{}; // DF2T delay line, zero-initialized

    constexpr Filter() = default;

    constexpr Filter(const std::array<T, NB> &b,
                     const std::array<T, NA> &a)
        : _b(b), _a(a)
    {
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
    constexpr std::array<T, N> operator()(
        const std::array<T, N> &input) const
    {
        std::array<T, N> output{};
        std::array<T, M> local_state{};

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

        return output;
    }

    // Accessors (for testing)
    constexpr const std::array<T, NB> &coeffs_b() const
    {
        return _b;
    }
    constexpr const std::array<T, NA> &coeffs_a() const
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

} // namespace constfilt

#endif // CONSTFILT_FILTER_HPP
