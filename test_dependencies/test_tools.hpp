#ifndef CONSTFILT_TEST_TOOLS_HPP
#define CONSTFILT_TEST_TOOLS_HPP

// Absolute tolerance for comparing b/a transfer-function coefficients
// against the Octave ZOH reference.
#ifndef CONSTFILT_COEFF_TOL
#define CONSTFILT_COEFF_TOL 1e-9
#endif

// Absolute tolerance for comparing step-response samples against
// Octave's filter() output.
#ifndef CONSTFILT_STEP_TOL
#define CONSTFILT_STEP_TOL 1e-7
#endif

#include <cmath>
#include <constfilt/filter.hpp>
#include <cstddef>

// Absolute-value comparison helper (no std::abs in constexpr context).
template <typename T> static constexpr bool withinTol(T a, T b, T tol)
{
    T diff = a - b;
    return (diff < static_cast<T>(0) ? -diff : diff) < tol;
}

// Checks every element of a ConstexprArray against a matching raw-array
// reference. Returns true only if all elements are within tol.
template <typename T, size_t N>
constexpr bool all_within_tol(const constfilt::ConstexprArray<T, N> &computed,
                              const T (&ref)[N], T tol)
{
    for (size_t i = 0; i < N; ++i)
        if (!withinTol(computed.data[i], ref[i], tol))
            return false;
    return true;
}

// Returns a constexpr unit-step array of length N (all elements == 1).
template <typename T, size_t N>
constexpr constfilt::ConstexprArray<T, N> make_step()
{
    constfilt::ConstexprArray<T, N> s{};
    for (size_t i = 0; i < N; ++i)
        s.data[i] = static_cast<T>(1);
    return s;
}

#endif // CONSTFILT_TEST_TOOLS_HPP
