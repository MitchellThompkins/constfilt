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

// Absolute-value comparison helper (no std::abs in constexpr context).
template <typename T>
static constexpr bool withinTol(T a, T b, T tol)
{
    T diff = a - b;
    return (diff < static_cast<T>(0) ? -diff : diff) < tol;
}

#endif // CONSTFILT_TEST_TOOLS_HPP
