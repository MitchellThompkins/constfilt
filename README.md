![ci](https://github.com/mitchellthompkins/constfilt/actions/workflows/main.yml/badge.svg)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://mitchellthompkins.github.io/constfilt/)

# constfilt

constfilt is a header-only C++17 `constexpr` template library for designing IIR
digital filters at compile time. Coefficients are computed by the compiler from
the filter specification (cutoff, sample rate, ripple, etc...) and stored as
`static constexpr` values, so no processor time is spent computing them at
runtime and no offline tool (MATLAB, Python, SciPy, Octave, etc...) is needed to
generate them. The coefficients live directly in the program binary. If the
parameters change, the compiler recomputes them.

constfilt depends only on two other header-only libraries both of which are
vendored with constfilt:
1. [consteig](https://github.com/MitchellThompkins/consteig) (compile-time
   linear algebra)
2. [gcem](https://github.com/MitchellThompkins/gcem) (compile-time math)

All at compile time, constfilt supports:

- Butterworth lowpass and highpass filters of arbitrary order, from cutoff
  frequency and sample rate.
- Elliptic (Cauer) lowpass and highpass filters of arbitrary order, from
  cutoff frequency, passband ripple, and stopband attenuation.
- ZOH and Matched-Z discretization of arbitrary continuous-time
  transfer functions, via matrix exponential through eigendecomposition.
- Direct Form II Transposed filter implementation, with both a real-time
  sample-by-sample interface and a `constexpr`-capable batch interface.

Full documentation can be found at
[documentation](https://mitchellthompkins.github.io/constfilt/).

# Usage

Simply `#include <constfilt/constfilt.hpp>` (optionally consume via CMake with
`add_subdirectory` or `FetchContent`). Requires C++17.

```cpp
#include <constfilt/constfilt.hpp>

static constexpr constfilt::Butterworth<double, 4> bw(100.0, 1000.0);
static constexpr constfilt::Elliptic<double, 4>    el(100.0, 0.5, 60.0, 1000.0);
```

For a quick reference on getting started, including examples, see
[getting-started](https://mitchellthompkins.github.io/constfilt/getting-started/).
