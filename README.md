![ci](https://github.com/mitchellthompkins/constfilt/actions/workflows/main.yml/badge.svg)

# constfilt

constfilt is a header-only C++17 `constexpr` template library for designing
IIR digital filters at compile time. Coefficients are computed by the compiler
from the filter specification (cutoff, sample rate, ripple, ...) and stored as
`static constexpr` values, so no processor time is spent computing them at
runtime and no offline tool (MATLAB, Python, SciPy, Octave, ...) is needed to
generate them.

This is particularly useful when a filter's parameters are fixed at compile
time. The coefficients live directly in the program binary; if the parameters
change, the compiler recomputes them. When the math is wrong, `static_assert`
can catch it at build time. constfilt is header-only and depends only on
[consteig](https://github.com/MitchellThompkins/consteig) (compile-time linear
algebra) and [gcem](https://github.com/MitchellThompkins/gcem) (compile-time
math), both fetched automatically.

All at compile time, constfilt supports:

* **Butterworth** lowpass and highpass filters of arbitrary order, from cutoff
  frequency and sample rate.
* **Elliptic (Cauer)** lowpass and highpass filters of arbitrary order, from
  cutoff frequency, passband ripple, and stopband attenuation.
* **ZOH** and **Matched-Z** discretization of arbitrary continuous-time
  transfer functions, via matrix exponential through eigendecomposition.
* **Direct Form II Transposed** filtering, with both a real-time
  sample-by-sample path and a `constexpr`-capable batch path.

Practical filter order is bounded only by `double` precision; the library
itself imposes no upper limit.

Full documentation — usage, discretization methods, the elliptic algorithm,
and verification methodology — lives in [`docs/`](docs/).

# Usage

```cpp
#include <constfilt/constfilt.hpp>

// All coefficient math happens at compile time.
static constexpr constfilt::Butterworth<double, 4> bw(100.0, 1000.0);
static constexpr constfilt::Elliptic<double, 4>    el(100.0, 0.5, 60.0, 1000.0);

// Compile-time batch filtering:
static constexpr consteig::Array<double, 32> input = { /* ... */ };
static constexpr auto output = bw(input);

// Real-time sample-by-sample (mutates internal state):
constfilt::Butterworth<double, 4> rt(100.0, 1000.0);
double y = rt(x);
```

Add to a CMake project via `add_subdirectory`, `FetchContent`, or installed
package — see [`docs/getting-started.md`](docs/getting-started.md). Requires
C++17.

# Documentation

* [Getting started](docs/getting-started.md) — installation, supported filters,
  configuration options, full examples.
* [Building and testing](docs/building.md) — dev container, host build,
  regenerating reference data.
* [Discretization methods](docs/methods.md) — continuous → state-space → ZOH /
  Matched-Z → `b`, `a`.
* [Elliptic filter design](docs/elliptic.md) — the nome-based algorithm.
* [Verification](docs/verification.md) — test suites, references, tolerances.
