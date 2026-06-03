---
title: constfilt
---

# constfilt

![CI](https://github.com/mitchellthompkins/constfilt/actions/workflows/main.yml/badge.svg?branch=main)

constfilt is a header-only C++17 `constexpr` template library for designing IIR
digital filters at compile time. Coefficients are computed by the compiler from
the filter specification (cutoff, sample rate, ripple, etc...) and stored as
`static constexpr` values, so no processor time is spent computing them at
runtime and no offline tool (MATLAB, Python, SciPy, Octave, etc...) is needed to
generate them. The coefficients live directly in the program binary. If the
parameters change, the compiler recomputes them.

constfilt is header-only and depends only on
[consteig](https://github.com/MitchellThompkins/consteig) (compile-time linear
algebra) and [gcem](https://github.com/MitchellThompkins/gcem) (compile-time
math), both of which are vendored with constfilt.

## Quick Start

```cpp
#include <constfilt/constfilt.hpp>

// All coefficient math happens at compile time.
static constexpr constfilt::Butterworth<double, 4> bw(100.0, 1000.0);
static constexpr constfilt::Elliptic<double, 4>    el(100.0, 0.5, 60.0, 1000.0);

// Batch filtering:
double input[32]{ /* ... */ };
double output[32]{};
bw(input, output);

// Real-time sample-by-sample (mutates internal state):
constfilt::Butterworth<double, 4> rt(100.0, 1000.0);
double y = rt(x);
```

See [Getting Started](getting-started.md) for installation and a full walkthrough.

## What constfilt Supports

All at compile time:

- Butterworth lowpass and highpass filters of arbitrary order, from cutoff
  frequency and sample rate.
- Elliptic (Cauer) lowpass and highpass filters of arbitrary order, from
  cutoff frequency, passband ripple, and stopband attenuation.
- ZOH and Matched-Z discretization of arbitrary continuous-time
  transfer functions, via matrix exponential through eigendecomposition.
- Direct Form II Transposed filter implementation, with both a real-time
  sample-by-sample interface and a `constexpr`-capable batch interface.

Practical filter order is bounded only by `double` precision; the library
imposes no upper limit.

## When To Use constfilt

- Filter coefficients need to be known at compile time (embedded systems, safety-critical code).
- The standard library is limited or unavailable and a freestanding solution is desired.

## Why Does This Exist

Embedded systems often need fixed digital filters whose parameters are known at
compile time.  The usual workflow is: compute coefficients offline in MATLAB or
Python, hard-code them, and update them by hand whenever the specification
changes.

constfilt attempts to eliminate that workflow. The compiler does the math,
guarantees the coefficients match the specification, and recomputes them
automatically when the specification changes.

Computing discrete filter coefficients can be reframed as solving the roots of a
polynomial, which can be reframed as an eigenvalue problem, the same approach
used by LAPACK, and by extension MATLAB and Octave.
[consteig](https://github.com/MitchellThompkins/consteig) provides the
eigendecomposition at compile time.
