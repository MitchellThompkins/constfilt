# Getting Started

## Requirements

* A C++17 compiler (tested on GCC 15 and Clang 21 in the dev container).
* CMake 3.14 or newer if consuming via CMake.

constfilt is header-only and depends on
[consteig](https://github.com/MitchellThompkins/consteig) and
[gcem](https://github.com/MitchellThompkins/gcem). Both are fetched
automatically by the CMake build.

## Consuming the library

### add_subdirectory

```cmake
add_subdirectory(third_party/constfilt)
target_link_libraries(your_target PRIVATE constfilt::constfilt)
```

### FetchContent

```cmake
include(FetchContent)
FetchContent_Declare(
    constfilt
    GIT_REPOSITORY https://github.com/MitchellThompkins/constfilt
    GIT_TAG        main
)
FetchContent_MakeAvailable(constfilt)

target_link_libraries(your_target PRIVATE constfilt::constfilt)
```

`constfilt::constfilt` is an INTERFACE target that propagates the include path
and the `consteig` / `gcem` dependencies transitively.

### Header-only, no CMake

Add the constfilt, consteig, and gcem include directories to your compiler's
include path, then `#include <constfilt/constfilt.hpp>`.

## Hello, filter

```cpp
#include <constfilt/constfilt.hpp>

int main()
{
    // 4th-order Butterworth lowpass, 100 Hz cutoff at 1 kHz sample rate.
    // All coefficient math happens at compile time.
    static constexpr constfilt::Butterworth<double, 4> bw(100.0, 1000.0);

    // Compile-time batch filtering on a known input.
    static constexpr consteig::Array<double, 8> input{
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static constexpr auto output = bw(input);

    // Real-time sample-by-sample filtering (mutates internal state).
    constfilt::Butterworth<double, 4> rt(100.0, 1000.0);
    double y = rt(0.5);
}
```

## Supported filters

| Filter      | Variants            | Order  | Discretization      |
|-------------|---------------------|--------|---------------------|
| Butterworth | lowpass, highpass   | `N ≥ 1`| ZOH (default), MatchedZ |
| Elliptic    | lowpass, highpass   | `N ≥ 1`| ZOH (default), MatchedZ |

The library imposes no upper bound on the filter order. The practical upper
bound is set by `double` precision: as `N` grows, the magnitudes of polynomial
coefficients and pole separations stretch over many orders of magnitude and
round-off begins to dominate. The committed regression tests cover orders 1–8
because that is the range that is comfortably verifiable against Octave's
reference implementations.

### Selecting variants

```cpp
// Lowpass / highpass (template parameter 4):
constfilt::Butterworth<double, 4>                                       lp(100.0, 1000.0);
constfilt::Butterworth<double, 4, constfilt::ZOH, constfilt::HighPass>  hp(100.0, 1000.0);

// Discretization method (template parameter 3):
constfilt::Butterworth<double, 4, constfilt::MatchedZ>                  mz(100.0, 1000.0);

// Elliptic takes (cutoff_hz, passband_ripple_dB, stopband_atten_dB, sample_rate_hz):
constfilt::Elliptic<double, 4>                                          el(100.0, 0.5, 60.0, 1000.0);
constfilt::Elliptic<double, 4, constfilt::ZOH, constfilt::HighPass>     el_hp(100.0, 0.5, 60.0, 1000.0);

// Convenience aliases for first-order RC equivalents:
constfilt::FirstOrderLowPass<double>  first_lp(100.0, 1000.0);
constfilt::FirstOrderHighPass<double> first_hp(100.0, 1000.0);
```

### Accessing coefficients

The discretized transfer function is exposed as `coeffs_b()` and `coeffs_a()`
on the underlying `Filter`, both `constexpr`:

```cpp
static constexpr auto bw = constfilt::Butterworth<double, 2>(100.0, 1000.0);
static constexpr auto b  = bw.coeffs_b();  // numerator
static constexpr auto a  = bw.coeffs_a();  // denominator, a[0] == 1
```

## Configuration

| Macro | Default | Description |
|-------|---------|-------------|
| `CONSTFILT_PI` | `3.14159265358979323846` | Value of π used in frequency conversions. Define before including `<constfilt/constfilt.hpp>` to override. |

## Stability checking

`AnalogFilter` (the base of `Butterworth` and `Elliptic`) checks the
continuous-time poles before discretization. If the analog system is
`Stability::Unstable`, the constructor throws a string literal at runtime — or
fails to compile if the filter is `constexpr`. `Stable` and `MarginallyStable`
both pass. Pass `CheckStab = false` as an extra template parameter to skip
the check.

See [methods.md](methods.md) for the algorithms used to construct and
discretize the analog prototype, and [elliptic.md](elliptic.md) for the nome
series used by `Elliptic`.
