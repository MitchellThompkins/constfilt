# AGENTS.md — constfilt

Header-only C++17 library for compile-time IIR digital filter design.
All coefficient math is `constexpr`, built on [consteig](https://github.com/mitchellthompkins/consteig).

## Repo layout

```
include/constfilt/
  constfilt.hpp          — umbrella include (include this in user code)
  filter.hpp             — Filter<T,NB,NA>: DF2T, real-time + constexpr batch;
                           ConstexprArray<T,N> + batch_filter() for compile-time use
  butterworth.hpp        — Butterworth<T,N>: lowpass, orders 1–8
  analog_filter.hpp      — AnalogFilter<T,N,Method,CheckStab>: s-domain TF → digital
  discretize.hpp         — zoh_discretize, matched_z_discretize, ss_to_tf, StateSpace
  stability.hpp          — check_stability, Stability enum
  constfilt_options.hpp  — CONSTFILT_PI override

tests/
  filter.test.cpp
  discretize.test.cpp
  butterworth.test.cpp
  analog_filter.test.cpp
  butterworth_reference.hpp  — Octave-generated ground truth (do not edit by hand)

octave/
  generate_butterworth_tests.m  — regenerates butterworth_reference.hpp
```

## Building and testing

All development happens inside the consteig dev container (GCC 15 / Clang 21).
The container image is pulled automatically on first use.

**Run tests (preferred):**

```sh
make container.make.test.gcc
make container.make.test.clang
```

`container.make.<target>` runs `make <target>` inside the container.
Any Makefile target can be passed this way.

**Other useful targets via the container:**

```sh
make container.make.build.gcc     # build without running tests
make container.make.build.clang
make container.make.check-format  # clang-format dry-run (CI enforces this)
make container.make.format        # apply clang-format in-place
```

**Interactive shell inside the container:**

```sh
make container.start
# then inside the container:
make test.gcc
make test.clang
```

**Direct host build** (requires GCC/Clang and CMake installed locally):

```sh
make test.gcc
make test.clang
```

## Constraints

- **C++17 only.** No C++20 features (no `std::is_constant_evaluated`, no concepts).
- **Header-only.** Do not add `.cpp` translation units to the library itself.
- **`constexpr` everywhere possible.** Coefficient computation must work at compile
  time. The batch `Filter::operator()` must remain `constexpr`.
- **Template parameter type:** constfilt uses `unsigned int` for order `N`;
  consteig uses `consteig::Size` (= `size_t`). Implicit conversion is fine in C++17
  template arguments — do not add casts unless a compiler error requires it.
- **`a[0] = 1` convention:** all discrete transfer functions are monic-denominator.
  `Filter` stores `_a` with `a[0]` included for uniformity but the DF2T equations
  assume `a[0] = 1`.
- **Stability check:** `AnalogFilter` throws a string literal at runtime (or is a
  compile-time error) if the analog system is `Stability::Unstable`. Both `Stable`
  and `MarginallyStable` are accepted. Use `CheckStab = false` to skip.
- **No new dependencies.** consteig and googletest are the only external deps;
  both are managed by CMake `FetchContent`/`add_subdirectory`.
- **Format:** clang-format is enforced by CI. Run `make container.make.format`
  before committing.

## Compile-time verification in tests

Tests verify correctness at **both** compile time and runtime. The pattern:

```cpp
// 1. Coefficient checks — static_assert directly on constexpr filter
static constexpr constfilt::Butterworth<double, 2> filt(100.0, 1000.0);
static_assert(withinTol(filt.coeffs_b()[0], Ref::b[0], 1e-9), "b[0]");
EXPECT_NEAR(filt.coeffs_b()[0], Ref::b[0], CONSTFILT_COEFF_TOL);

// 2. Batch filtering — use ConstexprArray + batch_filter to capture output
static constexpr auto STEP32 = make_step<double, 32>();
static constexpr auto out = constfilt::batch_filter(filt, STEP32.data);
static_assert(all_within_tol(out, Ref::step, 1e-7), "step response");
// runtime check follows...
```

`ConstexprArray<T,N>` and `batch_filter()` live in `filter.hpp` (library API, not
test-only) — users can apply the same pattern in their own code.
`withinTol`, `all_within_tol`, and `make_step` are in `test_dependencies/test_tools.hpp`.

Functions like `char_poly` that use output parameters can be wrapped in a
constexpr immediately-invoked lambda to capture results for `static_assert`:

```cpp
constexpr auto poly = [&A]() {
    constfilt::ConstexprArray<double, 3> r{};
    constfilt::char_poly(A, r.data);
    return r;
}();
static_assert(withinTol(poly.data[0], 1.0, 1e-12), "coeffs[0]");
```

## Adding a new test

1. Create `tests/<name>.test.cpp`.
2. Add `add_constfilt_test(<name> <name>.test.cpp)` to `tests/CMakeLists.txt`.
3. Verify with `make container.make.test.gcc` and `make container.make.test.clang`.

## Regenerating Butterworth reference values

```sh
# inside the container or with Octave installed:
octave octave/generate_butterworth_tests.m
# overwrites tests/butterworth_reference.hpp
```

## Key algorithms (brief)

| Component | Algorithm |
|---|---|
| `zoh_discretize` | Matrix exponential via eigendecomposition; Bd via LU solve |
| `Butterworth` | Hardcoded normalized poles → companion-form SS → ZOH |
| `AnalogFilter` | s-domain TF → controllable-canonical SS → discretize → TF |
| `Filter` | Direct Form II Transposed (DF2T) |
| `batch_filter` | Thin wrapper around `Filter::operator()(input, output)` returning `ConstexprArray`; enables constexpr batch evaluation |
| Characteristic poly | Faddeev-LeVerrier (compute `c_k` from `p_{k-1}` *before* updating `M`) |
