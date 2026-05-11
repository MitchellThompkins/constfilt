# AGENTS.md — constfilt

Header-only C++17 library for compile-time IIR digital filter design.
All coefficient math is `constexpr`, built on
[consteig](https://github.com/mitchellthompkins/consteig) and
[gcem](https://github.com/MitchellThompkins/gcem).

User-facing documentation lives in `README.md` and the `docs/` directory —
treat those as the source of truth for *what* the library does and how it is
used. This file is a guide to working on the code.

## Repo layout (by responsibility)

The library sources live under `include/constfilt/` as a flat collection of
headers. By responsibility:

- **Umbrella include** — pulls in everything; the entry point for user code.
- **Runtime filter** — Direct Form II Transposed implementation with both a
  sample-by-sample real-time path and a `constexpr`-capable batch path.
- **Analog front end** — converts a continuous-time transfer function to
  controllable-canonical state-space, then discretizes it. Used as a base
  class by the concrete filter designs.
- **Discretization** — `zoh_discretize`, `matched_z_discretize`, `ss_to_tf`,
  and the `StateSpace` struct.
- **Concrete designs** — Butterworth and Elliptic, both supporting lowpass
  and highpass, arbitrary order (`N >= 1`), and ZOH or Matched-Z
  discretization.
- **Stability classification** — `check_stability` and the `Stability` enum,
  used by the analog front end before discretization.
- **Options** — `CONSTFILT_PI` override.

Tests live under `tests/` (GoogleTest); numerical references are committed as
generated C++ headers next to the tests that consume them. The generators
live under `octave/`. See [docs/verification.md](docs/verification.md) for
the layered testing approach.

## Building and testing

All development happens inside the consteig dev container (GCC 15 / Clang 21).
The container image is pulled automatically on first use.

```sh
make container.make.test.gcc        # configure, build, run tests with GCC
make container.make.test.clang      # ... with Clang
make container.make.check-format    # clang-format dry-run (CI enforces this)
make container.make.format          # apply clang-format in-place
```

`container.make.<target>` runs `make <target>` inside the container, so any
host Makefile target can be invoked this way. For an interactive container
shell, use `make container.start`. A direct host build is supported (same
top-level targets without the `container.make.` prefix) but CI runs only in
the container.

See [docs/building.md](docs/building.md) for the full set of targets and the
reference-data regeneration workflow.

## Constraints

- **C++17 only.** No C++20 features (no `std::is_constant_evaluated`, no
  concepts).
- **Header-only.** Do not add `.cpp` translation units to the library itself.
- **`constexpr` everywhere possible.** Coefficient computation must work at
  compile time. The batch `Filter::operator()` must remain `constexpr`.
- **No filter-order ceiling.** Filter classes only `static_assert(N >= 1)`.
  The practical upper bound is `double` precision; do not introduce
  hard-coded upper limits.
- **Template parameter type:** constfilt uses `unsigned int` for order `N`;
  consteig uses `consteig::Size` (= `size_t`). Implicit conversion is fine in
  C++17 template arguments — do not add casts unless a compiler error
  requires it.
- **`a[0] = 1` convention:** all discrete transfer functions are
  monic-denominator. `Filter` stores `_a` with `a[0]` included for uniformity
  but the DF2T equations assume `a[0] = 1`.
- **Stability check:** the analog front end throws a string literal at
  runtime (or is a compile-time error) if the analog system is
  `Stability::Unstable`. Both `Stable` and `MarginallyStable` are accepted.
  The `CheckStab = false` template parameter skips the check.
- **No new dependencies.** consteig, gcem, and googletest are the only
  external deps; all are managed by CMake `FetchContent`/`add_subdirectory`.
- **Format:** clang-format is enforced by CI. Run `make container.make.format`
  before committing.

## Key algorithms (brief)

| Component | Algorithm |
|---|---|
| `zoh_discretize` | Matrix exponential via eigendecomposition; $B_d$ via LU solve |
| `matched_z_discretize` | Map continuous poles via $z = e^{sT_s}$; finite zeros mapped to $e^{sT_s}$, missing zeros placed at $z = -1$; gain set to match continuous DC gain |
| `Butterworth` | Analytic pole formula $\theta_k = \pi(2k+N-1)/(2N)$ → real polynomial coefficients via conjugate-pair multiply-out → controllable-canonical SS → discretize |
| `Elliptic` | Nome-series Cauer design (see [docs/elliptic.md](docs/elliptic.md)) → s-domain TF → controllable-canonical SS → discretize |
| `AnalogFilter` | s-domain TF → controllable-canonical SS → discretize → TF |
| `Filter` | Direct Form II Transposed (DF2T) |
| Characteristic poly | Eigenvalues from `consteig::eigenvalues` expanded as $\prod_k(z - \lambda_k)$ in complex arithmetic; real parts taken at the end |
