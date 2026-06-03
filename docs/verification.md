# Verification

constfilt is verified by a layered set of GoogleTest suites that together
cover the library from individual math primitives up through end-to-end
filter output. The goal is to catch numerical errors and regressions at the
level where they originate, without relying on any single reference being
infallible.

The layers, from the inside out:

| Test file | Scope |
|---|---|
| `filter.test.cpp` | The `Filter` class in isolation: stored coefficients, batch output, sample-by-sample output, state reset. |
| `discretize.test.cpp` | `expm`, `char_poly`, `ss_to_tf`, ZOH, Matched-Z, each checked against closed-form analytic results. |
| `analog_filter.test.cpp` | Continuous-time transfer function to state-space conversion and stability classification; full filter output checks (coefficients, step, impulse, chirp) against Octave references. |
| `butterworth.test.cpp` | Full Butterworth pipeline: pole computation, state-space, discretization, filtering, lowpass and highpass, ZOH and Matched-Z. |
| `elliptic.test.cpp` | Full Elliptic pipeline: nome-series design, discretization, filtering, lowpass and highpass, ZOH and Matched-Z. |

## Philosophy

Octave is the primary reference. Scripts generate coefficients and filter
outputs and write them to committed C++ headers, so CI needs no Octave
installation and any disagreement shows up as a hard failure. See
[building.md](building.md#regenerating-reference-data) for how references are
regenerated.

Filter correctness is checked against step, impulse, and chirp inputs.

Because constfilt coefficients are `constexpr`, filter behavior must be
verified at compile time as well as runtime. The batch `operator()` is
`constexpr` and is exercised as such in tests; the sample-by-sample path
runs at runtime and is tested separately. Both paths are also checked for
bit-identical output to confirm they implement the same computation.

## Tolerances

| Quantity | Tolerance | Rationale |
|---|---|---|
| Transfer function coefficients vs. Octave | $10^{-7}$ | High-order coefficient computation accumulates floating-point error; N=8 was observed near $3 \times 10^{-8}$. |
| Step / impulse / chirp response vs. Octave | $10^{-7}$ | Recursive filtering accumulates rounding error; this bound accommodates that without masking real discrepancies. |
| Analytic ZOH / Matched-Z checks | $10^{-10}$ | Closed-form references carry no approximation error beyond `double` representation. |
| DF2T hand-trace | $10^{-13}$ to $10^{-14}$ | Rational inputs yield nearly exact results; tight bounds catch any arithmetic mistake. |
