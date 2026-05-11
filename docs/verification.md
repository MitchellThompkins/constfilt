# Verification

constfilt is verified by a layered set of GoogleTest suites that together
cover the library from individual math primitives up through end-to-end
filter output. The goal is to catch numerical errors and regressions at the
level where they originate — without relying on any single reference being
infallible.

The layers, from the inside out:

| Layer | Scope |
|---|---|
| DF2T state machine | The runtime filter in isolation: stored coefficients, batch output, sample-by-sample output, state reset. |
| Math primitives & discretization | `expm`, `char_poly`, `ss_to_tf`, ZOH, Matched-Z — checked against closed-form analytic results. |
| Analog front end | Continuous-time transfer function → state-space conversion and stability classification. |
| End-to-end Butterworth | Pole computation → continuous state-space → discretization → DF2T filtering, lowpass and highpass, ZOH and Matched-Z. |
| End-to-end Elliptic | The nome-series algorithm (see [elliptic.md](elliptic.md)) → discretization → DF2T filtering, lowpass and highpass, ZOH and Matched-Z. |

---

## Philosophy

**Test against independent references.** The primary source of truth for
filter coefficients and step-response values is Octave. Octave scripts
generate the values and write them to committed C++ headers that the test
suite consumes directly; this means CI does not need an Octave installation,
and any disagreement between constfilt and Octave shows up as a hard test
failure. See [building.md](building.md#regenerating-reference-data) for how
the references are regenerated.

**Test at each level of the pipeline.** End-to-end suites confirm the full
result, but errors can be hidden when everything is tested only at the top.
The discretization-primitive suite tests `expm`, `char_poly`, ZOH, and
Matched-Z against analytic results that require no external tool. The DF2T
suite tests the filter state machine against a hand-traced reference. With
this layering, a bug in `expm` shows up as an `expm` failure, not as a
buried coefficient mismatch.

**Use analytic references where possible.** For simple first-order systems,
the correct ZOH and Matched-Z results have closed-form expressions. These are
the sharpest possible tests because there is no numerical reference error to
absorb — the expected values are exact.

**Separate coefficient accuracy from filter accuracy.** Transfer-function
coefficients and step-response samples carry different numerical sensitivities.
Two tolerances reflect this: `1e-9` for coefficients, `1e-7` for step-response
samples. The looser step-response tolerance accounts for accumulated
floating-point error over 32 samples of recursive filtering.

---

## What each layer checks

### DF2T state machine

A hand-chosen second-order IIR with rational coefficients
(`b = [0.25, 0.5, 0.25]`, `a = [1, -0.5, 0.0625]`) drives:

- Stored-coefficient round-trip via `coeffs_b()` / `coeffs_a()`.
- Batch output on a short unit step, compared against a DF2T recurrence
  traced by hand. Because all inputs and coefficients are exactly
  representable in `double`, agreement to ~1e-14 is achievable.
- Sample-by-sample output, compared against the batch path to ensure both
  evaluation paths produce identical numbers.
- `reset()` returns the filter to its freshly-constructed state.

### Math primitives & discretization

- **ZOH, first-order.** For $H(s) = 1/(s+a)$ the analytic ZOH result is
  $A_d = e^{-aT_s}$, $B_d = (1 - e^{-aT_s})/a$, $C_d = 1$, $D_d = 0$.
  Multiple $(a, T_s)$ pairs are checked.
- **Matrix exponential.** A 1×1 case ($e^{[-2]}$) is the simplest non-trivial
  test of the full eigendecomposition path.
- **Characteristic polynomial.** A diagonal matrix gives a closed-form
  characteristic polynomial against which `char_poly` is checked.
- **State-space → transfer function.** A first-order discrete system with
  hand-chosen $(A, B, C, D)$ is converted to its known $(b, a)$.
- **Matched-Z, first-order.** For $H(s) = 1/(s+1)$ at $T_s = 0.1$ the pole,
  zero, coefficients, and DC gain ($H_d(z=1) = 1$) are all checked
  independently.

### Analog front end

- Continuous-time transfer functions are converted to controllable-canonical
  state-space form. Strictly-proper and proper (non-zero `D`) cases are
  covered for both first- and second-order systems, including the case where
  the denominator is unnormalized ($a[0] \neq 1$).
- The stability classifier is exercised against stable, unstable,
  marginally-stable (simple imaginary-axis poles), and unstable
  (repeated imaginary-axis poles) systems.
- End-to-end analog → digital discretization is checked against Octave-
  generated references for several continuous-time transfer functions
  (some strictly proper, some proper) under both ZOH and Matched-Z.
- The `CheckStab = false` escape hatch is exercised on an unstable system.

### End-to-end Butterworth

The full pipeline — pole computation, continuous state-space construction,
discretization, and DF2T filtering — is checked against Octave references
for a matrix of $(N, f_c, f_s)$ cases spanning different orders and
$f_c/f_s$ ratios. For each case, three checks run:

- **Coefficients** — each $b[i]$ and $a[i]$ against the reference
  (tolerance `1e-9`).
- **Batch** — 32-sample unit-step response from the batch `operator()`
  against Octave's `filter()` output (tolerance `1e-7`).
- **Real-time** — the same 32-sample step processed sample-by-sample
  (tolerance `1e-7`).

Lowpass and highpass variants are covered, as are ZOH and Matched-Z
discretization.

### End-to-end Elliptic

The nome-series algorithm (see [elliptic.md](elliptic.md)) is exercised for
a matrix of cases spanning multiple orders, two passband-ripple / stopband-
attenuation specifications, and both lowpass/highpass and ZOH/Matched-Z
variants. Each case runs the same coefficient / batch / real-time triple as
Butterworth.

An order-3, $R_p=10\,$dB, $R_s=60\,$dB case is included specifically to
verify that non-power-of-two orders and aggressive stopband attenuation are
handled correctly.

---

## Tolerances

| Quantity | Tolerance | Rationale |
|---|---|---|
| Transfer function coefficients vs. Octave | `1e-9` | Both implementations use `double`; the constexpr path and Octave's C-based routines agree to near the limit of `double` precision |
| Step response (32 samples) vs. Octave | `1e-7` | Recursive filtering accumulates rounding error over 32 steps; the looser bound accommodates this without masking real discrepancies |
| Analytic ZOH / Matched-Z checks | `1e-10` | Closed-form references carry no approximation error beyond `double` representation |
| DF2T hand-trace | `1e-13` to `1e-14` | Rational inputs yield nearly exact results; tight bounds catch any arithmetic mistake |
