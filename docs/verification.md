# Verification

constfilt uses three GoogleTest suites that together cover the library from
individual math primitives up through end-to-end filter output. The goal is to
catch numerical errors and regressions at the level where they originate -
without relying on any single reference being infallible.

---

## Philosophy

**Test against independent references.** The primary source of truth for
Butterworth filter coefficients and step-response values is Octave's `c2d`
function. An Octave script (`octave/generate_butterworth_tests.m`) generates
these values and writes them to a committed C++ header
(`tests/butterworth_reference.hpp`). constfilt's output must agree with Octave
within the declared tolerances; any discrepancy indicates a bug in the
implementation.

**Test at each level of the pipeline.** `butterworth.test.cpp` confirms the
full end-to-end result, but errors can be hidden when everything is tested only
at the top. `discretize.test.cpp` tests `expm`, `char_poly`, ZOH, and
Matched-Z against analytic results that require no external tool. `filter.test.cpp`
tests the DF2T state machine against a hand-traced reference. This layering
means a bug in `expm` shows up in `discretize.test`, not buried in a
coefficient mismatch.

**Use analytic references where possible.** For simple first-order systems,
the correct ZOH and Matched-Z results have closed-form expressions. These are
the sharpest possible tests because there is no numerical reference error to
absorb - the expected values are exact.

**Separate coefficient accuracy from filter accuracy.** Transfer-function
coefficients and step-response samples carry different numerical sensitivities.
Two tolerances reflect this: `1e-9` for coefficients, `1e-7` for step-response
samples. The looser step-response tolerance accounts for accumulated
floating-point error over 32 samples of recursive filtering.

---

## Test Suites

### `filter.test` - Direct Form II Transposed state machine

Tests the `Filter<T, NB, NA>` implementation in isolation using a hand-chosen
second-order IIR with known rational coefficients:

```
b = [0.25, 0.5, 0.25]
a = [1.0, -0.5, 0.0625]
```

| Test | What it checks |
|---|---|
| `FilterCoeffs::StoresB` | `coeffs_b()` returns the stored numerator coefficients |
| `FilterCoeffs::StoresA` | `coeffs_a()` returns the stored denominator coefficients |
| `FilterBatch::StepResponse4` | Batch output on a 4-sample unit step matches a hand-traced DF2T calculation (tolerance 1e-14) |
| `FilterRealTime::MatchesBatch` | Sample-by-sample output over 8 samples matches the batch path (tolerance 1e-13) |
| `FilterReset::StateIsZeroedAfterReset` | After `reset()`, the first output equals the first output from a freshly constructed filter |

The step-response reference values in `FilterBatch::StepResponse4` are derived
by tracing the DF2T recurrences by hand. Because the inputs and coefficients
are all rational numbers with exact `double` representations, agreement to
1e-14 is achievable.

---

### `discretize.test` - Math primitives and discretization methods

Tests `expm`, `char_poly`, `zoh_discretize`, `ss_to_tf`, and
`matched_z_discretize` against analytic references.

**ZOH - first-order systems**

For $H(s) = 1/(s+a)$ the analytic ZOH result is:

$$A_d = e^{-aT_s}, \quad B_d = \frac{1 - e^{-aT_s}}{a}, \quad C_d = 1, \quad D_d = 0$$

Two cases are tested (tolerance 1e-10):

| Test | $a$ | $T_s$ |
|---|---|---|
| `ZOH::FirstOrder_a1_T0p1` | 1 | 0.1 |
| `ZOH::FirstOrder_a5_T0p01` | 5 | 0.01 |

**Matrix exponential**

`Expm::Scalar` checks $e^{[-2]} = [e^{-2}] \approx 0.13533528$ (tolerance 1e-10).
A 1x1 matrix is the simplest non-trivial case and exercises the full
eigendecomposition path.

**Characteristic polynomial**

`CharPoly::DiagonalMatrix` checks `char_poly(diag(2,3))` returns
`[1, -5, 6]`, which is $(\lambda-2)(\lambda-3)$ expanded (tolerance 1e-12).

**State-space to transfer function**

`SsToTf::FirstOrder` checks a first-order discrete system
`A=[0.9], B=[0.1], C=[1], D=0` produces `a=[1, -0.9]`, `b=[0, 0.1]`
(tolerance 1e-12).

**Matched-Z - first-order system**

For $H(s) = 1/(s+1)$ with $T_s = 0.1$, the expected result is:

$$\text{pole: } z = e^{-0.1}, \quad \text{zero: } z = -1$$

$$b = \left[\frac{1 - e^{-0.1}}{2},\ \frac{1 - e^{-0.1}}{2}\right], \quad a = [1,\ -e^{-0.1}]$$

| Test | What it checks |
|---|---|
| `MatchedZ::FirstOrder_a1_T0p1` | Pole, zero, and coefficient values (tolerance 1e-10) |
| `MatchedZ::FirstOrder_DCGain` | $H_d(z=1) = (b[0]+b[1]) / (a[0]+a[1])$ equals 1.0 (tolerance 1e-10) |

The DC gain test is independent of the coefficient test and verifies the
gain-matching constraint directly.

---

### `butterworth.test` - End-to-end Butterworth filter

Tests the full pipeline - pole computation, continuous state-space
construction, ZOH discretization, and DF2T filtering - against Octave reference
values. Four cases span different filter orders and frequency ratios:

| Case | Order | Cutoff | Sample rate | $f_c/f_s$ |
|---|---|---|---|---|
| 1 | 2 | 100 Hz | 1000 Hz | 0.10 |
| 2 | 4 | 100 Hz | 1000 Hz | 0.10 |
| 3 | 2 | 500 Hz | 8000 Hz | 0.0625 |
| 4 | 3 | 200 Hz | 4000 Hz | 0.05 |

For each case, up to three tests are run:

- **Coefficients** - each `b[i]` and `a[i]` compared against the Octave
  reference (tolerance `1e-9`).
- **Batch** - 32-sample unit step response computed via the batch
  `operator()` compared against Octave's `filter()` output (tolerance `1e-7`).
- **Real-time** - same 32-sample step processed sample-by-sample compared
  against the same Octave reference (tolerance `1e-7`).

Cases 3 and 4 omit the explicit batch test because the real-time test covers
both paths (the real-time and batch paths are independently verified to produce
identical output in `filter.test`).

---

## Reference Data Generation

`octave/generate_butterworth_tests.m` produces `tests/butterworth_reference.hpp`.
The script uses Octave's `buttap`, `zp2tf`, and `c2d(..., 'zoh')` to design
and discretize each filter, then writes `constexpr` C++ structs with the
coefficient arrays and step-response values. The generated file is committed so
tests can run without an Octave installation.

To regenerate after changing test cases, run the script in Octave and commit
the updated header.

---

## Tolerances

| Quantity | Tolerance | Rationale |
|---|---|---|
| Transfer function coefficients vs. Octave | `1e-9` | Both implementations use `double`; the constexpr path and Octave's C-based routines agree to near the limit of `double` precision |
| Step response (32 samples) vs. Octave | `1e-7` | Recursive filtering accumulates rounding error over 32 steps; the looser bound accommodates this without masking real discrepancies |
| Analytic ZOH / Matched-Z checks | `1e-10` | Closed-form references carry no approximation error beyond `double` representation |
| DF2T hand-trace | `1e-13` to `1e-14` | Rational inputs yield nearly exact results; tight bounds catch any arithmetic mistake |
