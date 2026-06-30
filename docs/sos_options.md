# Matched-Z SOS: options for extra zeros and gain matching

## Background

The Matched-Z transform maps each analog pole and zero independently:

$$p_s \to e^{p_s T_s}, \qquad z_s \to e^{z_s T_s}$$

For a strictly-proper system (more poles than finite zeros), the discrete-time
result is also strictly proper.  That is mathematically complete — no extra
zeros are needed.

Octave's `__c2d__` adds zeros at $z = -1$ (the Nyquist frequency) to pad the
numerator up to degree $N-1$.  This is a convention, not a requirement.  The
library copies it for the direct-form (`SOS=false`) path so that coefficient
outputs match Octave references.

Octave also chooses a single test frequency for the whole filter when matching
gain — the smallest non-negative $\omega$ that avoids every pole and zero.

## The SOS cascade and numerator factoring

The SOS cascade evaluates as a product of section transfer functions:

$$H(z) = \prod_{i} H_i(z)$$

so the numerator of the cascade is the product of section numerators, and the
denominator is the product of section denominators.  This means that for the
cascade to equal the full-order Matched-Z filter exactly, the section
numerators must multiply out to the full numerator — not each section must
independently have the same numerator structure as the full filter.

The cascade's numerator is the **product** of section numerators, not the sum.
So the question is not "does each section have the same zeros as the full
filter" but "do the section numerators multiply out to the full-filter
numerator."  Option A achieves this by distributing zeros across sections so
that their product matches.  Option B gives each section no zeros, so the
product has no zeros either — a genuinely different filter.

### Example: N=4 LP Butterworth

The full-order Matched-Z maps the four analog poles to four discrete poles
$z_1, z_1^*, z_2, z_2^*$ and, by Octave convention, adds $N-1 = 3$ zeros at
$z = -1$:

$$H_\text{full}(z) = \frac{k \cdot (z+1)^3}{(z-z_1)(z-z_1^*)(z-z_2)(z-z_2^*)}$$

An SOS cascade of two biquads has numerator $B_1(z) \cdot B_2(z)$ and
denominator $(z^2 + a_{1,1}z + a_{1,2})(z^2 + a_{2,1}z + a_{2,2})$.  The
denominator matches $H_\text{full}$ because the poles are mapped identically.
For the numerators to match:

$$B_1(z) \cdot B_2(z) = k \cdot (z+1)^3$$

One distribution that achieves this is $B_1(z) = k_1(z+1)^2$ and
$B_2(z) = k_2(z+1)$, provided $k_1 k_2 = k$.  This is Option A.

## Why SOS breaks with independent per-section conventions

If each biquad independently applies the Octave convention — $n_\text{extra} =
\text{order} - n_z - 1$ — then each 2-pole LP section (with 0 finite zeros)
gets one zero at $z=-1$:

$$H_i(z) = \frac{k_i(z+1)}{z^2 + a_{i,1}z + a_{i,2}}$$

Cascading two sections gives $(z+1)^2$ in the numerator, not $(z+1)^3$.  The
cascade has two zeros at $z=-1$ where the full-order filter has three.  These
are different filters.

The same problem applies to the gain test frequency: a section with a zero at
$s = 0$ shifts its test frequency away from $\omega = 0$ while other sections
in the same cascade do not, so the per-section gains do not compose to the
full-filter gain.

Both defects are pure implementation artifacts — the math is fine when the
conventions are applied globally.

## Option A: distribute extras globally (recommended)

Compute $n_\text{extra,total} = N - n_{z,\text{total}} - 1$ and $\omega_c$
from the full-filter `FactoredTF` before building any section.  Distribute
$n_\text{extra,total}$ zeros at $z = -1$ across sections (greedily: each
biquad takes up to 2, each real section up to 1).  Pass $\omega_c$ to every
section.

For N=4 LP this gives section 1 two zeros at $z=-1$ and section 2 one zero,
so the cascade numerator is $(z+1)^3$ — identical to the full-order filter.

**Result**: `SOS=true` and `SOS=false` produce identical output for all filter
types and all orders.

**Cost**: the SOS constructor must build a `FactoredTF<T, N>` for the full
filter (already computed by the Butterworth and Elliptic designs).  The section
builders need explicit $n_\text{extra}$ and $\omega_c$ parameters, bypassing
`AnalogFilter` for the Matched-Z path.  Small code addition; no new reference
data needed.

## Option B: pure Matched-Z per section (no padding, currently implemented)

Do not add extra zeros at $z = -1$ in SOS sections at all.  Each section's
transfer function is strictly proper.  A global test frequency is still used
for gain matching (fixing cross-section gain inconsistency), but
$n_\text{extra} = 0$ everywhere.

For N=4 LP each biquad numerator is a constant $k_i$, so the cascade numerator
is also a constant — a degree-0 polynomial.  The full-order Matched-Z
numerator is degree 3.  The two realizations have the same poles but different
zeros, so they are genuinely different filters that happen to share the same
analog prototype.

**Result**: `SOS=true` and `SOS=false` produce different output.  The SOS
filter has no notch at Nyquist; the direct form does.  For lowpass designs with
$f_c \ll f_s/2$ the difference is negligible in the passband but measurable.
Separate SOS-specific reference data is required.

**Cost**: new Octave-generated references (`case_sos_mz_*`) are needed; the
`ButterworthSOS_MZ_LP` test cases use those instead of the direct-form
references.

## Option C: static_assert MatchedZ + SOS

Treat Matched-Z as non-separable over the SOS decomposition and forbid the
combination at compile time.

**Result**: users must choose `SOS=false` for Matched-Z, or switch to Tustin.

**Cost**: breaks any downstream code using the default `SOS=true` with
MatchedZ.  Weakens the library — Matched-Z SOS is a valid and useful
combination, just requiring a little care to implement correctly.
