# Elliptic Filter Coefficient Computation

`Elliptic<T, N, Method, FilterType>` computes a continuous-time elliptic (Cauer)
transfer function in the s-domain, then discretizes it via the parent class
`AnalogFilter`. All arithmetic is `constexpr`.

The elliptic prototype is designed from passband ripple $R_p$ and stopband
attenuation $R_s$ specifications. Poles and zeros are placed using Jacobi
elliptic functions to achieve equiripple behavior in both the passband and
stopband, giving the steepest possible transition for a given order and ripple
budget.

The core parameter driving the design is the nome $q$, a small positive number
derived from the filter specifications that controls all of the theta-function
series used to place poles and zeros. The pole/zero placement series are taken
from Octave's `ncauer` algorithm, with one departure in how $q$ is obtained.
`ncauer` solves the degree equation iteratively via `fminbnd` to get the
stopband edge $w_s$, derives the design modulus $k = 1/w_s$, and computes $q$
from $k$. This implementation instead applies the modular identity
$q = q_1^{1/N}$, where $q_1$ is the nome of the selectivity modulus
$k_1 = \epsilon_p/\epsilon_s$. The identity is an exact algebraic consequence
of the nome definition and the degree equation, so no solver is needed, which
is the natural fit for a `constexpr` computation. From
that point on, the $\sigma_0$ series, the $w_i$ series, and the pole/zero
placement all follow `ncauer` implementation.

## Notation

| Symbol | Meaning |
|--------|---------|
| $N$ | filter order |
| $M$ | $\lfloor N/2 \rfloor$, number of complex-conjugate pole/zero pairs |
| $R_p$ | passband ripple in dB |
| $R_s$ | stopband attenuation in dB |
| $\omega_c$ | passband edge in rad/s ($2\pi f_c$) |
| $\epsilon_p$ | passband ripple factor: $\sqrt{10^{R_p/10} - 1}$ |
| $\epsilon_s$ | stopband ripple factor: $\sqrt{10^{R_s/10} - 1}$ |
| $k_1$ | selectivity modulus: $\epsilon_p / \epsilon_s$ |
| $k$ | design modulus recovered from $q$ via theta functions |
| $w_s$ | stopband edge of the normalized prototype: $1/k$ |
| $q$ | nome (elliptic parameter driving all series) |
| $\sigma_0$ | real part of the normalized prototype poles |

## Step 1: Nome via modular identity

The selectivity modulus is $k_1 = \epsilon_p / \epsilon_s$.

The nome of $k_1$ is computed by the q-series approximation used in `ncauer`:

$$k_1' = \sqrt{1 - k_1^2}$$

$$q_0 = \frac{1}{2} \cdot \frac{1 - \sqrt{k_1'}}{1 + \sqrt{k_1'}}$$

$$q_1 = q_0 + 2q_0^5 + 15q_0^9 + 150q_0^{13}$$

The design nome $q$ is then obtained by the modular identity:

$$q = q_1^{1/N} = \exp\!\left(\frac{\ln q_1}{N}\right)$$

The elliptic modular equation $K'(k)/K(k) = N \cdot K'(k_1)/K(k_1)$ is equivalent in
nome space to $q = q_1^{1/N}$. This is an exact identity, so no iteration or
convergence criterion is needed and no error accumulates across iteration steps.

## Step 2: Design modulus from nome

The design modulus $k$ is recovered from $q$ via Jacobi theta functions:

$$\vartheta_2(q) = 2q^{1/4} \sum_{n=0}^{\infty} q^{n(n+1)}$$

$$\vartheta_3(q) = 1 + 2\sum_{n=1}^{\infty} q^{n^2}$$

$$k = \left(\frac{\vartheta_2(q)}{\vartheta_3(q)}\right)^2$$

Both series are summed to $n = 30$, which gives full double precision for $q < 0.5$.
The stopband edge of the normalized prototype is $w_s = 1/k$.

## Step 3: Pole-shift parameter $\sigma_0$

$\sigma_0$ is the (normalized, real) common real part of all analog prototype poles.
It is computed from $R_p$ and $q$ by the theta-function series:

$$g = 10^{R_p/20}, \qquad \ell = \frac{1}{2N} \ln\!\left(\frac{g+1}{g-1}\right)$$

$$\sigma_{01} = \sum_{m=0}^{30} (-1)^m \, q^{m(m+1)} \sinh\!\bigl((2m+1)\ell\bigr)$$

$$\sigma_{02} = \sum_{m=1}^{30} (-1)^m \, q^{m^2} \cosh(2m\ell)$$

$$\sigma_0 = \left|\frac{2q^{1/4}\,\sigma_{01}}{1 + 2\sigma_{02}}\right|$$

Powers of $q$ are accumulated incrementally: $q^{m(m+1)}$ via ratio $q^{2m}$,
and $q^{m^2}$ via ratio $q^{2m-1}$.

## Step 4: Zero positions $w_i$

For each $i = 1, \ldots, M$, the normalized zero position $w_i$ on the real axis of the
unit circle (in Jacobi elliptic function space) is:

$$\mu = \begin{cases} i & \text{odd } N \\ i - \tfrac{1}{2} & \text{even } N \end{cases}$$

$$s_1 = 2q^{1/4} \sum_{m=0}^{30} (-1)^m \, q^{m(m+1)} \sin\!\left(\frac{(2m+1)\pi\mu}{N}\right)$$

$$s_2 = 2\sum_{m=1}^{30} (-1)^m \, q^{m^2} \cos\!\left(\frac{2m\pi\mu}{N}\right)$$

$$w_i = \frac{s_1}{1 + s_2}$$

## Step 5: Poles and zeros in the s-domain

For each pair index $i = 1, \ldots, M$:

$$V_i = \sqrt{(1 - kw_i^2)(1 - w_i^2/k)}, \qquad w = \sqrt{(1 + k\sigma_0^2)(1 + \sigma_0^2/k)}$$

$$\omega_z = \frac{\sqrt{w_s}}{w_i}, \qquad d = 1 + \sigma_0^2 w_i^2$$

$$p_\mathrm{re} = \frac{-\sigma_0 V_i \sqrt{w_s}}{d}, \qquad p_\mathrm{im} = \frac{w_i \, w \sqrt{w_s}}{d}$$

Zeros are at $\pm j\omega_z$ (purely imaginary, on the $j\omega$ axis).
Poles are at $p_\mathrm{re} \pm j p_\mathrm{im}$ (complex conjugate pairs, left half-plane).

For odd $N$, a single real pole is added at $-\sigma_0 \sqrt{w_s}$.

Polynomials $a(s)$ (denominator) and $b(s)$ (numerator) are built in
ascending-coefficient order by repeated application of `poly_mul_root`, which
multiplies the running polynomial by $(s - r)$ in-place. After all roots
are accumulated, coefficients are reversed to descending order and the imaginary
parts (which are zero by construction for real filter specs) are discarded.

## Step 6: Gain normalization

The unnormalized transfer function has arbitrary DC gain. It is normalized so that:

- Odd $N$: $H(0) = 1$
- Even $N$: $H(0) = G_p = 1/\sqrt{1 + \epsilon_p^2}$ (the passband edge gain)

This matches the convention used by Octave's `ellipap`. The normalization factor
applied to the numerator coefficients is:

$$\mathrm{gain} = H_0 \cdot \frac{a[N]}{b[N]}$$

where $a[N]$ and $b[N]$ are the constant terms ($s^0$ coefficients) of the
descending-order polynomials.

## Step 7: Cutoff frequency scaling

The normalized prototype has its passband edge at 1 rad/s. Scaling to $\omega_c$
rad/s is done by substituting $s \to s/\omega_c$, which multiplies the coefficient
of $s^i$ by $\omega_c^i$:

$$a[i] \mathrel{*}= \omega_c^i, \qquad b[i] \mathrel{*}= \omega_c^i$$

This is applied to both numerator and denominator.

## Highpass transform

The highpass overload computes the normalized lowpass prototype at $\omega_c = 1$
rad/s, then applies the LP-to-HP s-domain substitution $s \to \omega_c/s$. In
polynomial coefficient terms this is coefficient reversal combined with
$\omega_c^j$ scaling:

$$a_\mathrm{hp}[j] = a_\mathrm{lp}[N - j] \cdot \omega_c^j, \qquad b_\mathrm{hp}[j] = b_\mathrm{lp}[N - j] \cdot \omega_c^j$$

## Test reference and the `ncauer.m` shadow

`tests/elliptic_reference.hpp` is regenerated from Octave's `ellip()` via
`octave/generate_elliptic_tests.m`. The reference is independent of the
constfilt implementation: Octave computes elliptic coefficients via its own
numerical solver, and constfilt's modular-identity path is exercised separately.

One complication: `ncauer`'s internal `__ellip_ws` uses `fminbnd` with its
default `TolX`, which propagates ~1e-6 error to the analog pole locations, too
loose for the coefficient tolerance used in the assertions. `octave/ncauer.m` is therefore a
verbatim copy of upstream
[`ncauer.m`](https://github.com/gnu-octave/octave-signal/blob/main/inst/ncauer.m)
with a single-line patch that passes `optimset('TolX', 1e-15, 'MaxIter', 10000)`
to `fminbnd`. The generator script shadows the signal-package version with this
local copy; stock `ellip()` is called unchanged. With the tighter tolerance, both paths agree well within the assertion tolerance.

See [verification.md](verification.md) for the full tolerance rationale.
