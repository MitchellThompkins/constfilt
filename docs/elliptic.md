# Elliptic Filter Coefficient Computation

`Elliptic<T, N, Method, FilterType>` computes a continuous-time elliptic (Cauer)
transfer function in the s-domain, then discretizes it via the parent class
`AnalogFilter`. All arithmetic is `constexpr`.

The algorithm follows Octave's `ncauer` theta-function path. The classical
degree-equation approach solves `K'(k)/K(k) = N * K'(k1)/K(k1)` for the design
modulus `k` via iterative root finding on a transcendental scalar equation.
This implementation bypasses that entirely: the modular identity `q = q1^(1/N)`
is an analytic equivalence, not an approximation, so `q` is obtained in a single
closed-form expression at full machine precision with no iteration.

---

## Notation

| Symbol | Meaning |
|--------|---------|
| N      | filter order |
| M      | floor(N/2), number of complex-conjugate pole/zero pairs |
| Rp     | passband ripple in dB |
| Rs     | stopband attenuation in dB |
| wc     | passband edge in rad/s (2*pi*cutoff_hz) |
| ep     | passband ripple factor: sqrt(10^(Rp/10) - 1) |
| es     | stopband ripple factor: sqrt(10^(Rs/10) - 1) |
| k1     | selectivity modulus: ep / es |
| k      | design modulus recovered from q via theta functions |
| ws     | stopband edge of the normalized prototype: 1/k |
| q      | nome (elliptic parameter driving all series) |
| sig0   | real part of the normalized prototype poles |

---

## Step 1: Nome via modular identity

The selectivity modulus is `k1 = ep / es`.

The nome of `k1` is computed by the q-series approximation used in `ncauer`:

```
k1' = sqrt(1 - k1^2)
q0  = 0.5 * (1 - sqrt(k1')) / (1 + sqrt(k1'))
q1  = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13
```

The design nome `q` is then obtained by the modular identity:

```
q = q1^(1/N) = exp(log(q1) / N)
```

The elliptic modular equation `K'(k)/K(k) = N * K'(k1)/K(k1)` is equivalent in
nome space to `q = q1^(1/N)`. This is an exact identity, so no iteration or
convergence criterion is needed and no error accumulates across steps.

---

## Step 2: Design modulus from nome

The design modulus `k` is recovered from `q` via Jacobi theta functions:

```
theta2(q) = 2 * q^(1/4) * sum_{n=0}^{inf} q^(n*(n+1))
theta3(q) = 1 + 2 * sum_{n=1}^{inf} q^(n^2)
k = (theta2(q) / theta3(q))^2
```

Both series are summed to n=30, which gives full double precision for q < 0.5.
The stopband edge of the normalized prototype is `ws = 1/k`.

---

## Step 3: Pole-shift parameter sig0

`sig0` is the (normalized, real) common real part of all analog prototype poles.
It is computed from `Rp` and `q` by the theta-function series:

```
g  = 10^(Rp/20)       (power ratio, half the ripple dB)
l  = (1/(2N)) * log((g + 1) / (g - 1))

sig01 = sum_{m=0}^{30} (-1)^m * q^(m*(m+1)) * sinh((2m+1)*l)
sig02 = sum_{m=1}^{30} (-1)^m * q^(m^2)     * cosh(2*m*l)

sig0 = abs(2 * q^(1/4) * sig01 / (1 + 2*sig02))
```

Powers of q are accumulated incrementally: `q^(m*(m+1))` via ratio `q^(2m)`,
and `q^(m^2)` via ratio `q^(2m-1)`.

---

## Step 4: Zero positions wi

For each i = 1 ... M, the normalized zero position `wi` on the real axis of the
unit circle (in Jacobi elliptic function space) is:

```
mu    = i           (odd N)
mu    = i - 0.5     (even N)

soma1 = 2*q^(1/4) * sum_{m=0}^{30} (-1)^m * q^(m*(m+1)) * sin((2m+1)*pi*mu/N)
soma2 =             2 * sum_{m=1}^{30} (-1)^m * q^(m^2) * cos(2*m*pi*mu/N)

wi = soma1 / (1 + soma2)
```

---

## Step 5: Poles and zeros in the s-domain

For each pair index i = 1 ... M, the following intermediate quantities are
computed:

```
Vi      = sqrt((1 - k*wi^2) * (1 - wi^2/k))
omega_z = sqrt(ws) / wi
denom   = 1 + sig0^2 * wi^2
p_re    = sqrt(ws) * (-sig0 * Vi) / denom
p_im    = sqrt(ws) * (wi * w)     / denom
```

where `w = sqrt((1 + k*sig0^2) * (1 + sig0^2/k))`.

Zeros are at `+/- j*omega_z` (purely imaginary, on the jw axis).
Poles are at `p_re +/- j*p_im` (complex conjugate pairs, left half-plane).

For odd N, a single real pole is added at `-sig0 * sqrt(ws)`.

Polynomials `poly_a` (denominator) and `poly_b` (numerator) are built in
ascending-coefficient order by repeated application of `poly_mul_root`, which
multiplies the running polynomial by `(s - root)` in-place. After all roots
are accumulated, coefficients are reversed to descending order and the imaginary
parts (which are zero by construction for real filter specs) are discarded.

---

## Step 6: Gain normalization

The unnormalized transfer function has arbitrary DC gain. It is normalized so
that:

- Odd N: H(0) = 1
- Even N: H(0) = Gp = 1/sqrt(1 + ep^2)  (the passband edge gain)

This matches the convention used by Octave's `ellipap`. The normalization factor
applied to the numerator coefficients is:

```
gain = H0 * a[N] / b[N]
```

where `a[N]` and `b[N]` are the constant terms (s^0 coefficients) of the
descending-order polynomials.

---

## Step 7: Cutoff frequency scaling

The normalized prototype has its passband edge at 1 rad/s. Scaling to `wc`
rad/s is done by substituting `s -> s/wc`, which multiplies the coefficient of
`s^i` by `wc^i`:

```
a[i] *= wc^i
b[i] *= wc^i
```

This is applied to both numerator and denominator.

---

## Highpass transform

The highpass overload computes the normalized lowpass prototype at `wc = 1`
rad/s, then applies the LP-to-HP s-domain substitution `s -> wc/s`. In
polynomial coefficient terms this is coefficient reversal combined with
`wc^j` scaling:

```
a_hp[j] = a_lp[N - j] * wc^j
b_hp[j] = b_lp[N - j] * wc^j
```

---

## Discretization

After `elliptic_tf` returns the continuous-time coefficients, the `AnalogFilter`
base class constructs a state-space representation in controllable-canonical form
and discretizes it using the method selected by the `Method` template parameter
(ZOH by default). ZOH discretization uses matrix exponentiation via
eigendecomposition followed by an LU solve for the input matrix Bd. The discrete
transfer function is then extracted via Faddeev-LeVerrier and stored in the
`Filter` base.

---

## Test reference and the `ncauer.m` shadow

`tests/elliptic_reference.hpp` is regenerated from Octave's stock `ellip()` via
`octave/generate_elliptic_tests.m`. Since the reference needs to be independent
of the constfilt implementation (otherwise a bug common to both would not be
caught), the reference path uses Octave's own numerical elliptic-degree-equation
solver rather than constfilt's modular-identity path.

Inside the signal-package `ncauer.m`, the helper `__ellip_ws` solves the degree
equation via `fminbnd(@f, eps, 1-eps)` with no explicit options. `fminbnd`'s
default `TolX` (roughly 1e-4 on the search variable) propagates to the analog
pole locations at the ~1e-6 level, which is too loose to serve as a reference
for the `CONSTFILT_COEFF_TOL = 1e-9` assertions in `tests/elliptic.test.cpp`.

To be clear, this does **not** mean Octave's `ellip` is incorrect. It produces
a valid elliptic filter that meets the requested passband-ripple and
stopband-attenuation specs well within engineering tolerances, and no user of
`ellip` for actual filter design would ever notice the ~1e-6 discrepancy. The
upstream package simply picked `fminbnd`'s default `TolX` when solving the
degree equation — a reasonable speed/accuracy tradeoff for general use, but
not documented as an accuracy guarantee and not tight enough for a 1e-9
cross-check.

For this reason `octave/ncauer.m` is a verbatim copy of upstream
[`ncauer.m`](https://github.com/gnu-octave/octave-signal/blob/main/inst/ncauer.m)
with a single-line patch: `__ellip_ws` now passes
`optimset('TolX', 1e-15, 'MaxIter', 10000)` to `fminbnd`. The generator script
prepends `octave/` to Octave's function path so this local copy shadows the
signal-package version, and stock `ellip()` is called unchanged. With that
change, Octave's independent numerical-solver path and constfilt's
modular-identity path agree to roughly 1e-13 on every test case, and the
reference faithfully exercises the 1e-9 coefficient tolerance.
