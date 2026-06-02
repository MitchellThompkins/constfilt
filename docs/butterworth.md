# Butterworth Filter Coefficient Computation

`Butterworth<T, N, Method, FilterType>` computes a continuous-time Butterworth
transfer function in the s-domain, then discretizes it via the parent class
`AnalogFilter`. All arithmetic is `constexpr`.

The algorithm follows the standard analytic Butterworth pole formula directly:
poles lie on the unit circle in the left half of the s-plane at angles spaced
to give a maximally flat (no ripple) magnitude response. There is no iteration
and no numerical solver; the pole locations are exact closed-form expressions.

## Notation

| Symbol | Meaning |
|--------|---------|
| $N$ | filter order |
| $M$ | $\lfloor N/2 \rfloor$, number of complex-conjugate pole pairs |
| $\omega_c$ | cutoff frequency in rad/s ($2\pi f_c$) |
| $\theta_k$ | angle of the $k$-th pole on the unit circle |
| $p_k$ | $k$-th analog prototype pole (complex) |
| $p[i]$ | coefficient of $s^i$ in the normalized denominator polynomial |

## Step 1: Pole placement

The $N$ poles of the normalized ($\omega_c = 1$) Butterworth prototype lie on the unit
circle in the left half of the s-plane. Their angles are:

$$\theta_k = \frac{\pi(2k + N - 1)}{2N}, \qquad k = 1, \ldots, N$$

Each pole is:

$$p_k = \cos\theta_k + j\sin\theta_k$$

By construction, $\mathrm{Re}(p_k) < 0$ for all $k$ (left half-plane), so the prototype
is stable. For odd $N$, $\theta_N = \pi$ gives the real pole at $p_N = -1$. For
even $N$, all poles come in conjugate pairs, so the denominator polynomial has
real coefficients without any special-casing.

## Step 2: Denominator polynomial

The normalized denominator is the product of all $N$ linear factors $(s - p_k)$.
It is computed by sequential polynomial multiplication, accumulating the product
in ascending-coefficient order, starting from the constant polynomial 1:

```text
for k = 1..N:
    for j = k downto 1:
        poly[j] = poly[j-1] - p_k * poly[j]
    poly[0] = -p_k * poly[0]
```

The arithmetic is carried out in complex type. Because every complex pole either
has a conjugate partner in the set or is real, all imaginary parts of the
coefficients cancel exactly (to within floating-point rounding). The real parts
are extracted at the end and reversed to descending order, giving the standard
Butterworth polynomial of order $N$. For example:

$$N=1: \quad s + 1$$

$$N=2: \quad s^2 + \sqrt{2}\,s + 1$$

$$N=3: \quad s^3 + 2s^2 + 2s + 1$$

## Step 3: Cutoff frequency scaling

The normalized prototype has its $-3\,\mathrm{dB}$ point at 1 rad/s. Scaling to
$\omega_c$ rad/s is done by substituting $s \to s/\omega_c$, which multiplies the
coefficient of $s^k$ in the denominator by $\omega_c^k$. Writing the denominator in
descending order as $a[0]s^N + a[1]s^{N-1} + \cdots + a[N]$:

$$a[k] = p[k] \cdot \omega_c^k, \qquad k = 0, \ldots, N$$

where $p[k]$ are the normalized coefficients in descending order and $p[0] = 1$.

## Step 4: Numerator

The Butterworth lowpass prototype is all-pole (no finite zeros). The transfer
function is:

$$H(s) = \frac{\omega_c^N}{a[0]s^N + a[1]s^{N-1} + \cdots + a[N]}$$

so $b[N] = \omega_c^N$ and $b[k] = 0$ for $k < N$. This gives unity DC gain: $H(0) = 1$.

## Example: 3rd-order state-space construction

The normalized 3rd-order Butterworth denominator is $(s+1)(s^2+s+1) = s^3 + 2s^2 + 2s + 1$.
Scaling to cutoff $\omega_c$ (Steps 3 and 4):

$$H(s) = \frac{\omega_c^3}{s^3 + 2\omega_c s^2 + 2\omega_c^2 s + \omega_c^3}$$

`AnalogFilter` places this into controllable canonical form by reading off
$p_0 = \omega_c^3$, $p_1 = 2\omega_c^2$, $p_2 = 2\omega_c$, $b_0 = \omega_c^3$:

$$A_c = \begin{bmatrix}
0           & 1            & 0          \\
0           & 0            & 1          \\
-\omega_c^3 & -2\omega_c^2 & -2\omega_c
\end{bmatrix}, \qquad
B_c = \begin{bmatrix} 0 \\ 0 \\ \omega_c^3 \end{bmatrix}, \qquad
C_c = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix}, \qquad D_c = 0$$

The three eigenvalues of $A_c$ are the three Butterworth poles (one real at
$-\omega_c$, one conjugate pair), and $C_c(sI - A_c)^{-1}B_c$ recovers $H(s)$ exactly.
See [discretization.md](discretization.md) for the general controllable canonical form and the
discretization step that follows.

## Highpass transform

The highpass variant applies the LP-to-HP s-domain substitution $s \to \omega_c/s$
to the normalized lowpass prototype evaluated at $\omega_c = 1$. In polynomial
coefficient terms this reverses the denominator and applies $\omega_c^k$ scaling:

$$a_\mathrm{hp}[k] = p[N - k] \cdot \omega_c^k, \qquad k = 0, \ldots, N$$

Since $p[0] = 1$, the leading coefficient is $a_\mathrm{hp}[0] = p[N] = 1$,
keeping the denominator monic.

The numerator has $b[0] = 1$ and all other $b[k] = 0$, giving unity
high-frequency gain: $H(\infty) = 1$.

## Test reference

`tests/butterworth_reference.hpp` is regenerated from Octave's stock `butter()`
via `octave/generate_butterworth_tests.m`. The reference is independent of the
constfilt implementation: Octave computes Butterworth coefficients via its own
internal routines, and the constfilt pole-multiply-out path is exercised
separately.

Coefficients and filter responses are checked against the Octave reference
across a range of orders, cutoff-to-sample-rate ratios, lowpass and highpass
variants, and both discretization methods.

See [verification.md](verification.md) for the full tolerance rationale.
