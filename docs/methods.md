# Discretization Methods

constfilt converts continuous-time state-space models to discrete-time transfer
functions at compile time. This document describes each method and the
supporting math behind it.

---

## Zero-Order Hold (ZOH)

**Header:** `constfilt/discretize.hpp` - `zoh_discretize(sys_c, Ts, ZOH{})`

ZOH is the primary discretization method and is what `Butterworth` uses
internally. It assumes the input is held constant between samples - the
zero-order hold interpretation - and produces exact discrete-time poles and
zeros for that assumption.

### Formulas

Given a continuous state-space system $(A_c, B_c, C_c, D_c)$ and sample period
$T_s$:

$$A_d = e^{A_c T_s}$$

$$B_d = A_c^{-1}(A_d - I)B_c$$

$$C_d = C_c, \quad D_d = D_c$$

$C_c$ and $D_c$ pass through unchanged. $B_d$ is computed by solving the linear
system $A_c B_d = (A_d - I)B_c$ via LU decomposition rather than forming the
inverse explicitly.

### Matrix Exponential

$e^A$ is computed via eigendecomposition:

$$e^A = V \operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n}) V^{-1}
     = \sum_i e^{\lambda_i} v_i w_i^T$$

where $V$ is the matrix of eigenvectors (columns), $w_i^T$ is row $i$ of
$V^{-1}$, and the $\lambda_i$ are the (possibly complex) eigenvalues. The sum
is accumulated in complex arithmetic; the imaginary parts cancel for real $A$
and the real part is extracted at the end.

$V^{-1}$ is formed column-by-column by solving $V e_j = \text{col}_j$ via LU
decomposition for each standard basis vector $e_j$.

---

## Matched-Z

**Header:** `constfilt/discretize.hpp` - `matched_z_discretize(sys_c, Ts, MatchedZ{})`

Matched-Z maps continuous poles directly to discrete poles via $z = e^{sT_s}$
and returns a transfer function $(b, a)$ directly rather than a state-space
model. It is intended for **all-pole continuous systems** (no finite
continuous-domain zeros).

### Algorithm

1. $A_d = e^{A_c T_s}$ -- pole mapping (same as ZOH)
2. $a = \operatorname{char\_poly}(A_d)$ -- discrete denominator coefficients
3. $H_c(0) = D_c - C_c A_c^{-1} B_c$ -- continuous DC gain
4. $K = H_c(0) \cdot a(1) / 2^N$ -- scaling constant
5. $b[k] = K \binom{N}{k}$ -- numerator coefficients of $K(z+1)^N$

Step 3 evaluates the continuous transfer function at DC. The steady-state
response of $\dot{x} = Ax + Bu$, $y = Cx + Du$ at $s = 0$ gives
$H(0) = D - CA^{-1}B$, computed by solving $A_c x = B_c$ via LU then
$H_c(0) = D_c - C_c x$.

Step 4 chooses $K$ so that the discrete DC gain matches the continuous one:
$H_d(z=1) = H_c(0)$. The denominator evaluated at $z = 1$ is
$a(1) = \sum_k a[k]$, and the numerator $(z+1)^N$ evaluated at $z = 1$ is $2^N$.

Step 5 places $N$ zeros at $z = -1$, giving maximum attenuation at the Nyquist
frequency and preserving the all-pole character of the original system.

### Characteristic Polynomial

`char_poly(Ad)` returns the monic polynomial $[1, c_1, \ldots, c_n]$ whose
roots are the eigenvalues of $A_d$. It builds the polynomial in complex
arithmetic by multiplying in each eigenvalue $\lambda_k$:

$$p(z) = (z - \lambda_1)(z - \lambda_2) \cdots (z - \lambda_n)$$

The product is accumulated in-place working from highest to lowest degree at
each step. The imaginary parts cancel for real $A_d$; the real part is
extracted at the end.

---

## State-Space to Transfer Function (`ss_to_tf`)

**Header:** `constfilt/discretize.hpp` - `ss_to_tf(sys_d)`

Converts a discrete state-space model to $(b, a)$ form. This is a two-step
pipeline used internally by `Butterworth` after ZOH discretization:

1. **Denominator** - `char_poly(Ad)` gives the monic denominator.
2. **Numerator** - `markov_numerator` computes the numerator from the Markov
   parameters and the denominator.

### Markov Parameters

The Markov parameters are the impulse-response coefficients of the state-space
system:

$$h[0] = D, \qquad h[k] = C A^{k-1} B \quad \text{for } k = 1, \ldots, N$$

They are computed iteratively by maintaining a running power $A^{k-1}$ and
multiplying by $B$ and $C$ at each step.

The numerator is then formed by convolving the Markov parameters with the
denominator:

$$b[k] = \sum_{j=0}^{k} a[j]\, h[k-j]$$

This is equivalent to the $z$-transform identity $B(z) = A(z) H(z)$ truncated
to degree $N$.

---

## Butterworth Lowpass Filter

**Header:** `constfilt/butterworth.hpp` - `Butterworth<T, N>(cutoff_hz, sample_rate_hz)`

`Butterworth` is the only filter type currently provided. It constructs a
degree-$N$ Butterworth lowpass filter from a cutoff frequency and sample rate,
performing all computation at compile time via `constexpr`.

### Design Pipeline

```
Specification (fc, fs)
  -> Normalized Butterworth polynomial coefficients  [butterworth_poly_coeffs]
  -> Continuous state-space in controllable canonical form  [build_continuous_ss]
  -> ZOH discretization  [zoh_discretize]
  -> Transfer function (b, a)  [ss_to_tf]
  -> Filter  [Filter<T, N+1, N+1>]
```

### Normalized Prototype

The normalized Butterworth prototype has poles on the unit circle at:

$$s_k = e^{j\theta_k}, \quad \theta_k = \frac{\pi(2k + N - 1)}{2N}, \quad k = 1, \ldots, N$$

These are the left-half-plane roots ($\operatorname{Re}(s_k) < 0$), which give
a stable filter. The polynomial $\prod_k (s - s_k)$ is accumulated in complex
arithmetic; real coefficients emerge because poles come in conjugate pairs.

### Frequency Scaling

The normalized cutoff is $\omega_c = 1\,\text{rad/s}$. To shift the cutoff to
$\omega_c = 2\pi f_c$, the coefficient of $s^k$ is multiplied by
$\omega_c^{N-k}$. This is applied when constructing the companion matrix below.

### Continuous State-Space (Controllable Canonical Form)

For an $N$-th order all-pole system with denominator
$s^N + p_{N-1} s^{N-1} + \cdots + p_0$, the controllable canonical form is:

```
A[i][i+1] = 1          for i = 0..N-2
A[N-1][k]  = -p[k]     (last row holds negated coefficients)
B[N-1]     = wc^N      (all other B entries are 0)
C[0]       = 1         (all other C entries are 0)
D          = 0
```

where `p[k]` is the frequency-scaled Butterworth polynomial coefficient of
$s^k$.
