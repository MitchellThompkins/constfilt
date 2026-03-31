# Discretization Methods

A filter specification — order $N$, cutoff frequency $f_c$, sample rate $f_s$ — needs
to become a pair of coefficient arrays $b$ and $a$ that the runtime filter can use.
This document describes how that conversion works.

---

## Why continuous time first?

Butterworth filters have a clean analytical definition in continuous time: the
poles of the transfer function lie on the unit circle in the $s$-plane, spaced
to give a maximally flat magnitude response. Working in continuous time first
lets us exploit that structure directly.

The normalized (unit cutoff) Butterworth poles are:

$$s_k = e^{j\theta_k}, \quad \theta_k = \frac{\pi(2k + N - 1)}{2N}, \quad k = 1, \ldots, N$$

These are the left-half-plane roots ($\operatorname{Re}(s_k) < 0$), giving a stable
filter. The polynomial $\prod_k (s - s_k)$ is computed in complex arithmetic;
real coefficients emerge because the poles come in conjugate pairs (and the
real pole at $s = -1$ for odd $N$).

Scaling to the desired cutoff $\omega_c = 2\pi f_c$ is straightforward: the
coefficient of $s^k$ is multiplied by $\omega_c^{N-k}$.

---

## Why state-space?

The discretization methods used here (ZOH and Matched-Z) both work on a
continuous state-space model of the form:

$$\dot{x} = A_c x + B_c u, \quad y = C_c x + D_c u$$

The ZOH formula requires a matrix exponential $e^{A_c T_s}$, computed via
eigendecomposition of $A_c$. consteig provides that eigendecomposition at
compile time, and it operates on matrices, not polynomials. The state-space
$(A, B, C, D)$ struct is just how the system is organized; the $A_c$ matrix
itself is a companion matrix built directly from the denominator polynomial
coefficients (controllable canonical form). It could equally be described as
"construct the companion matrix from the polynomial and feed it to the
eigendecomposition." State-space is the framing, not the reason.

The continuous Butterworth system is expressed in controllable canonical form.
For a denominator $s^N + p_{N-1} s^{N-1} + \cdots + p_0$:

```
A[i][i+1] = 1          for i = 0..N-2
A[N-1][k]  = -p[k]     (last row holds negated, frequency-scaled coefficients)
B[N-1]     = wc^N      (all other B entries are 0)
C[0]       = 1         (all other C entries are 0)
D          = 0
```

---

## Discretization

With a continuous state-space model in hand, the goal is a discrete state-space
model $(A_d, B_d, C_d, D_d)$ valid at sample rate $f_s$.

### Zero-Order Hold (ZOH)

ZOH assumes the input is held constant between samples. Under that assumption,
the discrete model is exact:

$$A_d = e^{A_c T_s}$$

$$B_d = A_c^{-1}(A_d - I)B_c$$

$$C_d = C_c, \quad D_d = D_c$$

The matrix exponential is computed via eigendecomposition:

$$e^A = \sum_i e^{\lambda_i} v_i w_i^T$$

where $v_i$ are the eigenvectors of $A$ (columns of $V$) and $w_i^T$ are the
rows of $V^{-1}$. The computation is done in complex arithmetic; the imaginary
parts cancel because $A$ is real.

$B_d$ is obtained by solving $A_c B_d = (A_d - I)B_c$ via LU decomposition
rather than forming $A_c^{-1}$ explicitly.

### Matched-Z

Matched-Z maps continuous poles directly to discrete poles via $z = e^{s T_s}$
and constructs the transfer function directly, without producing a discrete
state-space model as an intermediate.

For an all-pole continuous system ($D = 0$, no finite zeros), the steps are:

1. Map poles: $A_d = e^{A_c T_s}$, then $a = \operatorname{char\_poly}(A_d)$
2. Place $N-1$ zeros at $z = -1$ (Nyquist frequency), with $b[0] = 0$
3. Choose the gain $K$ to match the continuous DC gain $H_c(0)$:

$$K = H_c(0) \cdot \frac{\sum_k a[k]}{2^{N-1}}, \quad b[k] = K \binom{N-1}{k-1} \text{ for } k = 1,\ldots,N$$

The continuous DC gain is $H_c(0) = D_c - C_c A_c^{-1} B_c$, computed by
solving $A_c x = B_c$ and evaluating $D_c - C_c x$.

---

## From state-space to $b$ and $a$

After ZOH discretization, the discrete state-space model is converted to a
transfer function. The denominator comes from the characteristic polynomial of
$A_d$:

$$a(z) = \det(zI - A_d) = \prod_k (z - \lambda_k)$$

The numerator is recovered from the Markov parameters — the impulse-response
coefficients of the discrete system:

$$h[0] = D_d, \qquad h[k] = C_d A_d^{k-1} B_d \quad \text{for } k = 1, \ldots, N$$

Convolving with the denominator gives the numerator:

$$b[k] = \sum_{j=0}^{k} a[j]\, h[k-j]$$

Matched-Z skips this step because it produces $b$ and $a$ directly.

---

## ZOH vs. Matched-Z

Both methods produce the same poles (the pole mapping $z = e^{s T_s}$ is
identical). They differ in how the zeros are placed and how the gain is set.

ZOH is exact for piecewise-constant inputs and is the standard choice. Matched-Z
preserves the pole locations and matches DC gain, which can give better
frequency-domain accuracy near Nyquist, but is less accurate for step and
impulse inputs. The difference is small when $f_c \ll f_s$ and grows as the
cutoff approaches Nyquist.
