# Discretization Methods

A filter specification (order $N$, cutoff frequency $f_c$, sample rate $f_s$)
needs to become a pair of coefficient arrays $b$ and $a$ that the runtime filter
can use. The pipeline has three steps: the filter is designed as a
continuous-time transfer function, converted to a state-space representation,
and then discretized.

## Continuous time

Classical IIR filter types (Butterworth, elliptic, Chebyshev, etc...) have clean
well understood analytical definitions in continuous time. Their poles and zeros are placed in
the s-plane according to closed-form transfer functions that produce the
desired magnitude and phase response. Describing filters in the s-domain is
intuitive and the library here intends to provide that functionality directly
(i.e. by providing a filter type, cut-off, and sample rate).
Working in continuous time first lets the design exploit that structure
directly, then convert to discrete time in a single well-defined step.

## State-space Representation

The discretization methods used here (ZOH, Matched-Z, and Tustin) all work on a
continuous state-space model of the form:

$$\dot{x} = A_c x + B_c u, \quad y = C_c x + D_c u$$

where $x \in \mathbb{R}^N$ is the state vector (internal memory of the system),
$u$ is the scalar input, and $y$ is the scalar output. $A_c$ is the system
matrix whose eigenvalues are the continuous-time poles (describes the internal
dynamics of the system); $B_c$ maps the input into the state derivatives
(describes how the inputs affect the states); $C_c$ reads the output from the
state (describes how the states relate to the outputs); and $D_c$ is the direct
feedthrough from input to output (describes the direct path from input to
output, bypassing the states), which is zero for strictly proper filters (i.e.
more poles than zeros).

ZOH requires $e^{A_c T_s}$ to compute the full discrete state-space.
Matched-Z maps continuous poles and zeros with scalar exponentials
($z = e^{sT_s}$) rather than forming the full matrix exponential.
Tustin does not use the matrix exponential. It operates on
$(A_c, B_c, C_c, D_c)$ directly via a matrix inversion. consteig provides
the eigendecomposition of $A_c$ at compile time. ZOH needs the full
decomposition (eigenvalues and eigenvectors) to compute $B_d$ accurately.
Matched-Z only needs the eigenvalues. Tustin needs neither.

The $A_c$ matrix itself is a companion matrix built directly from the
denominator polynomial coefficients in controllable canonical form, which is the
shared first step before either method is applied.

## Transfer function to state-space

This is that shared first step for ZOH, Matched-Z, and Tustin. Given the
continuous-time transfer function that describes the filter (Butterworth,
Elliptic, etc.), the $(A_c, B_c, C_c, D_c)$ matrices are constructed in
controllable canonical form. All three discretization methods then operate on this
same $A_c$.

For an all-pole filter with numerator gain $b_0$:

$$H(s) = \frac{b_0}{s^N + a_{N-1} s^{N-1} + \cdots + a_1 s + a_0}$$

**Controllable canonical form** assigns the denominator coefficients directly to the last row of $A_c$[^1]. Note that some references[^2][^3] use the transposed convention with coefficients in the first row and a sub-diagonal of 1s; both are valid but differ from the form used here:

$$A_c = \begin{bmatrix}
0      & 1      & 0      & \cdots & 0 \\
0      & 0      & 1      & \cdots & 0 \\
\vdots &        &        & \ddots & \vdots \\
0      & 0      & 0      & \cdots & 1 \\
-a_0   & -a_1   & -a_2   & \cdots & -a_{N-1}
\end{bmatrix}$$

$$B_c = \begin{bmatrix} 0 \\ \vdots \\ 0 \\ b_0 \end{bmatrix}, \qquad
C_c = \begin{bmatrix} 1 & 0 & \cdots & 0 \end{bmatrix}, \qquad
D_c = 0$$

The characteristic polynomial of $A_c$ is $\det(sI - A_c) = s^N + a_{N-1}s^{N-1} + \cdots + a_0$, so the poles of $H(s)$ are exactly the eigenvalues of $A_c$.

### Example: 2nd-order all-pole system

The standard 2nd-order all-pole form uses natural frequency $\omega_n$ and damping ratio $\zeta$:

$$H(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}$$

Reading off $a_0 = \omega_n^2$, $a_1 = 2\zeta\omega_n$, $b_0 = \omega_n^2$:

$$A_c = \begin{bmatrix} 0 & 1 \\ -\omega_n^2 & -2\zeta\omega_n \end{bmatrix}, \qquad
B_c = \begin{bmatrix} 0 \\ \omega_n^2 \end{bmatrix}, \qquad
C_c = \begin{bmatrix} 1 & 0 \end{bmatrix}, \qquad D_c = 0$$

All three methods take this $A_c$ as their starting point. ZOH computes
$e^{A_c T_s}$ to get $A_d$ and solves for $B_d$. Matched-Z uses the eigenvalues
of $A_c$ to map poles and evaluates the continuous transfer function at a test
frequency for gain matching. Tustin applies the bilinear substitution directly
to $(A_c, B_c, C_c, D_c)$ via a single matrix inversion.


## Discretization

### Choosing a method

ZOH and Matched-Z both map continuous poles to discrete poles via
$z = e^{s T_s}$. Tustin maps them via the bilinear substitution
$s = \frac{2}{T_s}\frac{z-1}{z+1}$, which places poles at different locations
than the other two.

ZOH asks: given that the input is held constant between samples, what is the
exact discrete equivalent? It derives the discrete state-space model from the
continuous one via the matrix exponential, then extracts $b$ and $a$ from the
Markov parameters. The result is exact for piecewise-constant inputs but does
not explicitly preserve the continuous frequency response shape.

Matched-Z asks: what discrete filter has the same poles and zeros as the
continuous filter, with gain matched at a reference frequency? It maps
continuous poles and zeros to discrete ones via $z = e^{s T_s}$ and sets the
gain by evaluating both filters at a test frequency. The result preserves the
shape of the continuous frequency response.

Tustin asks: what discrete filter results from substituting the bilinear
approximation $s \approx \frac{2}{T_s}\frac{z-1}{z+1}$ into the continuous
transfer function? The substitution is applied directly to the state-space
matrices. The result warps the frequency axis (frequencies near Nyquist are
compressed), but the discrete filter is stable whenever the continuous filter
is stable and preserves the frequency-domain shape up to that warping.

The choice comes down to what fidelity matters. ZOH is exact for
piecewise-constant (sample-and-hold) inputs, making it the natural choice when
time-domain behavior matters: step response, impulse response, control systems,
and signal reconstruction. Matched-Z explicitly places zeros to match the
continuous-time frequency response shape, which can give better stopband
attenuation when the cutoff is a significant fraction of the sample rate.
Tustin is a good choice when frequency-domain shape is important and
prewarping is acceptable. It is the most commonly used method in digital
control and audio signal processing. The difference is small when
$f_c \ll f_s$. If you are unsure, use ZOH (the default).

Given a continuous state-space model with $(A_c, B_c, C_c, D_c)$, the goal is a
discrete state-space model $(A_d, B_d, C_d, D_d)$ valid at sample rate $f_s$.

### Zero-Order Hold (ZOH)

ZOH assumes the input is held constant between samples. Under that assumption,
the discrete model is exact. The eigenvalues of $A_d = e^{A_c T_s}$ are the
discrete poles, related to the continuous poles by $z = e^{s T_s}$.

$$A_d = e^{A_c T_s}$$

$$B_d = A_c^{-1}(A_d - I)B_c$$

$$C_d = C_c, \quad D_d = D_c$$

The matrix exponential is computed via eigendecomposition:

$$e^A = \sum_i e^{\lambda_i} v_i w_i^T$$

where $v_i$ are the eigenvectors of $A$ (columns of $V$) and $w_i^T$ are the
rows of $V^{-1}$. The computation is done in complex arithmetic; the imaginary
parts cancel because $A$ is real.

$B_d$ is obtained by solving $A_c B_d = (A_d - I)B_c$ via LU decomposition
rather than forming $A_c^{-1}$ explicitly. Explicitly inverting a matrix
amplifies numerical errors in its entries; solving the equivalent linear system
via LU decomposition gives the same result with better numerical stability.

#### From state-space to $b$ and $a$

The denominator $a$ comes from expanding the characteristic polynomial of $A_d$.
The eigenvalues of $A_d$ are the discrete poles, so $\prod_k(z - \lambda_k)$
gives the denominator directly. The numerator $b$ requires one more step. Since
$H(z) = b(z)/a(z)$, rearranging gives $b(z) = H(z) \cdot a(z)$, which in the
time domain is a convolution of $a$ with the impulse response $h[k]$. The
impulse response is traced from the state-space matrices by feeding in a unit
impulse and tracking the output at each step; these samples are the Markov
parameters. Convolving $a$ with $h$ gives $b$. This is shown explicitly below.

The characteristic polynomial of $A_d$ is:

$$a(z) = \det(zI - A_d) = \prod_k (z - \lambda_k)$$

The numerator is recovered from the Markov parameters (the impulse-response
coefficients of the discrete system):

$$h[0] = D_d, \qquad h[k] = C_d A_d^{k-1} B_d \quad \text{for } k = 1, \ldots, N$$

Convolving with the denominator gives the numerator:

$$b[k] = \sum_{j=0}^{k} a[j]\, h[k-j]$$

Matched-Z skips this step because it produces $b$ and $a$ directly.

### Matched-Z

Matched-Z maps continuous poles directly to discrete poles via $z = e^{s T_s}$
and constructs the transfer function directly, without needing a complete
discrete state-space model as an intermediate.

Placing poles and zeros from roots alone leaves the overall gain as a free
parameter: specifying roots determines the shape of the frequency response but
not its magnitude, so any scalar multiple of $b$ gives the same pole/zero
locations. Matched-Z must therefore pin the gain explicitly by matching the continuous
and discrete transfer functions at a test frequency $\omega_c$ (chosen to avoid
poles and zeros). ZOH does not have this problem because $b$ and $a$
are derived from the full matrix computation involving $B_d$, which carries the
input scaling through from $B_c$; the gain is implicit in the matrices.

The steps are:

1. Map poles via eigendecomposition: continuous poles $s_k$ become discrete poles $z_k = e^{s_k T_s}$; build the denominator $a$ as $\prod_k (z - z_k)$.
2. Map finite continuous zeros via $z = e^{s T_s}$. For strictly proper systems, pad any remaining missing zeros at $z = -1$ (Nyquist) until the numerator degree reaches $N - 1$.
3. Match the gain by evaluating the continuous and discrete transfer functions at a test frequency $\omega_c$ (chosen to avoid poles and zeros) and scaling $b$ so that $|H_d(e^{j\omega_c T_s})| = |H_c(j\omega_c)|$.


### Tustin (Bilinear)

Tustin substitutes $s = \frac{2}{T_s}\frac{z-1}{z+1}$ into the continuous
state-space equations and applies a coordinate transformation to put the result
in standard discrete form[^4]. With $\alpha = 2 / T_s$:

$$M = I - \frac{1}{\alpha}A_c$$

$$A_d = \left(I + \frac{1}{\alpha}A_c\right) M^{-1}$$

$$B_d = \frac{1}{\alpha}(A_d + I) B_c$$

$$C_d = C_c M^{-1}$$

$$D_d = D_c + \frac{1}{\alpha} C_d B_c$$

$M^{-1}$ is computed once via LU decomposition and reused for both $A_d$ and
$C_d$. Unlike ZOH, no matrix exponential is required. Unlike Matched-Z, no
companion-matrix eigendecomposition is needed for the zeros. The resulting
discrete state-space is then converted to $(b, a)$ via the same characteristic
polynomial and Markov parameter steps used by ZOH.

## References

[^1]: Swarthmore LPSA, [Transfer Function to State Space](https://lpsa.swarthmore.edu/Representations/SysRepTransformations/TF2SS.html). Uses the same convention as this library (coefficients in the last row, super-diagonal of 1s).
[^2]: R. Murray, [State Feedback](https://www.cds.caltech.edu/~murray/amwiki/State_Feedback.html), CDS 110b Lecture Notes, Caltech.
[^3]: MathWorks, [`tf2ss`](https://www.mathworks.com/help/signal/ref/tf2ss.html).
[^4]: DSP Stack Exchange, [Bilinear transformation of continuous-time state-space system](https://dsp.stackexchange.com/questions/45042/bilinear-transformation-of-continuous-time-state-space-system).
