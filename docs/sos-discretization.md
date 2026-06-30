# SOS discretization separability

To avoid the numerical accuracy errors introduced by operating on large
state-space matrices that are produced by high-order filters, a cascade of
second-order sections (SOS) can be used where each section is discretized
independently. That is the default method used by this library to improve
numerical accuracy. However whether the result is equivalent to discretizing the
full-order filter depends on the discretization method.

## Tustin (bilinear)

Tustin applies the algebraic substitution

$$s = \alpha \frac{z - 1}{z + 1}$$

Since this is a direct algebraic substitution, it distributes over products:

$$\text{Tustin}(H_1 \cdot H_2) = \text{Tustin}(H_1) \cdot \text{Tustin}(H_2)$$

SOS Tustin produces an identical filter to full-order Tustin.

## Matched-Z

Matched-Z maps each pole and zero independently:

$$p_s \;\to\; e^{p_s T_s}, \qquad z_s \;\to\; e^{z_s T_s}$$

The pole and zero mapping itself distributes over products. However, the
full-order Matched-Z implementation (following Octave's convention) also pads
the numerator with $n_p - n_z - 1$ zeros at $z = -1$ for strictly-proper
systems, and matches the scalar gain for the full padded numerator at a single
test frequency.

For the SOS cascade to equal the full-order result, these two global quantities
must be distributed consistently across sections rather than applied
independently per section:

- **Padding zeros**: the $n_p - n_z - 1$ total zeros at $z = -1$ are
  distributed greedily across sections (each biquad takes up to 2, each
  real section takes up to 1), so their product equals the full numerator.
- **Test frequency**: a single $\omega_c$ derived from the full analog filter
  is used for all sections, so per-section gains compose to the full-filter
  gain.

With this global distribution the SOS cascade produces an identical filter to
full-order Matched-Z.

## ZOH (zero-order hold)

ZOH is **not** separable over products of transfer functions. The ZOH transfer
function of a continuous-time system $H(s)$ is:

$$H_d(z) = (1 - z^{-1}) \cdot \mathcal{Z}\!\left\{ \mathcal{L}^{-1}\!\left\{ \frac{H(s)}{s} \right\}(nT_s) \right\}$$

The $1/s$ factor means it is the step response of the **entire** system that gets
sampled as a unit. This makes ZOH non-separable over products:

$$\text{ZOH}(H_1 \cdot H_2) \neq \text{ZOH}(H_1) \cdot \text{ZOH}(H_2)$$

### Counterexample

Let $H(s) = 1/(s+1)^2$ with $H_1(s) = H_2(s) = 1/(s+1)$, and let $a = e^{-T_s}$.

The step response of $H$ is $h_\text{step}(t) = 1 - e^{-t} - te^{-t}$. Sampling
and applying the ZOH formula gives:

$$\text{ZOH}(H_1 \cdot H_2) = \frac{\bigl(1 - a(1 + T_s)\bigr)\,z + a(a - 1 + T_s)}{(z - a)^2}$$

Cascading the individually ZOH-discretized first-order sections gives:

$$\text{ZOH}(H_1) \cdot \text{ZOH}(H_2) = \frac{(1 - a)^2}{(z - a)^2}$$

Both have the same denominator but different numerators so the ZOH of a cascaded
set of filters is _not_ equilvalent to the ZOH of a single filter of the same
effective order.

### Consequences for SOS

ZOH SOS is never mathematically equivalent to full-order ZOH, for any filter
type. While it _can_ be done, the error introduced depends on the filter poles,
zeros, and the ratio $T_s \omega_c$.

For Butterworth LP at typical test parameters ($f_c = 100\,\text{Hz}$, $f_s =
1000\,\text{Hz}$) the error happens to fall below $10^{-7}$. That's not bounded
though so at different parameters it may have more error. For Elliptic filters
the finite imaginary zeros produce larger errors.
