# SOS discretization separability

When a filter is realized as a cascade of second-order sections (SOS), each
section is discretized independently. Whether the result is equivalent to
discretizing the full-order filter depends on the discretization method.

## Tustin (bilinear)

Tustin applies the algebraic substitution

$$s = \alpha \frac{z - 1}{z + 1}$$

Since this is a direct algebraic substitution, it distributes over products:

$$\text{Tustin}(H_1 \cdot H_2) = \text{Tustin}(H_1) \cdot \text{Tustin}(H_2)$$

SOS Tustin produces an identical filter to full-order Tustin.

## Matched-Z

Matched-Z maps each pole and zero independently:

$$p_s \;\to\; e^{p_s T_s}, \qquad z_s \;\to\; e^{z_s T_s}$$

Because the transform operates pole-by-pole and zero-by-zero, cascading
independently matched-Z-discretized sections preserves the same poles and zeros
as discretizing the full filter. The result is equivalent.

## ZOH (zero-order hold)

ZOH is **not** separable over products of transfer functions. The ZOH transfer
function of a continuous-time system $H(s)$ is:

$$H_d(z) = (1 - z^{-1}) \cdot \mathcal{Z}\!\left\{ \mathcal{L}^{-1}\!\left\{ \frac{H(s)}{s} \right\}(nT_s) \right\}$$

The $1/s$ factor means it is the step response of the **entire** system that gets
sampled, not the step responses of the individual sections. Because

$$\mathcal{L}^{-1}\!\left\{ \frac{H_1(s) H_2(s)}{s} \right\} \neq
\mathcal{L}^{-1}\!\left\{ \frac{H_1(s)}{s} \right\} *
\mathcal{L}^{-1}\!\left\{ \frac{H_2(s)}{s} \right\}$$

the ZOH of a product is not the product of the ZOHs:

$$\text{ZOH}(H_1 \cdot H_2) \neq \text{ZOH}(H_1) \cdot \text{ZOH}(H_2)$$

### Counterexample

Let $H(s) = 1/(s+1)^2$, $T_s = 0.1$:

- $\text{ZOH}(H)$ samples the double-pole step response $t e^{-t}$, producing a
  numerator that depends on $(z-1)$ terms.
- $\text{ZOH}(1/(s+1)) \cdot \text{ZOH}(1/(s+1)) = (1 - e^{-0.1})^2 / (z - e^{-0.1})^2$

These are numerically different transfer functions.

### Consequences for SOS

ZOH SOS is never mathematically equivalent to full-order ZOH, for any filter
type. The error depends on the filter poles, zeros, and the ratio $T_s \omega_c$.
For Butterworth LP at typical test parameters ($f_c = 100\,\text{Hz}$,
$f_s = 1000\,\text{Hz}$) the error happens to fall below $10^{-7}$; at different
parameters it need not. For Elliptic filters the finite imaginary zeros produce
larger coupling errors ($10^{-6}$ to $10^{-3}$ observed) that exceed the test
tolerance for the same parameters.

### Summary

ZOH is not appropriate for SOS. Use `SOS = false` when `Method = ZOH`.
Tustin and Matched-Z are separable and work correctly with SOS for any filter
type.
