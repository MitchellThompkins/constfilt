# SOS discretization separability

When a filter is realized as a cascade of second-order sections (SOS), each
section is discretized independently. Whether the result is equivalent to
discretizing the full-order filter depends on the discretization method.

## Tustin (bilinear)

Tustin applies the algebraic substitution

```
s = alpha * (z - 1) / (z + 1)
```

Since this is a direct substitution, it distributes over products:

```
Tustin(H1 * H2) = Tustin(H1) * Tustin(H2)
```

SOS Tustin produces an identical filter to full-order Tustin.

## Matched-Z

Matched-Z maps each pole and zero independently:

```
p_s  ->  z = exp(p_s * Ts)
z_s  ->  w = exp(z_s * Ts)
```

Because the transform operates pole-by-pole and zero-by-zero, cascading
independently matched-Z-discretized sections preserves the same poles and zeros
as discretizing the full filter. The result is equivalent.

## ZOH (zero-order hold)

ZOH is **not** separable over products of transfer functions. The ZOH transfer
function of a continuous-time system H(s) is:

```
H_d(z) = (1 - z^{-1}) * Z{ L^{-1}{ H(s)/s }(n*Ts) }
```

The 1/s factor means it is the step response of the **entire** system that gets
sampled, not the step responses of the individual sections. Because

```
L^{-1}{ H1(s)*H2(s) / s } != L^{-1}{ H1(s)/s } * L^{-1}{ H2(s)/s }
```

the ZOH of a product is not the product of the ZOHs:

```
ZOH(H1 * H2) != ZOH(H1) * ZOH(H2)
```

### Counterexample

Let H(s) = 1/(s+1)^2, Ts = 0.1:

- ZOH(H) samples the double-pole step response t*exp(-t), producing a
  numerator that depends on (z-1) terms.
- ZOH(1/(s+1)) * ZOH(1/(s+1)) = (1 - exp(-0.1))^2 / (z - exp(-0.1))^2

These are numerically different transfer functions.

### Consequences for SOS

ZOH SOS is never mathematically equivalent to full-order ZOH, for any filter
type. The error depends on the filter poles, zeros, and the ratio Ts*wc. For
Butterworth LP at typical test parameters (fc=100Hz, fs=1000Hz) the error
happens to fall below 1e-7; at different parameters it need not. For Elliptic
filters the finite imaginary zeros produce larger coupling errors (1e-6 to
1e-3 observed) that exceed the test tolerance for the same parameters.

### Summary

ZOH is not appropriate for SOS. Use `SOS = false` when `Method = ZOH`.
Tustin and Matched-Z are separable and work correctly with SOS for any filter
type.
