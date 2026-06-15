# Profiling

constfilt includes a profiling suite that measures compile-time cost, runtime
throughput, and numerical accuracy. Results are committed under
`profiling/results/` so that it is easy to track regressions and improvements.

## Running

```sh
make container.make.profile.gcc
make container.make.profile.clang
```

The Octave accuracy reference must be generated before the first run:

```sh
make generate-accuracy-reference
```

## What is measured

### Compile-time cost

Each file under `profiling/compile_time/` forces a `static constexpr` filter
instantiation; the compiler must evaluate all coefficient math at compile time.
Wall time and peak RSS are recorded by the shell.

![Compile times](../profiling/results/compile_times_gcc_15.2.0.png)

Compile time grows roughly exponentially with order. Orders 1 and 2 for
Butterworth and order 2 for Elliptic complete successfully under default
compiler limits; higher orders hit the constexpr recursion depth limit. Tustin
runs slightly slower than ZOH and MatchedZ at each order.

### Runtime throughput

Benchmarks process 5 million samples and report the best time across 5
repetitions, in nanoseconds per sample.

![Runtime throughput](../profiling/results/runtime_gcc_15.2.0.png)

constfilt and iir1 scale linearly with filter order and are within a few
nanoseconds of each other. Elliptic throughput matches Butterworth at the same
order because both use the same Direct Form II Transposed runtime path. See the
Comparison libraries section for notes on the kfr results.

### Numerical accuracy

Coefficient and step response accuracy are measured against an Octave reference
for orders 1 through 12.

![Accuracy](../profiling/results/accuracy_gcc_15.2.0.png)

ZOH accuracy degrades sharply above order 8 because the matrix exponential uses
eigendecomposition, which loses precision as eigenvalues spread at high orders.
MatchedZ and Tustin remain near machine precision through order 12. kfr's
elliptic step responses differ from the Octave reference by roughly 1e-5 to
1e-6 across all tested orders, reflecting a different pole placement algorithm
rather than a numerical precision problem.

## Comparison libraries

[iir1](https://github.com/berndporr/iir1) does runtime coefficient design and
Direct Form II Transposed filtering with plain scalar doubles, covering
Butterworth only.

[KFR](https://github.com/kfrlib/kfr) is a SIMD-accelerated DSP library. Its
`iir_filter<T>` is a type-erased convenience class: each `apply(dest, src, N)`
call heap-allocates a wrapper for the input pointer, substitutes it into a
type-erased expression tree, then dispatches each chunk through two levels of
virtual function pointers. The benchmarks are built without optimization flags,
so none of this inlines and the overhead dominates. At `-O2` or `-O3` the
compiler collapses the expression chain and the gap with scalar implementations
narrows substantially. The kfr runtime numbers reflect abstraction cost at zero
optimization, not kfr's peak IIR throughput.

## Interpreting results

Compile-time cost scales with constexpr instantiation depth, not runtime cost.
A filter that takes 1 second to compile runs at the same ns/sample as one that
compiled in 80 ms.

Runtime numbers are measured without optimization flags, so they reflect
per-sample call overhead rather than peak throughput. constfilt and iir1 are
built under the same conditions and are meant to be an apples-to-apples
comparison.

Accuracy is measured against Octave's Signal Processing Toolbox. constfilt ZOH
diverges from the reference above order 8 due to eigendecomposition
conditioning; constfilt MatchedZ and Tustin remain accurate through order 12.
