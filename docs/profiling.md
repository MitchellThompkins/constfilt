# Profiling

constfilt includes a profiling suite that measures compile-time cost, runtime
throughput, and numerical accuracy. Results are committed under
`profiling/results/` so regressions are visible in code review.

## Running

```sh
make container.make.profile.gcc    # GCC
make container.make.profile.clang  # Clang
```

The Octave accuracy reference must be generated before the first run:

```sh
make generate-accuracy-reference   # requires Octave with signal package
```

## What is measured

<!-- TODO: fill out after profiling effort is complete -->

### Compile-time cost

### Runtime throughput

### Numerical accuracy

## Comparison libraries

**iir1** does runtime coefficient design and DF2T filtering with plain scalar doubles —
a direct apples-to-apples comparison for constfilt's runtime path.

**KFR** is a SIMD-accelerated DSP library. Its `iir_filter<T>` is a type-erased
convenience class: each `apply(dest, src, N)` call heap-allocates a wrapper for the
input pointer, substitutes it into a type-erased expression tree, then dispatches every
SIMD-width chunk through two levels of virtual function pointers. At `-O0` (the default
when no `CMAKE_BUILD_TYPE` is set) none of this inlines and the overhead dominates. At
`-O2`/`-O3` the compiler collapses the expression chain and the gap with scalar
implementations narrows substantially. The KFR runtime numbers reflect abstraction cost
at zero optimisation, not KFR's peak IIR throughput.

## Interpreting results

<!-- TODO -->
