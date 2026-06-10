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

<!-- TODO: document iir1 and KFR comparison methodology -->

## Interpreting results

<!-- TODO -->
