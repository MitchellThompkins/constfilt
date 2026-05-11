# Building and Testing

constfilt is header-only — there is nothing to build to *use* it. This page
covers building and running the test suite.

## Dev container (preferred)

All development happens inside the consteig dev container (GCC 15 / Clang 21).
The image is pulled automatically on first use. CI runs in the same container.

Common Makefile entry points:

```sh
make container.make.test.gcc        # configure, build, and run tests with GCC
make container.make.test.clang      # ... with Clang
make container.make.build.gcc       # build without running tests
make container.make.build.clang
make container.make.check-format    # clang-format dry-run (enforced by CI)
make container.make.format          # apply clang-format in-place
```

`container.make.<target>` runs `make <target>` inside the container, so any
host Makefile target can be invoked this way.

For an interactive shell:

```sh
make container.start
```

…then inside the container the same top-level targets (`make test.gcc`,
`make test.clang`, etc.) are available without the `container.make.` prefix.

## Host build

If GCC/Clang and CMake are installed locally, the same Makefile targets work
directly on the host:

```sh
make test.gcc
make test.clang
```

This is convenient for quick iteration but less reproducible than the
container — CI uses the container exclusively.

## Adding a test

The test directory is registered via CMake; adding a new test means:

1. Drop a new `*.test.cpp` file into the test directory.
2. Register it with the project's `add_constfilt_test(...)` helper from the
   test `CMakeLists.txt` alongside it.
3. Verify with `make container.make.test.gcc` and
   `make container.make.test.clang`.

See [verification.md](verification.md) for the philosophy and structure of
the existing test layers.

## Regenerating reference data

Numerical references for the regression tests are produced by Octave scripts
checked into the repo and emitted as committed C++ headers. Committing the
generated headers means the test suite runs without an Octave installation
and the references are version-controlled alongside the code that consumes
them.

To regenerate, run the Octave scripts in the `octave/` directory and commit
the updated headers. Each script is self-contained and prints its output
path; running it under `octave` (or `octave-cli`) overwrites the
corresponding header in place. CI does not regenerate references — they are
considered ground truth until a maintainer changes the cases and regenerates
deliberately.
