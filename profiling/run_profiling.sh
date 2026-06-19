#!/bin/sh
#
# Compile-time and runtime profiling for constfilt.
#
# HOW IT WORKS
# ------------
# Compile-time: each profiling/compile_time/profile_*.cpp forces static
# constexpr filter construction; compilation time IS the benchmark.
#
# Runtime: builds and runs profiling/bench, which measures operator()
# throughput for constfilt, iir1, and KFR.
#
# The compile-time loop avoids per-file cmake overhead by:
#   1. Running cmake once to produce compile_commands.json (includes flags for
#      all profile_*.cpp targets).
#   2. Extracting the exact compiler invocation from compile_commands.json.
#   3. Invoking the compiler directly per file, timed with /usr/bin/time.
#
# USAGE
# -----
#   ./profiling/run_profiling.sh [compiler=g++] [timeout_sec=300]
#
# REQUIREMENTS
#   g++ or clang++, cmake (>=3.13), /usr/bin/time, timeout, python3, uv (optional)

set -eu

COMPILER="${1:-g++}"
TIMEOUT="${2:-300}"

OS=$(uname -s)
if [ "$OS" = "Darwin" ]; then
    TIME_FLAVOR="bsd"
else
    TIME_FLAVOR="gnu"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="$REPO_ROOT/docs/profiling/results"
COMPILE_DIR="$SCRIPT_DIR/compile_time"

mkdir -p "$RESULTS_DIR"

# Compiler identification
COMPILER_VERSION=$("$COMPILER" --version | head -1)
COMPILER_VER=$("$COMPILER" --version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)

PROBE=$(echo "" | "$COMPILER" -E -dM -x c++ - 2>/dev/null)
if echo "$PROBE" | grep -q "^#define __clang__"; then
    COMPILER_ID="clang"
elif echo "$PROBE" | grep -q "^#define __GNUC__"; then
    COMPILER_ID="gcc"
else
    COMPILER_ID=$(basename "$COMPILER")
    echo "Warning: unrecognized compiler family"
fi

BUILD_DIR="$SCRIPT_DIR/build-$COMPILER_ID"


CT_RESULTS_FILE="$RESULTS_DIR/compile_times_${COMPILER_ID}_${COMPILER_VER}.csv"
RT_RESULTS_FILE="$RESULTS_DIR/runtime_${COMPILER_ID}_${COMPILER_VER}.csv"

echo "Compiler: $COMPILER_VERSION"
echo "Family:   $COMPILER_ID"
echo "Timeout:  ${TIMEOUT}s per file"
echo "Output:   $CT_RESULTS_FILE"
echo "          $RT_RESULTS_FILE"
echo ""

# Step 1: CMake configure
echo "Configuring CMake (to extract compiler flags and build bench)..."
cmake -S "$REPO_ROOT" -B "$BUILD_DIR" \
    -DCMAKE_CXX_COMPILER="$COMPILER" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DCONSTFILT_BUILD_PROFILING=ON \
    -DCONSTFILT_BUILD_TESTS=OFF \
    > /dev/null 2>&1 \
    || { echo "CMake configuration failed"; exit 1; }

# Step 2: Extract compiler flags from compile_commands.json
# Find one profile_*.cpp entry and strip the compiler binary, -c <src>, -o <out>
# to get just the flags (e.g. -std=c++17 -O2 ...).
COMPILE_FLAGS=$(python3 - "$BUILD_DIR/compile_commands.json" "$COMPILE_DIR" <<'EOF'
import json, sys, re

db_path = sys.argv[1]
src_dir = sys.argv[2]

with open(db_path) as f:
    db = json.load(f)

entry = next(e for e in db if src_dir in e["file"])
cmd = entry["command"]

cmd = cmd.split(None, 1)[1]          # strip compiler binary
cmd = re.sub(r'-c\s+\S+', '', cmd)   # strip -c <file>
cmd = re.sub(r'-o\s+\S+', '', cmd)   # strip -o <file>
print(cmd.strip())
EOF
)

echo "Flags from CMake: $COMPILE_FLAGS"
echo ""

# Step 3: OS name for metadata
if [ -f /etc/os-release ]; then
    OS_NAME=$(. /etc/os-release && echo "$NAME")
elif sw_vers > /dev/null 2>&1; then
    OS_NAME="macOS $(sw_vers -productVersion)"
else
    OS_NAME="unknown"
fi

# Step 4: Compile-time profiling loop
printf '# family: %s\n'   "$COMPILER_ID"      > "$CT_RESULTS_FILE"
printf '# compiler: %s\n' "$COMPILER_VERSION" >> "$CT_RESULTS_FILE"
printf '# os: %s\n'       "$OS_NAME"          >> "$CT_RESULTS_FILE"
echo "filter_type,order,method,compile_time_sec,max_rss_kb,exit_code" >> "$CT_RESULTS_FILE"

TOTAL=$(find "$COMPILE_DIR" -name 'profile_*.cpp' | wc -l)
CURRENT=0

for src in "$COMPILE_DIR"/profile_*.cpp; do
    CURRENT=$((CURRENT + 1))
    basename_noext=$(basename "$src" .cpp)

    # Parse filter_type, order, method from profile_<type>_<order>_<method>.cpp
    without_prefix="${basename_noext#profile_}"
    method="${without_prefix##*_}"
    rest="${without_prefix%_*}"
    order="${rest##*_}"
    filter_type="${rest%_*}"

    printf "[%d/%d] %-12s order=%-2s method=%-8s ... " \
        "$CURRENT" "$TOTAL" "$filter_type" "$order" "$method"

    TIME_OUTPUT=$(mktemp)
    EXIT_CODE=0
    if [ "$TIME_FLAVOR" = "gnu" ]; then
        /usr/bin/time -f "%e %M" -o "$TIME_OUTPUT" \
            timeout "$TIMEOUT" \
            "$COMPILER" $COMPILE_FLAGS \
                -c "$src" -o /tmp/profile_out.o \
            2>/dev/null \
            || EXIT_CODE=$?
        WALL_SEC=$(awk 'END{print $1}' "$TIME_OUTPUT")
        MAX_RSS=$(awk 'END{print $2}' "$TIME_OUTPUT")
    else
        /usr/bin/time -l -o "$TIME_OUTPUT" \
            timeout "$TIMEOUT" \
            "$COMPILER" $COMPILE_FLAGS \
                -c "$src" -o /tmp/profile_out.o \
            2>/dev/null \
            || EXIT_CODE=$?
        WALL_SEC=$(awk '/real/{print $1}' "$TIME_OUTPUT")
        MAX_RSS=$(awk '/maximum resident set size/{printf "%.0f\n", $1/1024}' "$TIME_OUTPUT")
    fi

    if [ "$EXIT_CODE" -eq 0 ]; then
        printf "%ss %sKB\n" "$WALL_SEC" "$MAX_RSS"
    elif [ "$EXIT_CODE" -eq 124 ]; then
        WALL_SEC="$TIMEOUT"
        MAX_RSS=0
        printf "TIMEOUT\n"
    else
        printf "FAILED (exit %d) %ss %sKB\n" "$EXIT_CODE" "$WALL_SEC" "$MAX_RSS"
    fi

    echo "$filter_type,$order,$method,$WALL_SEC,$MAX_RSS,$EXIT_CODE" >> "$CT_RESULTS_FILE"
    rm -f "$TIME_OUTPUT"
done

echo ""
echo "Compile-time results: $CT_RESULTS_FILE"

# Step 5: Build bench executables
echo ""
echo "Building bench executables..."
cmake --build "$BUILD_DIR" --target bench_constfilt -- -j "$(getconf _NPROCESSORS_ONLN)" \
    || { echo "bench_constfilt build failed"; exit 1; }
cmake --build "$BUILD_DIR" --target bench_iir1               -- -j "$(getconf _NPROCESSORS_ONLN)" \
    2>/dev/null || true
cmake --build "$BUILD_DIR" --target bench_kfr                -- -j "$(getconf _NPROCESSORS_ONLN)" \
    2>/dev/null || true
cmake --build "$BUILD_DIR" --target bench_accuracy_constfilt -- -j "$(getconf _NPROCESSORS_ONLN)" \
    2>/dev/null || true
cmake --build "$BUILD_DIR" --target bench_accuracy_iir1      -- -j "$(getconf _NPROCESSORS_ONLN)" \
    2>/dev/null || true
cmake --build "$BUILD_DIR" --target bench_accuracy_kfr       -- -j "$(getconf _NPROCESSORS_ONLN)" \
    2>/dev/null || true

# Step 6: Runtime benchmark
echo ""
echo "Running runtime benchmark..."
printf '# family: %s\n'   "$COMPILER_ID"                                 > "$RT_RESULTS_FILE"
printf '# compiler: %s\n' "$COMPILER_VERSION"                           >> "$RT_RESULTS_FILE"
printf '# os: %s\n'       "$OS_NAME"                                    >> "$RT_RESULTS_FILE"
echo "library,filter_type,order,method,ns_per_sample,msa_per_s,dc_gain" >> "$RT_RESULTS_FILE"

"$BUILD_DIR/bin/bench_constfilt" >> "$RT_RESULTS_FILE"
[ -x "$BUILD_DIR/bin/bench_iir1" ] && "$BUILD_DIR/bin/bench_iir1" >> "$RT_RESULTS_FILE"
[ -x "$BUILD_DIR/bin/bench_kfr"  ] && "$BUILD_DIR/bin/bench_kfr"  >> "$RT_RESULTS_FILE"

echo ""
echo "Runtime results: $RT_RESULTS_FILE"

# Step 7: Accuracy benchmark (requires accuracy_reference.hpp)
ACC_RESULTS_FILE="$RESULTS_DIR/accuracy_${COMPILER_ID}_${COMPILER_VER}.csv"

if [ -x "$BUILD_DIR/bin/bench_accuracy_constfilt" ]; then
    echo ""
    echo "Running accuracy benchmark..."
    printf '# family: %s\n'   "$COMPILER_ID"                                       > "$ACC_RESULTS_FILE"
    printf '# compiler: %s\n' "$COMPILER_VERSION"                                 >> "$ACC_RESULTS_FILE"
    printf '# os: %s\n'       "$OS_NAME"                                          >> "$ACC_RESULTS_FILE"
    echo "library,filter_type,order,method,max_b_err,max_a_err,max_step_err"      >> "$ACC_RESULTS_FILE"

    "$BUILD_DIR/bin/bench_accuracy_constfilt" >> "$ACC_RESULTS_FILE"
    [ -x "$BUILD_DIR/bin/bench_accuracy_iir1" ] && "$BUILD_DIR/bin/bench_accuracy_iir1" >> "$ACC_RESULTS_FILE"
    [ -x "$BUILD_DIR/bin/bench_accuracy_kfr"  ] && "$BUILD_DIR/bin/bench_accuracy_kfr"  >> "$ACC_RESULTS_FILE"

    echo "Accuracy results: $ACC_RESULTS_FILE"
else
    echo ""
    echo "bench_accuracy_constfilt not built, skipping accuracy check."
    echo "Generate accuracy_reference.hpp first:"
    echo "  octave --no-gui profiling/octave/generate_accuracy_reference.m"
fi

# Step 8: Analysis
echo ""
echo "=== Compile-time summary ==="
uv run --project "$SCRIPT_DIR" "$SCRIPT_DIR/analyze_results.py" "$CT_RESULTS_FILE" 2>/dev/null || \
    echo "(uv not available, run: uv run profiling/analyze_results.py <csv>)"

echo ""
echo "=== Runtime summary ==="
uv run --project "$SCRIPT_DIR" "$SCRIPT_DIR/analyze_results.py" "$RT_RESULTS_FILE" 2>/dev/null || \
    echo "(uv not available, run: uv run profiling/analyze_results.py <csv>)"

if [ -f "$ACC_RESULTS_FILE" ]; then
    echo ""
    echo "=== Accuracy summary ==="
    uv run --project "$SCRIPT_DIR" "$SCRIPT_DIR/analyze_results.py" "$ACC_RESULTS_FILE" 2>/dev/null || \
        echo "(uv not available, run: uv run profiling/analyze_results.py <csv>)"
fi
