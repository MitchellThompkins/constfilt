THIS_DIR    := $(shell pwd)
UID          = $(shell id -u)
GID          = $(shell id -g)

BUILD_PREFIX ?= $(THIS_DIR)/build
CMAKE_GENERATOR = "Unix Makefiles"
JOB_FLAG     := -j 4

.PHONY: generate-reference
generate-reference:
	octave --no-gui octave/generate_butterworth_tests.m
	octave --no-gui octave/generate_elliptic_tests.m
	octave --no-gui octave/generate_continuous_tf_tests.m

.PHONY: remove
remove:
	rm -rf build/ build-gcc/ build-clang/ profiling/build/

.PHONY: format
format:
	find . \( -path "./test_dependencies/googletest" -o -path "./build*" -o -path "./profiling/build" -o -path "./include/constfilt/vendor/*" -o -path "./.worktree/*" \) -prune \
		-o -type f \( -name "*.hpp" -o -name "*.cpp" \) ! -name "*_reference.hpp" -print \
		| xargs clang-format -i

.PHONY: check-format
check-format:
	find . \( -path "./test_dependencies/googletest" -o -path "./build*" -o -path "./profiling/build" -o -path "./include/constfilt/vendor/*" -o -path "./.worktree/*" \) -prune \
		-o -type f \( -name "*.hpp" -o -name "*.cpp" \) ! -name "*_reference.hpp" -print \
		| xargs clang-format --dry-run --Werror

################################################################################
# Compiler-specific targets
################################################################################

.PHONY: build.gcc
build.gcc:
	cmake -S . -B $(BUILD_PREFIX)-gcc -G $(CMAKE_GENERATOR) \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
		-DCONSTFILT_BUILD_TESTS=ON
	cmake --build $(BUILD_PREFIX)-gcc --target all -- $(JOB_FLAG)

.PHONY: test.gcc
test.gcc: build.gcc
	ctest --test-dir $(BUILD_PREFIX)-gcc -j$$(getconf _NPROCESSORS_ONLN)

.PHONY: build.clang
build.clang:
	cmake -S . -B $(BUILD_PREFIX)-clang -G $(CMAKE_GENERATOR) \
		-DCMAKE_C_COMPILER=clang \
		-DCMAKE_CXX_COMPILER=clang++ \
		-DCONSTFILT_BUILD_TESTS=ON
	cmake --build $(BUILD_PREFIX)-clang --target all -- $(JOB_FLAG)

.PHONY: test.clang
test.clang: build.clang
	ctest --test-dir $(BUILD_PREFIX)-clang -j$$(getconf _NPROCESSORS_ONLN)

################################################################################
# Cross-compile targets (compile-only, static_assert IS the test)
################################################################################

.PHONY: cross.arm-gcc
cross.arm-gcc:
	cmake -S . -B $(BUILD_PREFIX)-arm-gcc -G $(CMAKE_GENERATOR) \
		-DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi-gcc.cmake \
		-DCONSTFILT_BUILD_TESTS=ON \
		-DCONSTFILT_COMPILE_ONLY=ON
	cmake --build $(BUILD_PREFIX)-arm-gcc --target all -- $(JOB_FLAG)

.PHONY: cross.arm-clang
cross.arm-clang:
	cmake -S . -B $(BUILD_PREFIX)-arm-clang -G $(CMAKE_GENERATOR) \
		-DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-none-eabi-clang.cmake \
		-DCONSTFILT_BUILD_TESTS=ON \
		-DCONSTFILT_COMPILE_ONLY=ON
	cmake --build $(BUILD_PREFIX)-arm-clang --target all -- $(JOB_FLAG)

################################################################################
# Profiling targets
################################################################################

.PHONY: generate-accuracy-reference
generate-accuracy-reference:
	octave --no-gui profiling/octave/generate_accuracy_reference.m

.PHONY: profile.gcc
profile.gcc:
	sh profiling/run_profiling.sh g++

.PHONY: profile.clang
profile.clang:
	sh profiling/run_profiling.sh clang++

################################################################################
# Container targets
################################################################################

.PHONY: container.pull
container.pull:
	docker pull ghcr.io/mitchellthompkins/consteig_dev_image:latest

.PHONY: container.start
container.start:
	touch $(THIS_DIR)/.ash_history
	MY_UID=$(UID) MY_GID=$(GID) \
		docker compose -f docker-compose.yml run --rm dev_env 'sh -x'

container.make.%:
	MY_UID=$(UID) MY_GID=$(GID) \
		docker compose -f docker-compose.yml run --rm dev_env \
		'make CC=$(CC) CXX=$(CXX) $*'

Makefile:
	;
