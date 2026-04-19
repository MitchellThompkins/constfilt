THIS_DIR    := $(shell pwd)
UID          = $(shell id -u)
GID          = $(shell id -g)

BUILD_PREFIX ?= $(THIS_DIR)/build
BUILD_TOOL   ?= make
BUILD_FILE    = Makefile
CMAKE_GENERATOR = "Unix Makefiles"
JOB_FLAG     := -j 4

INSTALL_PREFIX ?= $(THIS_DIR)/build

ifneq "$(CC)" ""
    CMAKE_OPTIONS += -DCMAKE_C_COMPILER=$(CC)
endif
ifneq "$(CXX)" ""
    CMAKE_OPTIONS += -DCMAKE_CXX_COMPILER=$(CXX)
endif

ifeq "$(CMAKE_OPTIONS)" ""
    CMAKE_OPTIONS := -G $(CMAKE_GENERATOR) -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX)
else
    $(shell rm -f $(BUILD_PREFIX)/$(BUILD_FILE))
endif

.PHONY: build
build: $(BUILD_PREFIX)/$(BUILD_FILE)
	@set -o xtrace; \
	export CTEST_OUTPUT_ON_FAILURE=1; \
	cmake --build $(BUILD_PREFIX) --target all -- $(JOB_FLAG) ${a}; \

.PHONY: test
test: $(BUILD_PREFIX)/$(BUILD_FILE)
	@set -o xtrace; \
	export CTEST_OUTPUT_ON_FAILURE=1; \
	cmake --build $(BUILD_PREFIX) --target all -- $(JOB_FLAG); \
	ctest --test-dir $(BUILD_PREFIX) -j$$(getconf _NPROCESSORS_ONLN); \

$(BUILD_PREFIX)/$(BUILD_FILE):
	mkdir -p $(BUILD_PREFIX)
	touch -c $@
	ln -sf $(BUILD_PREFIX)/compile_commands.json compile_commands.json; \
	cd $(BUILD_PREFIX) && \
	cmake .. $(CMAKE_OPTIONS); \

.PHONY: deps
deps:
	./scripts/fetch_consteig.sh
	./scripts/fetch_gcem.sh

.PHONY: generate-reference
generate-reference:
	octave --no-gui octave/generate_butterworth_tests.m

.PHONY: remove
remove:
	rm -rf build/ build-gcc/ build-clang/

.PHONY: format
format:
	find . \( -path "./test_dependencies/googletest" \
		-o -path "./include/constfilt/dependencies" \
		-o -path "./build*" \) -prune \
		-o -type f \( -name "*.hpp" -o -name "*.cpp" \) -print \
		| xargs clang-format -i

.PHONY: check-format
check-format:
	find . \( -path "./test_dependencies/googletest" \
		-o -path "./include/constfilt/dependencies" \
		-o -path "./build*" \) -prune \
		-o -type f \( -name "*.hpp" -o -name "*.cpp" \) -print \
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

%: $(BUILD_PREFIX)/$(BUILD_FILE)
	@set -o xtrace; \
	export CTEST_OUTPUT_ON_FAILURE=1; \
	cmake --build $(BUILD_PREFIX) --target $@ -- $(JOB_FLAG) ${a}; \

Makefile:
	;
