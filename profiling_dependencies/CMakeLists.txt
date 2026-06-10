include(FetchContent)

# KFR requires C++20 for both its internal CMake tests and its headers.
# Setting the global standard here is safe: by this point constfilt and its
# vendored targets are already defined (they won't be re-configured), and the
# profiling OBJECT libs + bench target are defined later — C++17 code compiles
# cleanly under C++20.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# iir1: runtime IIR filter design + execution; MIT license
FetchContent_Declare(
    iir1
    GIT_REPOSITORY https://github.com/berndporr/iir1.git
    GIT_TAG        1.9.4
    GIT_SHALLOW    TRUE
)
set(IIR1_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(IIR1_BUILD_DEMO    OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(iir1)

# KFR: full filter design + SIMD-accelerated execution; GPL v2+
# Requires C++20; bench target is compiled as C++20 to accommodate this.
FetchContent_Declare(
    kfr
    GIT_REPOSITORY https://github.com/kfrlib/kfr.git
    GIT_TAG        7.0.1
    GIT_SHALLOW    TRUE
)
set(KFR_ENABLE_TESTS    OFF CACHE BOOL "" FORCE)
set(KFR_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(kfr)
