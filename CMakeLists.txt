cmake_minimum_required(VERSION 3.13)

project(constfilt)

################################################################################
# Top-level-only settings
################################################################################
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

    option(CONSTFILT_BUILD_TESTS "Build constfilt tests" ON)

    add_library(constfilt_warnings INTERFACE)
    target_compile_options(constfilt_warnings INTERFACE
        -Wall
        -Wextra
        -pedantic
        -Werror
        -Wmisleading-indentation
        -Wshadow
        -Wconversion
        -Wsign-conversion
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Wundef
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(constfilt_warnings INTERFACE
            -Wduplicated-cond
            -Wduplicated-branches
            -Wlogical-op
        )
    endif()
else()
    option(CONSTFILT_BUILD_TESTS "Build constfilt tests" OFF)
endif()

################################################################################
# consteig dependency
################################################################################
if(NOT TARGET consteig)
    include(FetchContent)
    FetchContent_Declare(
        consteig
        GIT_REPOSITORY https://github.com/MitchellThompkins/consteig
        GIT_TAG        1.0.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(consteig)
endif()

################################################################################
# constfilt — header-only INTERFACE library
################################################################################
add_library(constfilt INTERFACE)
target_include_directories(constfilt INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(constfilt INTERFACE consteig::consteig)
target_compile_features(constfilt INTERFACE cxx_std_17)

################################################################################
# Tests
################################################################################
if(CONSTFILT_BUILD_TESTS)
    add_subdirectory(test_dependencies)
    enable_testing()
    add_subdirectory(tests)
endif()
