cmake_minimum_required(VERSION 3.13)

# Version detection: Git tag with fallback to 0.0.0
# conan portfile can inject version via -DCONSTFILT_VERSION=<version>
if(DEFINED CONSTFILT_VERSION)
    set(_constfilt_version "${CONSTFILT_VERSION}")
else()
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE _constfilt_git_tag
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE _git_result
        )
        if(_git_result EQUAL 0)
            string(REGEX REPLACE "^v" "" _constfilt_version "${_constfilt_git_tag}")
        else()
            set(_constfilt_version "0.0.0")
        endif()
    else()
        set(_constfilt_version "0.0.0")
    endif()
endif()

project(constfilt VERSION ${_constfilt_version} LANGUAGES CXX)

include(GNUInstallDirs)

################################################################################
# Public INTERFACE target (the only target external consumers need)
################################################################################
add_library(constfilt INTERFACE)
target_compile_features(constfilt INTERFACE cxx_std_17)
add_library(constfilt::constfilt ALIAS constfilt)

################################################################################
# Bundled dependencies
# Run scripts/fetch_consteig.sh and scripts/fetch_gcem.sh to populate.
################################################################################
set(_deps "${CMAKE_CURRENT_SOURCE_DIR}/include/constfilt/dependencies")

option(CONSTFILT_USE_GCEM
    "Use gcem as the math backend for consteig instead of built-in constexpr math" OFF)
option(CONSTFILT_GCEM_USE_STDLIB
    "Allow gcem to use stdlib headers instead of builtins (only with CONSTFILT_USE_GCEM)" OFF)

# consteig: bundled under include/constfilt/dependencies/consteig/
# Include as SYSTEM to suppress third-party warnings.
add_library(constfilt_consteig INTERFACE)
target_include_directories(constfilt_consteig SYSTEM INTERFACE
    $<BUILD_INTERFACE:${_deps}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
if(CONSTFILT_USE_GCEM)
    target_compile_definitions(constfilt_consteig INTERFACE CONSTEIG_USE_GCEM)
    if(CONSTFILT_GCEM_USE_STDLIB)
        target_compile_definitions(constfilt_consteig INTERFACE CONSTEIG_GCEM_USE_STDLIB)
    endif()
endif()

# gcem: bundled under include/constfilt/dependencies/gcem/
# Include as SYSTEM to suppress third-party warnings.
add_library(constfilt_gcem INTERFACE)
target_include_directories(constfilt_gcem SYSTEM INTERFACE
    $<BUILD_INTERFACE:${_deps}/gcem>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
if(NOT CONSTFILT_GCEM_USE_STDLIB)
    target_compile_definitions(constfilt_gcem INTERFACE GCEM_TRAITS_BUILTIN)
endif()

target_include_directories(constfilt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(constfilt INTERFACE constfilt_consteig constfilt_gcem)

################################################################################
# Development mode (skipped when consumed via add_subdirectory/FetchContent)
################################################################################
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(CONSTFILT_BUILD_TESTS "Build constfilt tests" OFF)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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

    if(CONSTFILT_BUILD_TESTS)
        add_subdirectory(test_dependencies)
        enable_testing()
        add_subdirectory(tests)
    endif()
endif()

################################################################################
# Install rules (package distribution)
################################################################################
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CMakePackageConfigHelpers)

    install(TARGETS constfilt constfilt_consteig constfilt_gcem EXPORT constfiltTargets)
    install(DIRECTORY include/constfilt DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    # Bundle bundled dependencies so installed consumers need no extra find_package.
    install(DIRECTORY ${_deps}/consteig DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES ${_deps}/gcem/gcem.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(DIRECTORY ${_deps}/gcem/gcem_incl DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT constfiltTargets
        FILE constfiltTargets.cmake
        NAMESPACE constfilt::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/constfilt
    )

    configure_package_config_file(
        cmake/constfiltConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/constfiltConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/constfilt
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/constfiltConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/constfiltConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/constfiltConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/constfilt
    )
endif()

