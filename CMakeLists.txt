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
target_include_directories(constfilt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(constfilt INTERFACE cxx_std_17)
add_library(constfilt::constfilt ALIAS constfilt)

################################################################################
# consteig dependency
################################################################################
if(NOT TARGET consteig::consteig)
    include(FetchContent)
    FetchContent_Declare(
        consteig
        GIT_REPOSITORY https://github.com/MitchellThompkins/consteig
        GIT_TAG        1.0.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(consteig)

    set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
    FetchContent_Declare(
        gcem
        GIT_REPOSITORY https://github.com/MitchellThompkins/gcem
        GIT_TAG        v.1.19.0-consteig.1
    )
    FetchContent_MakeAvailable(gcem)
    unset(CMAKE_POLICY_VERSION_MINIMUM CACHE)

    # Use gcem in no-stdlib (builtin) mode and suppress third-party warnings.
    target_compile_definitions(gcem INTERFACE GCEM_TRAITS_BUILTIN)
    get_target_property(_gcem_inc gcem INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(gcem PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_gcem_inc}"
    )
endif()

target_link_libraries(constfilt INTERFACE
    $<BUILD_INTERFACE:consteig::consteig>
    $<BUILD_INTERFACE:gcem>
)

################################################################################
# Development mode (skipped when consumed via add_subdirectory/FetchContent)
################################################################################
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(CONSTFILT_BUILD_TESTS "Build constfilt tests" OFF)
    option(CONSTFILT_COMPILE_ONLY
           "Build tests as object libraries only (no run, no link to gtest); used for ARM cross-compile constexpr verification."
           OFF)

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

    install(TARGETS constfilt EXPORT constfiltTargets)
    install(DIRECTORY include/constfilt DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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

