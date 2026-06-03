if(CONSTFILT_COMPILE_ONLY)
    # Compile-only mode: skip GoogleTest fetch entirely, provide stub headers
    # that turn TEST/EXPECT_*/ASSERT_* into no-ops. Used for ARM cross-compile
    # where there is no runtime to execute tests.
    add_library(gtest_main INTERFACE)
    target_include_directories(gtest_main INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/gtest_stubs)
else()
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
        GIT_SHALLOW    TRUE
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()
