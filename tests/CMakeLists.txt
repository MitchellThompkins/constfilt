if(CONSTFILT_COMPILE_ONLY)
    # OBJECT library: produce .o files only, no link, no test registration.
    # Successful compilation IS the test (any failed static_assert or
    # non-constexpr code path errors out at compile time).
    function(add_constfilt_test name source)
        add_library(${name}.test OBJECT ${source})
        target_include_directories(${name}.test PRIVATE
            ${PROJECT_SOURCE_DIR}/test_dependencies)
        target_link_libraries(${name}.test
            constfilt gtest_main constfilt_warnings)
    endfunction()
else()
    macro(create_test target)
        add_test(NAME ${target} COMMAND ${target})
    endmacro()

    function(add_constfilt_test name source)
        add_executable(${name}.test ${source})
        target_include_directories(${name}.test PRIVATE
            ${PROJECT_SOURCE_DIR}/test_dependencies)
        target_link_libraries(${name}.test constfilt gtest_main constfilt_warnings)
        create_test(${name}.test)
    endfunction()
endif()

add_constfilt_test(filter            filter.test.cpp)
add_constfilt_test(discretize        discretize.test.cpp)
add_constfilt_test(butterworth       butterworth.test.cpp)
add_constfilt_test(analog_filter     analog_filter.test.cpp)
add_constfilt_test(elliptic          elliptic.test.cpp)
add_constfilt_test(constexpr_asserts constexpr_asserts.cpp)
