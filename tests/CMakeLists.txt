macro(create_test target)
    add_test(NAME ${target} COMMAND ${target})
endmacro()

function(add_constfilt_test name source)
    add_executable(${name}.test ${source})
    target_include_directories(${name}.test PRIVATE
        ${PROJECT_SOURCE_DIR}/test_dependencies
        ${PROJECT_SOURCE_DIR})
    target_link_libraries(${name}.test constfilt gtest_main constfilt_warnings)
    create_test(${name}.test)
endfunction()

add_constfilt_test(filter      filter.test.cpp)
add_constfilt_test(discretize  discretize.test.cpp)
add_constfilt_test(butterworth butterworth.test.cpp)
