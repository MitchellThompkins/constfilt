# Re-use the googletest that ships with consteig.
# consteig_build target was added by the parent CMakeLists, so the
# googletest source is already available on disk.
add_subdirectory(
    ${CMAKE_SOURCE_DIR}/../consteig/test_dependencies/googletest
    gtest_build
    EXCLUDE_FROM_ALL
)
