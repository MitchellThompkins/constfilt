#ifndef CONSTFILT_GTEST_STUBS_HPP
#define CONSTFILT_GTEST_STUBS_HPP

// gtest stubs for CONSTFILT_COMPILE_ONLY mode (e.g. ARM cross-compile).
//
// Replaces the entire <gtest/gtest.h> with no-op macros. Test bodies still
// compile (so any static_assert or constexpr-evaluation embedded in them
// fires), but EXPECT_*, ASSERT_*, and TEST() generate functions that link to
// nothing. Successful compilation IS the test.
//
// All EXPECT_*/ASSERT_* macros return a Sink that swallows any subsequent
// `<<` streaming so test code like `EXPECT_NEAR(a, b, t) << "msg";` works
// without pulling in any iostream code.

namespace gtest_stub
{
struct Sink
{
    template <typename T> constexpr Sink &operator<<(const T &) { return *this; }
};
inline Sink make_sink() { return {}; }
} // namespace gtest_stub

#define CONSTFILT_GTEST_SINK_ ::gtest_stub::make_sink()

#define TEST(suite, name)                                                      \
    [[maybe_unused]] static void suite##_##name##_stub_body()

#define EXPECT_TRUE(a)               (((void)(a)), CONSTFILT_GTEST_SINK_)
#define EXPECT_FALSE(a)              (((void)(a)), CONSTFILT_GTEST_SINK_)
#define EXPECT_EQ(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_NE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_LT(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_LE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_GT(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_GE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_DOUBLE_EQ(a, b)       (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_FLOAT_EQ(a, b)        (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define EXPECT_NEAR(a, b, t)         (((void)(a)), ((void)(b)), ((void)(t)), CONSTFILT_GTEST_SINK_)

#define ASSERT_TRUE(a)               (((void)(a)), CONSTFILT_GTEST_SINK_)
#define ASSERT_FALSE(a)              (((void)(a)), CONSTFILT_GTEST_SINK_)
#define ASSERT_EQ(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_NE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_LT(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_LE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_GT(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_GE(a, b)              (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_DOUBLE_EQ(a, b)       (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_FLOAT_EQ(a, b)        (((void)(a)), ((void)(b)), CONSTFILT_GTEST_SINK_)
#define ASSERT_NEAR(a, b, t)         (((void)(a)), ((void)(b)), ((void)(t)), CONSTFILT_GTEST_SINK_)

#define FAIL()                       CONSTFILT_GTEST_SINK_
#define SUCCEED()                    CONSTFILT_GTEST_SINK_

#endif // CONSTFILT_GTEST_STUBS_HPP
