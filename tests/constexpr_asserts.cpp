// constexpr_asserts.cpp
//
// Compile-time verification that the batch operator() is fully constexpr.
//
// Each block runs the batch operator inside a constexpr function and
// static_asserts that its output matches the Octave reference within
// tolerance. The file produces no runtime tests (no TEST(), no main());
// successful compilation is the test.
//
// In compile-only mode (CONSTFILT_COMPILE_ONLY=ON) the gtest stubs make this
// file's <gtest/gtest.h> include a no-op, so it links to nothing.

#include <gtest/gtest.h>

#include "butterworth_reference.hpp"
#include "continuous_tf_reference.hpp"
#include "elliptic_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

namespace
{

constexpr double kCoeffTol = 1e-7;
constexpr double kStepTol = 1e-7;

// Wrapper so constexpr functions can return fixed-size arrays.
template <unsigned N> struct Arr
{
    double v[N];
};

// Compile-time absolute-value-within-tolerance comparison.
constexpr bool near_eq(double a, double b, double tol)
{
    double d = a - b;
    return (d < 0 ? -d : d) < tol;
}

// Full-array constexpr comparison — returns false on first out-of-tolerance
// element.
template <unsigned N>
constexpr bool array_near_eq(const double (&a)[N], const double (&b)[N],
                             double tol)
{
    for (unsigned i = 0; i < N; ++i)
    {
        double d = a[i] - b[i];
        if ((d < 0 ? -d : d) >= tol)
        {
            return false;
        }
    }
    return true;
}

// =============================================================================
// Butterworth low-pass (ZOH), N=2, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> bw_lp2_step()
{
    constfilt::Butterworth<double, 2> f(100.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwLp2Step = bw_lp2_step();
static_assert(array_near_eq(kBwLp2Step.v, bw_ref::case_2_100Hz_1000Hz::step,
                            kStepTol),
              "Butterworth LP N=2: batch step disagrees with Octave");

constexpr Arr<256> bw_lp2_impulse()
{
    constfilt::Butterworth<double, 2> f(100.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwLp2Impulse = bw_lp2_impulse();
static_assert(array_near_eq(kBwLp2Impulse.v,
                            bw_ref::case_2_100Hz_1000Hz::impulse, kStepTol),
              "Butterworth LP N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> bw_lp2_chirp()
{
    using Ref = bw_ref::case_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2> f(100.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kBwLp2Chirp = bw_lp2_chirp();
static_assert(array_near_eq(kBwLp2Chirp.v, bw_ref::case_2_100Hz_1000Hz::chirp,
                            kStepTol),
              "Butterworth LP N=2: batch chirp disagrees with Octave");

// =============================================================================
// Butterworth high-pass (ZOH), N=2, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> bw_hp2_step()
{
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwHp2Step = bw_hp2_step();
static_assert(array_near_eq(kBwHp2Step.v, bw_ref::case_hp_2_100Hz_1000Hz::step,
                            kStepTol),
              "Butterworth HP N=2: batch step disagrees with Octave");

constexpr Arr<256> bw_hp2_impulse()
{
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwHp2Impulse = bw_hp2_impulse();
static_assert(array_near_eq(kBwHp2Impulse.v,
                            bw_ref::case_hp_2_100Hz_1000Hz::impulse, kStepTol),
              "Butterworth HP N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> bw_hp2_chirp()
{
    using Ref = bw_ref::case_hp_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kBwHp2Chirp = bw_hp2_chirp();
static_assert(array_near_eq(kBwHp2Chirp.v,
                            bw_ref::case_hp_2_100Hz_1000Hz::chirp, kStepTol),
              "Butterworth HP N=2: batch chirp disagrees with Octave");

// =============================================================================
// Butterworth low-pass (Matched-Z), N=2, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> bw_mz2_step()
{
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> f(100.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwMz2Step = bw_mz2_step();
static_assert(array_near_eq(kBwMz2Step.v, bw_ref::case_mz_2_100Hz_1000Hz::step,
                            kStepTol),
              "Butterworth MZ N=2: batch step disagrees with Octave");

constexpr Arr<256> bw_mz2_impulse()
{
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> f(100.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwMz2Impulse = bw_mz2_impulse();
static_assert(array_near_eq(kBwMz2Impulse.v,
                            bw_ref::case_mz_2_100Hz_1000Hz::impulse, kStepTol),
              "Butterworth MZ N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> bw_mz2_chirp()
{
    using Ref = bw_ref::case_mz_2_100Hz_1000Hz;
    constfilt::Butterworth<double, 2, constfilt::MatchedZ> f(100.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kBwMz2Chirp = bw_mz2_chirp();
static_assert(array_near_eq(kBwMz2Chirp.v,
                            bw_ref::case_mz_2_100Hz_1000Hz::chirp, kStepTol),
              "Butterworth MZ N=2: batch chirp disagrees with Octave");

// =============================================================================
// Elliptic low-pass (ZOH), N=2, Rp=0.5dB, Rs=40dB, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> el_lp2_step()
{
    constfilt::Elliptic<double, 2> f(100.0, 0.5, 40.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kElLp2Step = el_lp2_step();
static_assert(array_near_eq(kElLp2Step.v,
                            el_ref::lp_2_5rp_40rs_100Hz_1000Hz::step, kStepTol),
              "Elliptic LP N=2: batch step disagrees with Octave");

constexpr Arr<256> el_lp2_impulse()
{
    constfilt::Elliptic<double, 2> f(100.0, 0.5, 40.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kElLp2Impulse = el_lp2_impulse();
static_assert(array_near_eq(kElLp2Impulse.v,
                            el_ref::lp_2_5rp_40rs_100Hz_1000Hz::impulse,
                            kStepTol),
              "Elliptic LP N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> el_lp2_chirp()
{
    using Ref = el_ref::lp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2> f(100.0, 0.5, 40.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kElLp2Chirp = el_lp2_chirp();
static_assert(array_near_eq(kElLp2Chirp.v,
                            el_ref::lp_2_5rp_40rs_100Hz_1000Hz::chirp,
                            kStepTol),
              "Elliptic LP N=2: batch chirp disagrees with Octave");

// =============================================================================
// Elliptic high-pass (ZOH), N=2, Rp=0.5dB, Rs=40dB, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> el_hp2_step()
{
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 0.5, 40.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kElHp2Step = el_hp2_step();
static_assert(array_near_eq(kElHp2Step.v,
                            el_ref::hp_2_5rp_40rs_100Hz_1000Hz::step, kStepTol),
              "Elliptic HP N=2: batch step disagrees with Octave");

constexpr Arr<256> el_hp2_impulse()
{
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 0.5, 40.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kElHp2Impulse = el_hp2_impulse();
static_assert(array_near_eq(kElHp2Impulse.v,
                            el_ref::hp_2_5rp_40rs_100Hz_1000Hz::impulse,
                            kStepTol),
              "Elliptic HP N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> el_hp2_chirp()
{
    using Ref = el_ref::hp_2_5rp_40rs_100Hz_1000Hz;
    constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass> f(
        100.0, 0.5, 40.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kElHp2Chirp = el_hp2_chirp();
static_assert(array_near_eq(kElHp2Chirp.v,
                            el_ref::hp_2_5rp_40rs_100Hz_1000Hz::chirp,
                            kStepTol),
              "Elliptic HP N=2: batch chirp disagrees with Octave");

// =============================================================================
// AnalogFilter (ZOH), H(s) = 1/(s^2+3s+2), fs=10Hz
// =============================================================================
constexpr Arr<256> af_zoh_step()
{
    using Ref = ctf_ref::case_2_zoh_fs10;
    constfilt::AnalogFilter<double, 2u> f(Ref::b_s, Ref::a_s,
                                          Ref::sample_rate_hz);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kAfZohStep = af_zoh_step();
static_assert(array_near_eq(kAfZohStep.v, ctf_ref::case_2_zoh_fs10::step,
                            kStepTol),
              "AnalogFilter ZOH N=2: batch step disagrees with Octave");

constexpr Arr<256> af_zoh_impulse()
{
    using Ref = ctf_ref::case_2_zoh_fs10;
    constfilt::AnalogFilter<double, 2u> f(Ref::b_s, Ref::a_s,
                                          Ref::sample_rate_hz);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kAfZohImpulse = af_zoh_impulse();
static_assert(array_near_eq(kAfZohImpulse.v, ctf_ref::case_2_zoh_fs10::impulse,
                            kStepTol),
              "AnalogFilter ZOH N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> af_zoh_chirp()
{
    using Ref = ctf_ref::case_2_zoh_fs10;
    constfilt::AnalogFilter<double, 2u> f(Ref::b_s, Ref::a_s,
                                          Ref::sample_rate_hz);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kAfZohChirp = af_zoh_chirp();
static_assert(array_near_eq(kAfZohChirp.v, ctf_ref::case_2_zoh_fs10::chirp,
                            kStepTol),
              "AnalogFilter ZOH N=2: batch chirp disagrees with Octave");

// =============================================================================
// AnalogFilter (Matched-Z), H(s) = 1/(s^2+3s+2), fs=10Hz
// =============================================================================
constexpr Arr<256> af_mz_step()
{
    using Ref = ctf_ref::case_4_mz_fs10;
    constfilt::AnalogFilter<double, 2u, constfilt::MatchedZ> f(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kAfMzStep = af_mz_step();
static_assert(array_near_eq(kAfMzStep.v, ctf_ref::case_4_mz_fs10::step,
                            kStepTol),
              "AnalogFilter MZ N=2: batch step disagrees with Octave");

constexpr Arr<256> af_mz_impulse()
{
    using Ref = ctf_ref::case_4_mz_fs10;
    constfilt::AnalogFilter<double, 2u, constfilt::MatchedZ> f(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kAfMzImpulse = af_mz_impulse();
static_assert(array_near_eq(kAfMzImpulse.v, ctf_ref::case_4_mz_fs10::impulse,
                            kStepTol),
              "AnalogFilter MZ N=2: batch impulse disagrees with Octave");

constexpr Arr<4096> af_mz_chirp()
{
    using Ref = ctf_ref::case_4_mz_fs10;
    constfilt::AnalogFilter<double, 2u, constfilt::MatchedZ> f(
        Ref::b_s, Ref::a_s, Ref::sample_rate_hz);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kAfMzChirp = af_mz_chirp();
static_assert(array_near_eq(kAfMzChirp.v, ctf_ref::case_4_mz_fs10::chirp,
                            kStepTol),
              "AnalogFilter MZ N=2: batch chirp disagrees with Octave");

// =============================================================================
// Butterworth low-pass (ZOH), N=8, fc=100Hz, fs=1000Hz
// =============================================================================
constexpr Arr<256> bw_lp8_step()
{
    constfilt::Butterworth<double, 8> f(100.0, 1000.0);
    double in[256]{};
    for (unsigned i = 0; i < 256; ++i)
        in[i] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwLp8Step = bw_lp8_step();
static_assert(array_near_eq(kBwLp8Step.v, bw_ref::case_8_100Hz_1000Hz::step,
                            kStepTol),
              "Butterworth LP N=8: batch step disagrees with Octave");

constexpr Arr<256> bw_lp8_impulse()
{
    constfilt::Butterworth<double, 8> f(100.0, 1000.0);
    double in[256]{};
    in[0] = 1.0;
    Arr<256> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<256> kBwLp8Impulse = bw_lp8_impulse();
static_assert(array_near_eq(kBwLp8Impulse.v,
                            bw_ref::case_8_100Hz_1000Hz::impulse, kStepTol),
              "Butterworth LP N=8: batch impulse disagrees with Octave");

constexpr Arr<4096> bw_lp8_chirp()
{
    using Ref = bw_ref::case_8_100Hz_1000Hz;
    constfilt::Butterworth<double, 8> f(100.0, 1000.0);
    double in[4096]{};
    for (unsigned i = 0; i < 4096; ++i)
        in[i] = Ref::chirp_in[i];
    Arr<4096> out{};
    f(in, out.v);
    return out;
}
constexpr Arr<4096> kBwLp8Chirp = bw_lp8_chirp();
static_assert(array_near_eq(kBwLp8Chirp.v, bw_ref::case_8_100Hz_1000Hz::chirp,
                            kStepTol),
              "Butterworth LP N=8: batch chirp disagrees with Octave");

// =============================================================================
// Coefficient spot-checks: prove `static constexpr` filter construction
// produces the right b/a at compile time.
// =============================================================================
namespace
{
static constexpr constfilt::Butterworth<double, 2> kBwLp2(100.0, 1000.0);
static_assert(near_eq(kBwLp2.coeffs_b()[1], bw_ref::case_2_100Hz_1000Hz::b[1],
                      kCoeffTol),
              "Butterworth LP N=2: constexpr ctor b[1] disagrees");
static_assert(near_eq(kBwLp2.coeffs_a()[1], bw_ref::case_2_100Hz_1000Hz::a[1],
                      kCoeffTol),
              "Butterworth LP N=2: constexpr ctor a[1] disagrees");

static constexpr constfilt::Elliptic<double, 2> kElLp2(100.0, 0.5, 40.0,
                                                       1000.0);
static_assert(near_eq(kElLp2.coeffs_b()[0],
                      el_ref::lp_2_5rp_40rs_100Hz_1000Hz::b[0], kCoeffTol),
              "Elliptic LP N=2: constexpr ctor b[0] disagrees");
static_assert(near_eq(kElLp2.coeffs_a()[1],
                      el_ref::lp_2_5rp_40rs_100Hz_1000Hz::a[1], kCoeffTol),
              "Elliptic LP N=2: constexpr ctor a[1] disagrees");
} // namespace

} // namespace
