#include <gtest/gtest.h>

#include "butterworth_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

using namespace bw_ref;

// =============================================================================
// Butterworth low-pass, ZOH discretization
// =============================================================================

FULL_MATRIX(Butterworth, case_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2>(500.0, 8000.0))

FULL_MATRIX(Butterworth, case_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3>(200.0, 4000.0))

// =============================================================================
// Butterworth high-pass, ZOH discretization
// =============================================================================

FULL_MATRIX(ButterworthHPF, case_hp_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthHPF, case_hp_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::ZOH,
                                   constfilt::HighPass>(500.0, 8000.0))

FULL_MATRIX(ButterworthHPF, case_hp_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::ZOH,
                                   constfilt::HighPass>(200.0, 4000.0))

// =============================================================================
// Butterworth high-pass, Matched-Z discretization
// =============================================================================
//
// TODO(issue #16): HP Butterworth has N zeros at s=0 which map to z=1 under
// matched-Z, but matched_z_discretize hardcodes zeros at z=-1 and ignores
// finite analog zeros. Strip the DISABLED_ prefix once #16 is fixed.

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ,
                                   constfilt::HighPass>(500.0, 8000.0))

FULL_MATRIX(DISABLED_ButterworthMatchedZHPF, case_mz_hp_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::MatchedZ,
                                   constfilt::HighPass>(200.0, 4000.0))

// =============================================================================
// Butterworth low-pass, Matched-Z discretization
// =============================================================================

FULL_MATRIX(ButterworthMatchedZ, case_mz_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::MatchedZ>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ>(500.0,
                                                                   8000.0))

FULL_MATRIX(ButterworthMatchedZ, case_mz_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::MatchedZ>(200.0,
                                                                   4000.0))
