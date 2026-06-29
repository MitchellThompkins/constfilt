#include <gtest/gtest.h>

#include "butterworth_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

using namespace bw_ref;

// =============================================================================
// Butterworth low-pass, ZOH discretization
// =============================================================================

FULL_MATRIX(Butterworth, case_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::ZOH>(100.0, 1000.0))

FULL_MATRIX(Butterworth, case_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::ZOH>(500.0, 8000.0))

FULL_MATRIX(Butterworth, case_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::ZOH>(200.0, 4000.0))

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

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ,
                                   constfilt::HighPass>(500.0, 8000.0))

FULL_MATRIX(ButterworthMatchedZHPF, case_mz_hp_3_200Hz_4000Hz,
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

// =============================================================================
// Butterworth high-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinNW,
                                   constfilt::HighPass>(500.0, 8000.0))

FULL_MATRIX(ButterworthTustinHPF, case_tu_hp_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::TustinNW,
                                   constfilt::HighPass>(200.0, 4000.0))

// =============================================================================
// Butterworth low-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(ButterworthTustin, case_tu_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::TustinNW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustin, case_tu_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinNW>(500.0,
                                                                   8000.0))

FULL_MATRIX(ButterworthTustin, case_tu_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::TustinNW>(200.0,
                                                                   4000.0))

// =============================================================================
// Butterworth high-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinPW,
                                   constfilt::HighPass>(500.0, 8000.0))

FULL_MATRIX(ButterworthTustinPWHPF, case_tupw_hp_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::TustinPW,
                                   constfilt::HighPass>(200.0, 4000.0))

// =============================================================================
// Butterworth low-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(ButterworthTustinPW, case_tupw_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_4_100Hz_1000Hz,
            constfilt::Butterworth<double, 4, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_5_100Hz_1000Hz,
            constfilt::Butterworth<double, 5, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_6_100Hz_1000Hz,
            constfilt::Butterworth<double, 6, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_7_100Hz_1000Hz,
            constfilt::Butterworth<double, 7, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_8_100Hz_1000Hz,
            constfilt::Butterworth<double, 8, constfilt::TustinPW>(100.0,
                                                                   1000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_2_500Hz_8000Hz,
            constfilt::Butterworth<double, 2, constfilt::TustinPW>(500.0,
                                                                   8000.0))

FULL_MATRIX(ButterworthTustinPW, case_tupw_3_200Hz_4000Hz,
            constfilt::Butterworth<double, 3, constfilt::TustinPW>(200.0,
                                                                   4000.0))

// =============================================================================
// Butterworth default discretization (must be TustinPW)
// =============================================================================

FULL_MATRIX(ButterworthDefault, case_tupw_1_100Hz_1000Hz,
            constfilt::Butterworth<double, 1>(100.0, 1000.0))

FULL_MATRIX(ButterworthDefault, case_tupw_2_100Hz_1000Hz,
            constfilt::Butterworth<double, 2>(100.0, 1000.0))

FULL_MATRIX(FirstOrderLowPassDefault, case_tupw_1_100Hz_1000Hz,
            constfilt::FirstOrderLowPass<double>(100.0, 1000.0))

FULL_MATRIX(FirstOrderHighPassDefault, case_tupw_hp_1_100Hz_1000Hz,
            constfilt::FirstOrderHighPass<double>(100.0, 1000.0))

// =============================================================================
// Butterworth uniform-zeta low-pass, ZOH discretization
//
// N>=4 with uniform zeta produces repeated complex eigenvalue pairs. consteig's
// QR iteration cannot split the resulting 2x2 blocks, so ZOH is tested only
// for N<=3 here. MatchedZ and Tustin cover N=4,5 below.
// =============================================================================

FULL_MATRIX(ButterworthZetaLP, case_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::ZOH>(100.0, 1000.0,
                                                              0.5))

FULL_MATRIX(ButterworthZetaLP, case_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::ZOH>(100.0, 1000.0,
                                                              0.5))

// =============================================================================
// Butterworth uniform-zeta low-pass, Matched-Z discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZ, case_mz_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::MatchedZ>(100.0,
                                                                   1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta low-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustin, case_tu_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::TustinNW>(100.0,
                                                                   1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta low-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPW, case_tupw_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::TustinPW>(100.0,
                                                                   1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta high-pass, ZOH discretization
//
// Same repeated-eigenvalue limitation as LP: tested only for N<=3.
// =============================================================================

FULL_MATRIX(ButterworthZetaHP, case_hp_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaHP, case_hp_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::ZOH,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta high-pass, Matched-Z discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaMatchedZHP, case_mz_hp_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::MatchedZ,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta high-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinHP, case_tu_hp_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::TustinNW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

// =============================================================================
// Butterworth uniform-zeta high-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_2_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 2, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_3_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 3, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_4_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 4, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_5_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 5, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_6_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 6, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_7_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 7, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))

FULL_MATRIX(ButterworthZetaTustinPWHP, case_tupw_hp_zeta_8_100Hz_1000Hz_z50,
            constfilt::Butterworth<double, 8, constfilt::TustinPW,
                                   constfilt::HighPass>(100.0, 1000.0, 0.5))
