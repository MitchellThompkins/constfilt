#include <gtest/gtest.h>

#include "elliptic_reference.hpp"
#include "test_tools.hpp"
#include <constfilt/constfilt.hpp>

using namespace el_ref;

// =============================================================================
// Elliptic low-pass, ZOH discretization
// =============================================================================

FULL_MATRIX(EllipticLP, lp_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::ZOH, constfilt::LowPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticLP, lp_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::ZOH, constfilt::LowPass,
                                false>(200.0, 1.0, 60.0, 4000.0))

// =============================================================================
// Elliptic high-pass, ZOH discretization
// =============================================================================

FULL_MATRIX(EllipticHP, hp_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::ZOH, constfilt::HighPass,
                                false>(100.0, 0.5, 40.0, 1000.0))

FULL_MATRIX(EllipticHP, hp_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::ZOH, constfilt::HighPass,
                                false>(200.0, 1.0, 60.0, 4000.0))

// =============================================================================
// Elliptic + Matched-Z discretization
// =============================================================================

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::MatchedZ,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPMatchedZ, lp_mz_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::MatchedZ,
                                constfilt::LowPass, false>(200.0, 1.0, 60.0,
                                                           4000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::MatchedZ,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPMatchedZ, hp_mz_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::MatchedZ,
                                constfilt::HighPass, false>(200.0, 1.0, 60.0,
                                                            4000.0))

// =============================================================================
// Elliptic low-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(EllipticLPTustin, lp_tu_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::TustinNW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustin, lp_tu_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::TustinNW,
                                constfilt::LowPass, false>(200.0, 1.0, 60.0,
                                                           4000.0))

// =============================================================================
// Elliptic high-pass, TustinNW discretization
// =============================================================================

FULL_MATRIX(EllipticHPTustin, hp_tu_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::TustinNW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustin, hp_tu_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::TustinNW,
                                constfilt::HighPass, false>(200.0, 1.0, 60.0,
                                                            4000.0))

// =============================================================================
// Elliptic low-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::TustinPW,
                                constfilt::LowPass, false>(100.0, 0.5, 40.0,
                                                           1000.0))

FULL_MATRIX(EllipticLPTustinPW, lp_tupw_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::TustinPW,
                                constfilt::LowPass, false>(200.0, 1.0, 60.0,
                                                           4000.0))

// =============================================================================
// Elliptic high-pass, TustinPW discretization
// =============================================================================

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_2_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 2, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_4_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 4, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_5_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 5, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_6_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 6, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_7_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 7, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_8_5rp_40rs_100Hz_1000Hz,
            constfilt::Elliptic<double, 8, constfilt::TustinPW,
                                constfilt::HighPass, false>(100.0, 0.5, 40.0,
                                                            1000.0))

FULL_MATRIX(EllipticHPTustinPW, hp_tupw_3_10rp_60rs_200Hz_4000Hz,
            constfilt::Elliptic<double, 3, constfilt::TustinPW,
                                constfilt::HighPass, false>(200.0, 1.0, 60.0,
                                                            4000.0))

// =============================================================================
// Elliptic default discretization (must be TustinPW)
// =============================================================================

OUTPUT_MATRIX(EllipticDefault, lp_tupw_2_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 2>(100.0, 0.5, 40.0, 1000.0))

// =============================================================================
// Elliptic SOS low-pass, MatchedZ
// =============================================================================

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_2_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 2, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_4_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 4, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_5_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 5, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_6_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 6, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_7_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 7, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_8_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 8, constfilt::MatchedZ>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_MZ_LP, lp_mz_3_10rp_60rs_200Hz_4000Hz,
              constfilt::Elliptic<double, 3, constfilt::MatchedZ>(200.0, 1.0,
                                                                  60.0, 4000.0))

// =============================================================================
// Elliptic SOS high-pass, MatchedZ
// =============================================================================

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_2_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 2, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_4_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 4, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_5_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 5, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_6_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 6, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_7_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 7, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_8_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 8, constfilt::MatchedZ, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_MZ_HP, hp_mz_3_10rp_60rs_200Hz_4000Hz,
    constfilt::Elliptic<double, 3, constfilt::MatchedZ, constfilt::HighPass>(
        200.0, 1.0, 60.0, 4000.0))

// =============================================================================
// Elliptic SOS low-pass, TustinNW
// =============================================================================

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_2_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 2, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_4_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 4, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_5_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 5, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_6_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 6, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_7_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 7, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_8_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 8, constfilt::TustinNW>(100.0, 0.5,
                                                                  40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_TU_LP, lp_tu_3_10rp_60rs_200Hz_4000Hz,
              constfilt::Elliptic<double, 3, constfilt::TustinNW>(200.0, 1.0,
                                                                  60.0, 4000.0))

// =============================================================================
// Elliptic SOS high-pass, TustinNW
// =============================================================================

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_2_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 2, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_4_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 4, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_5_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 5, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_6_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 6, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_7_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 7, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_8_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 8, constfilt::TustinNW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_TU_HP, hp_tu_3_10rp_60rs_200Hz_4000Hz,
    constfilt::Elliptic<double, 3, constfilt::TustinNW, constfilt::HighPass>(
        200.0, 1.0, 60.0, 4000.0))

// =============================================================================
// Elliptic SOS low-pass, TustinPW (default method + default realization)
// =============================================================================

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_2_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 2>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_4_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 4>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_5_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 5>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_6_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 6>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_7_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 7>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_8_5rp_40rs_100Hz_1000Hz,
              constfilt::Elliptic<double, 8>(100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(EllipticSOS_LP, lp_tupw_3_10rp_60rs_200Hz_4000Hz,
              constfilt::Elliptic<double, 3>(200.0, 1.0, 60.0, 4000.0))

// =============================================================================
// Elliptic SOS high-pass, TustinPW
// =============================================================================

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_2_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 2, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_4_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 4, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_5_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 5, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_6_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 6, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_7_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 7, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_8_5rp_40rs_100Hz_1000Hz,
    constfilt::Elliptic<double, 8, constfilt::TustinPW, constfilt::HighPass>(
        100.0, 0.5, 40.0, 1000.0))

OUTPUT_MATRIX(
    EllipticSOS_HP, hp_tupw_3_10rp_60rs_200Hz_4000Hz,
    constfilt::Elliptic<double, 3, constfilt::TustinPW, constfilt::HighPass>(
        200.0, 1.0, 60.0, 4000.0))
