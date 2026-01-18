/**
 * @file buckets.h
 * @brief HMM profile size buckets
 * 
 * Contains size ranges for all HMM profiles based on LENG values.
 * Size range: -15% to +25% of LENG value
 */

#ifndef INCLUDE_BUCKETS_H_
#define INCLUDE_BUCKETS_H_

#include <string>
#include <vector>
#include <cstdint>

namespace HMMProfiles {

struct ProfileSize {
    std::string filename;
    int leng;       // LENG value from HMM file
    int min_aa;     // Minimum amino acids (-15%)
    int max_aa;     // Maximum amino acids (+25%)
    int min_bp;     // Minimum base pairs (min_aa * 3)
    int max_bp;     // Maximum base pairs (max_aa * 3)
};

// Sorted by LENG
static const std::vector<ProfileSize> ALL_PROFILES = {
    {"MerR_1_V-K.hmm", 69, 58, 86, 174, 258},
    {"Cas2_11_CAS-I-II-III-IV-V-VI.hmm", 71, 60, 88, 180, 264},
    {"Cas2_2_CAS-I-II-III-IV-V-VI.hmm", 79, 67, 98, 201, 294},
    {"Cas2_0_I-II-III-V.hmm", 83, 70, 103, 210, 309},
    {"2OG_1_CAS-I-II-III-IV-V-VI.hmm", 94, 79, 117, 237, 351},
    {"Cas6_15_CAS-I-II-III-IV-V-VI.hmm", 112, 95, 140, 285, 420},
    {"Csa3_1_CAS-I-II-III-IV-V-VI.hmm", 118, 100, 147, 300, 441},
    {"Cse2gr11_7_CAS-I-E.hmm", 144, 122, 180, 366, 540},
    {"Cas6_13_CAS-I-II-III-IV-V-VI.hmm", 151, 128, 188, 384, 564},
    {"DEDDh_0_CAS-I-II-III-IV-V-VI.hmm", 160, 136, 200, 408, 600},
    {"WYL_3_CAS-I-II-III-IV-V-VI.hmm", 172, 146, 215, 438, 645},
    {"Cas3_1_I.hmm", 178, 151, 222, 453, 666},
    {"Csm3_0_IIID.hmm", 178, 151, 222, 453, 666},
    {"Csm3gr7_2_CAS-III-A-III-D.hmm", 180, 153, 225, 459, 675},
    {"Csm3_1_IIIAD.hmm", 189, 160, 236, 480, 708},
    {"Csm3gr7_12_CAS-III-A-III-D.hmm", 199, 169, 248, 507, 744},
    {"Cas5_0_CAS-I.hmm", 201, 170, 251, 510, 753},
    {"Cas6e_2_CAS-I-II-III-IV-V-VI.hmm", 202, 171, 252, 513, 756},
    {"Cas3HD_0_CAS-I.hmm", 203, 172, 253, 516, 759},
    {"Cas7b_0_CAS-I-B-I-C.hmm", 254, 215, 317, 645, 951},
    {"Cas1_8_CAS-I-II-III-IV-V-VI.hmm", 283, 240, 353, 720, 1059},
    {"Cas1_0_I-II-III-V.hmm", 315, 267, 393, 801, 1179},
    {"Cas1_0_CAS-I-II-III-IV-V-VI.hmm", 317, 269, 396, 807, 1188},
    {"Cas6_10_CAS-I-II-III-IV-V-VI.hmm", 322, 273, 402, 819, 1206},
    {"TnpB_0_Competitive.hmm", 400, 340, 500, 1020, 1500},
    {"TnpB_1_Competitive.hmm", 409, 347, 511, 1041, 1533},
    {"Cas12f1_2_CAS-V-F.hmm", 429, 364, 536, 1092, 1608},
    {"Cas8e_4_CAS-I-E.hmm", 465, 395, 581, 1185, 1743},
    {"Cas8c_3_CAS-I-C.hmm", 559, 475, 698, 1425, 2094},
};

constexpr int MIN_PROFILE_BP = 174;
constexpr int MAX_PROFILE_BP = 2094;
constexpr int TOTAL_PROFILES = 29;

}  // namespace HMMProfiles

#endif  // INCLUDE_BUCKETS_H_
