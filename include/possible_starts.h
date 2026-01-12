/**
 * @file possible_starts.h
 * @brief CRISPR-Cas gene size ranges and filtering bins
 * 
 * Contains comprehensive data for all possible CRISPR-Cas start genes with:
 * - Individual gene size ranges (typical, min, max in both AA and BP)
 * - Size bins for fast filtering (SMALL, MEDIUM, LARGE, etc.)
 * 
 * Size range: -15% to +25% of typical size (accounts for natural variation and indels)
 */

#ifndef INCLUDE_POSSIBLE_STARTS_H_
#define INCLUDE_POSSIBLE_STARTS_H_

#include <string>
#include <vector>
#include <cstdint>

namespace CrisprCas {

/**
 * @brief Gene size information for CRISPR-Cas genes
 */
struct GeneSize {
    std::string gene_name;
    int typical_aa;  // Typical length in amino acids
    int min_aa;      // Minimum length in amino acids (-15%)
    int max_aa;      // Maximum length in amino acids (+25%)
    int min_bp;      // Minimum length in base pairs
    int max_bp;      // Maximum length in base pairs
};

/**
 * @brief Size bin for fast gene filtering
 */
struct SizeBin {
    std::string bin_name;        // e.g., "SMALL", "MEDIUM", "LARGE"
    int min_aa;                  // Minimum amino acids in this bin
    int max_aa;                  // Maximum amino acids in this bin
    int min_bp;                  // Minimum base pairs in this bin
    int max_bp;                  // Maximum base pairs in this bin
    std::vector<std::string> genes;  // Genes in this bin
};

/**
 * @brief Distance ranges for ORF finding from CRISPR repeats
 */
enum class DistanceRange {
    SMALL,          // 100-500 bp
    SMALL_MEDIUM,   // 100-1000 bp
    MEDIUM,         // 500-2000 bp
    LARGE,          // 1000-5000 bp
    VERY_LARGE      // 2000-10000 bp
};

/**
 * @brief Distance range boundaries in base pairs
 */
struct DistanceBounds {
    DistanceRange range_type;
    int min_distance;
    int max_distance;
};

/**
 * @brief All CRISPR-Cas gene size data (TABLE1)
 */
static const std::vector<GeneSize> ALL_CAS_GENES = {
    {"Cas1", 300, 255, 375, 765, 1125},
    {"Cas2", 100, 85, 125, 255, 375},
    {"Cas3", 800, 680, 1000, 2040, 3000},
    {"Cas3-Cas2", 900, 765, 1125, 2295, 3375},
    {"Cas3HD", 600, 510, 750, 1530, 2250},
    {"Cas4", 200, 170, 250, 510, 750},
    {"Cas5", 250, 212, 312, 636, 936},
    {"Cas5d", 250, 212, 312, 636, 936},
    {"Cas5f", 250, 212, 312, 636, 936},
    {"Cas5u", 250, 212, 312, 636, 936},
    {"Cas6", 200, 170, 250, 510, 750},
    {"Cas6f", 200, 170, 250, 510, 750},
    {"Cas7", 300, 255, 375, 765, 1125},
    {"Cas7d", 300, 255, 375, 765, 1125},
    {"Cas7f", 300, 255, 375, 765, 1125},
    {"Cas7u", 300, 255, 375, 765, 1125},
    {"Cas8a", 500, 425, 625, 1275, 1875},
    {"Cas8b", 500, 425, 625, 1275, 1875},
    {"Cas8c", 500, 425, 625, 1275, 1875},
    {"Cas8e", 500, 425, 625, 1275, 1875},
    {"Cas8f", 500, 425, 625, 1275, 1875},
    {"Cas8u", 500, 425, 625, 1275, 1875},
    {"Cas9", 1100, 935, 1375, 2805, 4125},
    {"Cas10", 500, 425, 625, 1275, 1875},
    {"Cas10-like", 400, 340, 500, 1020, 1500},
    {"Cas10d", 400, 340, 500, 1020, 1500},
    {"Cas11", 100, 85, 125, 255, 375},
    {"Cas11e", 120, 102, 150, 306, 450},
    {"Cas12a", 1300, 1105, 1625, 3315, 4875},
    {"Cas12b1", 1100, 935, 1375, 2805, 4125},
    {"Cas12b2", 1100, 935, 1375, 2805, 4125},
    {"Cas12c", 1200, 1020, 1500, 3060, 4500},
    {"Cas12d", 1100, 935, 1375, 2805, 4125},
    {"Cas12e", 1000, 850, 1250, 2550, 3750},
    {"Cas12f1", 500, 425, 625, 1275, 1875},
    {"Cas12f2", 500, 425, 625, 1275, 1875},
    {"Cas12f3", 500, 425, 625, 1275, 1875},
    {"Cas12g", 800, 680, 1000, 2040, 3000},
    {"Cas12h", 900, 765, 1125, 2295, 3375},
    {"Cas12i", 1000, 850, 1250, 2550, 3750},
    {"Cas12j", 750, 637, 937, 1911, 2811},
    {"Cas12k", 600, 510, 750, 1530, 2250},
    {"Cas12l", 1100, 935, 1375, 2805, 4125},
    {"Cas12m", 600, 510, 750, 1530, 2250},
    {"Cas13a", 1100, 935, 1375, 2805, 4125},
    {"Cas13b", 1100, 935, 1375, 2805, 4125},
    {"Cas13c", 1000, 850, 1250, 2550, 3750},
    {"Cas13d", 900, 765, 1125, 2295, 3375},
    {"CasR", 350, 297, 437, 891, 1311},
    {"Cmr1", 350, 297, 437, 891, 1311},
    {"Cmr3", 300, 255, 375, 765, 1125},
    {"Cmr4", 250, 212, 312, 636, 936},
    {"Cmr5", 150, 127, 187, 381, 561},
    {"Cmr6", 350, 297, 437, 891, 1311},
    {"Csa5", 150, 127, 187, 381, 561},
    {"Csf1", 400, 340, 500, 1020, 1500},
    {"Csf2", 300, 255, 375, 765, 1125},
    {"Csf3", 250, 212, 312, 636, 936},
    {"Csf3-Csf1", 600, 510, 750, 1530, 2250},
    {"Csf4", 200, 170, 250, 510, 750},
    {"Csm2", 150, 127, 187, 381, 561},
    {"Csm3", 250, 212, 312, 636, 936},
    {"Csm4", 300, 255, 375, 765, 1125},
    {"Csm5", 300, 255, 375, 765, 1125},
    {"Csm6", 400, 340, 500, 1020, 1500},
    {"Csn2", 200, 170, 250, 510, 750},
    {"Csx1", 400, 340, 500, 1020, 1500},
    {"Csx10", 400, 340, 500, 1020, 1500},
    {"Csx19", 150, 127, 187, 381, 561},
    {"Csx27", 150, 127, 187, 381, 561},
    {"Csx28", 200, 170, 250, 510, 750},
    {"CysH", 250, 212, 312, 636, 936},
    {"gRAMP", 1500, 1275, 1875, 3825, 5625},
    {"MerR", 150, 127, 187, 381, 561},
    {"PIWI", 800, 680, 1000, 2040, 3000},
    {"RecD", 600, 510, 750, 1530, 2250},
    {"RT", 450, 382, 562, 1146, 1686},
    {"SSgr11", 300, 255, 375, 765, 1125},
    {"TniQ", 300, 255, 375, 765, 1125},
    {"TnsB", 550, 467, 687, 1401, 2061},
    {"TnsC", 400, 340, 500, 1020, 1500}
};

/**
 * @brief Size bins for fast filtering
 */
static const std::vector<SizeBin> SIZE_BINS = {
    {"SMALL", 85, 187, 255, 561,
        {"Cas2", "Cas11", "Cas11e", "Cmr5", "Csa5", "Csm2", "Csx19", "Csx27", "MerR"}},
    
    {"SMALL-MEDIUM", 170, 312, 510, 936,
        {"Cas4", "Cas5", "Cas5d", "Cas5f", "Cas5u", "Cas6", "Cas6f", "Cmr4", 
         "Csf3", "Csf4", "Csm3", "Csn2", "Csx28", "CysH"}},
    
    {"MEDIUM", 255, 437, 765, 1311,
        {"Cas1", "Cas7", "Cas7d", "Cas7f", "Cas7u", "CasR", "Cmr1", "Cmr3", 
         "Cmr6", "Csf2", "Csm4", "Csm5", "RT", "SSgr11", "TniQ"}},
    
    {"MEDIUM-LARGE", 340, 625, 1020, 1875,
        {"Cas8a", "Cas8b", "Cas8c", "Cas8e", "Cas8f", "Cas8u", "Cas10", 
         "Cas10-like", "Cas10d", "Cas12f1", "Cas12f2", "Cas12f3", "Csf1", 
         "Csm6", "Csx1", "Csx10", "TnsC"}},
    
    {"LARGE", 467, 750, 1401, 2250,
        {"Cas3HD", "Cas12k", "Cas12m", "Csf3-Csf1", "RecD", "TnsB"}},
    
    {"LARGE", 637, 1000, 1911, 3000,
        {"Cas3", "Cas12g", "Cas12j", "PIWI"}},
    
    {"VERY-LARGE", 765, 1375, 2295, 4125,
        {"Cas3-Cas2", "Cas9", "Cas12b1", "Cas12b2", "Cas12d", "Cas12e", 
         "Cas12h", "Cas12i", "Cas12l", "Cas13a", "Cas13b", "Cas13c", "Cas13d"}},
    
    {"HUGE", 1020, 1625, 3060, 4875,
        {"Cas12a", "Cas12c"}},
    
    {"LARGEST", 1275, 1875, 3825, 5625,
        {"gRAMP"}}
};

/**
 * @brief Distance ranges for ORF finding from CRISPR repeats
 */
static const std::vector<DistanceBounds> DISTANCE_RANGES = {
    {DistanceRange::SMALL, 100, 500},
    {DistanceRange::SMALL_MEDIUM, 100, 1000},
    {DistanceRange::MEDIUM, 500, 2000},
    {DistanceRange::LARGE, 1000, 5000},
    {DistanceRange::VERY_LARGE, 2000, 10000}
};

/**
 * @brief Overall ORF range to consider
 */
constexpr int MIN_ORF_AA = 85;
constexpr int MAX_ORF_AA = 1875;
constexpr int MIN_ORF_BP = 255;
constexpr int MAX_ORF_BP = 5625;

/**
 * @brief Total number of known Cas genes
 */
constexpr int TOTAL_CAS_GENES = 80;

}  // namespace CrisprCas

#endif  // INCLUDE_POSSIBLE_STARTS_H_
