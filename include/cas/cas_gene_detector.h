#ifndef CAS_GENE_DETECTOR_H
#define CAS_GENE_DETECTOR_H

#include <vector>
#include <string>
#include <cstdint>
#include "sdbg/sdbg.h"
#include "profile.h"

// Forward declaration
enum class SearchDirection;

struct ViterbiColumn {
    std::vector<double> M;
    std::vector<double> I;
    std::vector<double> D;
    double E = -1e9;
    double best_score = -1e9;
    int best_hmm_pos = 0;
    int seq_length = 0;
    int band_center = 0;  // Expected HMM position (diagonal)

    void Initialize(int hmm_length) {
        M.assign(hmm_length + 1, -1e9);
        I.assign(hmm_length + 1, -1e9);
        D.assign(hmm_length + 1, -1e9);
        M[0] = 0.0;
        E = -1e9;
        best_score = -1e9;
        best_hmm_pos = 0;
        seq_length = 0;
        band_center = 0;
    }
};

struct AminoAcidPathInfo {
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    std::string dna_sequence;
    double total_score = 0.0;
    double normalized_score = 0.0;
    double evalue = -1.0;  // E-value from Gumbel EVD (-1 = not computed)
    int hmm_position = 0;
    bool is_complete = false;
};

// Result from multi-profile beam search
struct MultiProfileResult {
    std::vector<uint64_t> node_path;
    std::string dna_sequence;
    std::string amino_acids;
    bool is_complete = false;
    
    // Per-profile results: profile_index -> (total_score, normalized_score, hmm_position)
    std::vector<std::tuple<double, double, int>> profile_scores;
};

class CasGeneDetector {
public:
    static constexpr int BAND_WIDTH = 32;  // States above/below diagonal
    static constexpr double DEFAULT_MIN_NORMALIZED_SCORE = 0.05;  // For early termination
    
    CasGeneDetector(SDBG& sdbg);
    CasGeneDetector(SDBG& sdbg, const Profile* profile);
    
    std::vector<AminoAcidPathInfo> BeamSearchAminoAcids(
        uint64_t start_node, int beam_width, int max_depth, int start_codon_offset,
        SearchDirection direction, double min_normalized_score = DEFAULT_MIN_NORMALIZED_SCORE);
    
    // Multi-profile beam search: ONE graph traversal, scores ALL profiles
    // Returns best path(s) with scores for each profile
    // profiles: vector of (profile_ptr, profile_index) pairs
    std::vector<MultiProfileResult> BeamSearchMultiProfile(
        uint64_t start_node, int beam_width, int max_depth, int start_codon_offset,
        const std::vector<std::pair<const Profile*, size_t>>& profiles,
        double min_normalized_score = DEFAULT_MIN_NORMALIZED_SCORE);
    
private:
    SDBG& sdbg;
    const Profile* profile_;
    
    char CodonToAminoAcid(const std::string& codon);
    std::string GetNodeSequence(uint64_t node_id);
    
    ViterbiColumn InitializeViterbi();
    ViterbiColumn ExtendViterbiBanded(const ViterbiColumn& prev, char aa);
    
    // Multi-profile Viterbi helpers
    ViterbiColumn InitializeViterbi(int hmm_length);
    ViterbiColumn ExtendViterbiBanded(const ViterbiColumn& prev, char aa, const Profile* profile);
};

#endif
