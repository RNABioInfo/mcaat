#ifndef CAS_GENE_DETECTOR_H
#define CAS_GENE_DETECTOR_H

#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include "sdbg/sdbg.h"
#include "profile.h"

/**
 * @brief Viterbi DP column for incremental scoring
 * 
 * Stores the DP state after processing each amino acid,
 * allowing O(M) extension instead of O(L*M) recalculation.
 */
struct ViterbiColumn {
    std::vector<double> M;      // Match state scores [0..hmm_length]
    std::vector<double> I;      // Insert state scores [0..hmm_length]
    std::vector<double> D;      // Delete state scores [0..hmm_length]
    double best_score;          // Best score seen so far
    int best_hmm_pos;           // HMM position with best score
    int seq_length;             // Number of amino acids processed
    
    ViterbiColumn() : best_score(-1e9), best_hmm_pos(0), seq_length(0) {}
    
    void Initialize(int hmm_length) {
        M.assign(hmm_length + 1, -1e9);
        I.assign(hmm_length + 1, -1e9);
        D.assign(hmm_length + 1, -1e9);
        // Start state
        M[0] = 0.0;
        best_score = 0.0;
        best_hmm_pos = 0;
        seq_length = 0;
    }
};

/**
 * @brief Path information for detected amino acid sequences
 */
struct AminoAcidPathInfo {
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    std::string dna_sequence;
    double total_score;
    int hmm_position;
    bool is_complete;
    
    AminoAcidPathInfo() : total_score(0.0), hmm_position(0), is_complete(false) {}
};

/**
 * @brief CAS gene detector using beam search with incremental Viterbi scoring
 */
class CasGeneDetector {
public:
    CasGeneDetector(SDBG& sdbg);
    CasGeneDetector(SDBG& sdbg, Profile* profile);
    
    /**
     * @brief Beam search for amino acid paths with incremental HMM scoring
     * @param start_node Starting node (containing start codon)
     * @param beam_width Number of paths to keep at each step
     * @param max_depth Maximum traversal depth in bp
     * @param start_codon_offset Offset of start codon within the k-mer
     * @return Vector of amino acid paths with scores
     */
    std::vector<AminoAcidPathInfo> BeamSearchAminoAcids(
        uint64_t start_node,
        int beam_width,
        int max_depth,
        int start_codon_offset
    );
    
private:
    SDBG& sdbg;
    Profile* profile_;
    
    /**
     * @brief Convert DNA codon to amino acid
     */
    char CodonToAminoAcid(const std::string& codon);
    
    /**
     * @brief Get nucleotide sequence for a node
     */
    std::string GetNodeSequence(uint64_t node_id);
    
    /**
     * @brief Initialize Viterbi column for start of sequence
     */
    ViterbiColumn InitializeViterbi();
    
    /**
     * @brief Extend Viterbi by one amino acid - O(M) instead of O(L*M)
     * @param prev Previous DP column
     * @param aa Single amino acid character
     * @return New DP column
     */
    ViterbiColumn ExtendViterbi(const ViterbiColumn& prev, char aa);
    
    /**
     * @brief Core beam traversal with cycle detection and incremental Viterbi
     */
    void BeamTraverse(
        uint64_t start_node,
        int beam_width,
        int max_depth,
        int start_codon_offset,
        std::vector<AminoAcidPathInfo>& result_paths
    );
};

#endif // CAS_GENE_DETECTOR_H
