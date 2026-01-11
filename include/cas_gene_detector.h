#ifndef CAS_GENE_DETECTOR_H
#define CAS_GENE_DETECTOR_H

#include <vector>
#include <cstdint>
#include <string>
#include <sdbg/sdbg.h>
#include "profile.h"

struct AminoAcidPathInfo {
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    std::vector<double> scores;
    double total_score;
    std::string dna_sequence;
    int hmm_position;  // Current position in HMM profile
};

class CasGeneDetector {
private:
    SDBG& sdbg;
    Profile* profile_;  // Optional HMM profile for scoring
    
    // Convert DNA codon (triplet) to amino acid
    char CodonToAminoAcid(const std::string& codon);
    
    // Extract sequence from a node
    std::string GetNodeSequence(uint64_t node_id);
    
    // Beam search helper
    void BeamTraverse(
        uint64_t start_node,
        int beam_width,
        int max_depth,
        std::vector<AminoAcidPathInfo>& result_paths
    );

public:
    CasGeneDetector(SDBG& sdbg);
    CasGeneDetector(SDBG& sdbg, Profile* profile);
    
    // Main beam search that converts triplets to amino acids
    std::vector<AminoAcidPathInfo> BeamSearchAminoAcids(
        uint64_t start_node,
        int beam_width,
        int max_depth
    );
};

#endif // AMINO_ACIDATOR_H
