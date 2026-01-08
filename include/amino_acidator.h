#ifndef AMINO_ACIDATOR_H
#define AMINO_ACIDATOR_H

#include <vector>
#include <cstdint>
#include <string>
#include <sdbg/sdbg.h>

struct AminoAcidPathInfo {
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    std::vector<double> scores;
    double total_score;
    std::string dna_sequence;
};

class AminoAcidator {
private:
    SDBG& sdbg;
    
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
    AminoAcidator(SDBG& sdbg);
    
    // Main beam search that converts triplets to amino acids
    std::vector<AminoAcidPathInfo> BeamSearchAminoAcids(
        uint64_t start_node,
        int beam_width,
        int max_depth
    );
};

#endif // AMINO_ACIDATOR_H
