#include "cas_gene_detector.h"
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>
#include <cmath>

CasGeneDetector::CasGeneDetector(SDBG& sdbg) : sdbg(sdbg), profile_(nullptr) {}

CasGeneDetector::CasGeneDetector(SDBG& sdbg, Profile* profile) : sdbg(sdbg), profile_(profile) {}

char CasGeneDetector::CodonToAminoAcid(const std::string& codon) {
    static const std::unordered_map<std::string, char> codon_table = {
        // Standard genetic code
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
    
    if (codon.length() != 3) {
        return 'X'; // Unknown
    }
    
    auto it = codon_table.find(codon);
    if (it != codon_table.end()) {
        return it->second;
    }
    return 'X'; // Unknown codon
}

std::string CasGeneDetector::GetNodeSequence(uint64_t node_id) {
    const uint32_t k = sdbg.k();
    std::vector<uint8_t> seq(k);
    sdbg.GetLabel(node_id, seq.data());
    
    // Convert numeric encoding to nucleotides
    std::string result;
    result.reserve(k);
    for (uint32_t i = 0; i < k; ++i) {
        switch (seq[i]) {
            case 1: result += 'A'; break;
            case 2: result += 'C'; break;
            case 3: result += 'G'; break;
            case 4: result += 'T'; break;
            default: result += 'N'; break;
        }
    }
    return result;
}

void CasGeneDetector::BeamTraverse(
    uint64_t start_node,
    int beam_width,
    int max_depth,
    std::vector<AminoAcidPathInfo>& result_paths) {
    
    struct BeamState {
        std::vector<uint64_t> node_path;
        std::vector<std::string> amino_acids;
        std::vector<double> scores;
        double total_score;
        std::string accumulated_sequence;
        uint64_t current_node;
        int depth;
        int hmm_position;  // Current position in HMM profile
    };
    
    auto comp = [](const BeamState& a, const BeamState& b) {
        return a.total_score < b.total_score;
    };
    
    std::vector<BeamState> current_layer;
    
    // Initialize with start node
    BeamState initial;
    initial.node_path.push_back(start_node);
    initial.current_node = start_node;
    initial.total_score = 0.0;  // Start with zero score
    initial.scores.push_back(0.0);
    initial.accumulated_sequence = GetNodeSequence(start_node);
    initial.depth = 0;
    initial.hmm_position = 0;  // Start at beginning of HMM
    initial.depth = 0;
    
    current_layer.push_back(initial);
    
    // Breadth-first expansion with beam pruning
    for (int depth = 0; depth < max_depth; ++depth) {
        std::vector<BeamState> next_layer;
        
        for (const auto& state : current_layer) {
            // Get outgoing edges
            uint64_t outgoings[4];
            int outdegree = sdbg.OutgoingEdges(state.current_node, outgoings);
            
            if (outdegree <= 0) {
                // Dead end - save current path
                AminoAcidPathInfo path_info;
                path_info.amino_acids = state.amino_acids;
                path_info.node_path = state.node_path;
                path_info.scores = state.scores;
                path_info.total_score = state.total_score;
                path_info.dna_sequence = state.accumulated_sequence;
                path_info.hmm_position = state.hmm_position;
                result_paths.push_back(path_info);
                continue;
            }
            
            // Check if HMM is complete
            if (profile_ != nullptr && state.hmm_position >= profile_->GetLength()) {
                // HMM fully aligned, save this path and don't expand further
                AminoAcidPathInfo path_info;
                path_info.amino_acids = state.amino_acids;
                path_info.node_path = state.node_path;
                path_info.scores = state.scores;
                path_info.total_score = state.total_score;
                path_info.dna_sequence = state.accumulated_sequence;
                path_info.hmm_position = state.hmm_position;
                result_paths.push_back(path_info);
                continue;
            }
            
            // Expand to ALL outgoing edges
            for (int i = 0; i < outdegree; ++i) {
                uint64_t next_node = outgoings[i];
                
                if (!sdbg.IsValidEdge(next_node)) {
                    continue;
                }
                
                // Skip if we've already visited this node (avoid loops)
                bool already_visited = false;
                for (const auto& visited : state.node_path) {
                    if (visited == next_node) {
                        already_visited = true;
                        break;
                    }
                }
                if (already_visited) {
                    continue;
                }
                
                BeamState new_state;
                new_state.node_path = state.node_path;
                new_state.node_path.push_back(next_node);
                new_state.current_node = next_node;
                new_state.depth = depth + 1;
                
                // Get sequence and extend
                std::string node_seq = GetNodeSequence(next_node);
                
                // For de Bruijn graph, only the last character is new
                std::string extended_seq = state.accumulated_sequence;
                if (!node_seq.empty()) {
                    extended_seq += node_seq.back();
                }
                new_state.accumulated_sequence = extended_seq;
                
                // Process triplets into amino acids
                new_state.amino_acids = state.amino_acids;
                new_state.hmm_position = state.hmm_position;
                size_t seq_len = extended_seq.length();
                size_t num_complete_codons = seq_len / 3;
                
                // Convert complete triplets to amino acids
                for (size_t j = state.amino_acids.size(); j < num_complete_codons; ++j) {
                    std::string codon = extended_seq.substr(j * 3, 3);
                    char aa = CodonToAminoAcid(codon);
                    new_state.amino_acids.push_back(std::string(1, aa));
                }
                
                // Score entire amino acid sequence with Viterbi alignment
                double viterbi_score = 0.0;
                if (profile_ != nullptr && !new_state.amino_acids.empty()) {
                    // Run Viterbi on complete sequence to get optimal alignment score
                    auto [bit_score, alignment_path, hmm_end_pos] = profile_->ViterbiAlign(new_state.amino_acids);
                    viterbi_score = bit_score;
                    
                    // Update HMM position to the ending position from Viterbi
                    new_state.hmm_position = hmm_end_pos;
                    
                    // Early stopping criteria: check insertion/deletion ratios
                    // Count insertions and deletions from alignment path
                    int insertion_count = 0;
                    int deletion_count = 0;
                    for (char state_char : alignment_path) {
                        if (state_char == 'I') insertion_count++;
                        if (state_char == 'D') deletion_count++;
                    }
                    
                    // Get HMM length for threshold calculation
                    int hmm_length = profile_->GetLength();
                    
                    // Skip this path if insertions > 25% or deletions > 15% of HMM length
                    double insertion_ratio = static_cast<double>(insertion_count) / hmm_length;
                    double deletion_ratio = static_cast<double>(deletion_count) / hmm_length;
                    
                    if (insertion_ratio > 0.25 || deletion_ratio > 0.15) {
                        // Skip this path - too many indels
                        continue;
                    }
                } else {
                    // Fallback to simple node degree scoring
                    viterbi_score = static_cast<double>(sdbg.EdgeOutdegree(next_node));
                }
                
                // Use Viterbi score directly as total score (not incremental)
                new_state.scores = state.scores;
                new_state.scores.push_back(viterbi_score);
                new_state.total_score = viterbi_score;  // Total score IS the Viterbi score
                
                next_layer.push_back(new_state);
            }
        }
        
        if (next_layer.empty()) {
            break; // No more paths to explore
        }
        
        // Sort by score and keep top beam_width
        std::sort(next_layer.begin(), next_layer.end(), 
                  [](const BeamState& a, const BeamState& b) {
                      return a.total_score > b.total_score;
                  });
        
        if (next_layer.size() > static_cast<size_t>(beam_width)) {
            next_layer.resize(beam_width);
        }
        
        current_layer = std::move(next_layer);
    }
    
    // Add remaining paths
    for (const auto& state : current_layer) {
        AminoAcidPathInfo path_info;
        path_info.amino_acids = state.amino_acids;
        path_info.node_path = state.node_path;
        path_info.scores = state.scores;
        path_info.total_score = state.total_score;
        path_info.dna_sequence = state.accumulated_sequence;
        path_info.hmm_position = state.hmm_position;
        result_paths.push_back(path_info);
    }
}

std::vector<AminoAcidPathInfo> CasGeneDetector::BeamSearchAminoAcids(
    uint64_t start_node,
    int beam_width,
    int max_depth) {
    
    std::vector<AminoAcidPathInfo> result_paths;
    
    if (!sdbg.IsValidEdge(start_node)) {
        return result_paths;
    }
    
    BeamTraverse(start_node, beam_width, max_depth, result_paths);
    
    return result_paths;
}
