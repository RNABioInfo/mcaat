#include "amino_acidator.h"
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>

AminoAcidator::AminoAcidator(SDBG& sdbg) : sdbg(sdbg) {}

char AminoAcidator::CodonToAminoAcid(const std::string& codon) {
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

std::string AminoAcidator::GetNodeSequence(uint64_t node_id) {
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

void AminoAcidator::BeamTraverse(
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
    };
    
    auto comp = [](const BeamState& a, const BeamState& b) {
        return a.total_score < b.total_score;
    };
    
    std::vector<BeamState> current_layer;
    
    // Initialize with start node
    BeamState initial;
    initial.node_path.push_back(start_node);
    initial.current_node = start_node;
    initial.total_score = static_cast<double>(start_node); 
    initial.scores.push_back(initial.total_score);
    initial.accumulated_sequence = GetNodeSequence(start_node);
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
                size_t seq_len = extended_seq.length();
                size_t num_complete_codons = seq_len / 3;
                
                // Convert complete triplets to amino acids
                for (size_t j = state.amino_acids.size(); j < num_complete_codons; ++j) {
                    std::string codon = extended_seq.substr(j * 3, 3);
                    char aa = CodonToAminoAcid(codon);
                    new_state.amino_acids.push_back(std::string(1, aa));
                }
                
                // Score is node ID for now
                double node_score = static_cast<double>(sdbg.EdgeOutdegree(next_node));
                new_state.scores = state.scores;
                new_state.scores.push_back(node_score);
                new_state.total_score = state.total_score + node_score;
                
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
        result_paths.push_back(path_info);
    }
}

std::vector<AminoAcidPathInfo> AminoAcidator::BeamSearchAminoAcids(
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
