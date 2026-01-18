#include <iostream>
#include <string>
#include <queue>
#include <set>
#include <vector>
#include <cstdint>
#include <climits>
#include <iomanip>
#include "sdbg/sdbg.h"
#include "profile.h"
#include "buckets.h"
#include "cas_gene_detector.h"

std::string GetNodeSeq(SDBG& sdbg, uint64_t node_id) {
    const uint32_t k = sdbg.k();
    std::vector<uint8_t> seq(k);
    sdbg.GetLabel(node_id, seq.data());
    std::string result;
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

int GetStartCodonOffset(const std::string& seq) {
    for (size_t i = 0; i + 2 < seq.length(); ++i) {
        std::string codon = seq.substr(i, 3);
        if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

struct StartCodon {
    uint64_t node;
    int distance;
    int offset;
};

std::vector<StartCodon> FindAllStartCodons(SDBG& sdbg, uint64_t repeat_node, int min_dist, int max_dist) {
    std::vector<StartCodon> results;
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    q.push({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!q.empty()) {
        auto [node, dist] = q.front();
        q.pop();
        
        if (dist <= max_dist) {
            std::string seq = GetNodeSeq(sdbg, node);
            // Check all positions in the k-mer for start codons
            for (size_t i = 0; i + 2 < seq.length(); ++i) {
                std::string codon = seq.substr(i, 3);
                if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
                    int bp_dist = dist + static_cast<int>(i);
                    if (bp_dist >= min_dist && bp_dist <= max_dist) {
                        results.push_back({node, bp_dist, static_cast<int>(i)});
                    }
                }
            }
        }
        
        if (dist < max_dist) {
            uint64_t out[4];
            int outdeg = sdbg.OutgoingEdges(node, out);
            for (int i = 0; i < outdeg; ++i) {
                if (sdbg.IsValidEdge(out[i]) && !visited.count(out[i])) {
                    q.push({out[i], dist + 1});
                    visited.insert(out[i]);
                }
            }
        }
    }
    return results;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <graph_prefix> <profiles_dir> <repeat_sequence> [beam_width] [max_depth]" << std::endl;
        return 1;
    }
    
    std::string graph_prefix = argv[1];
    std::string profiles_dir = argv[2];
    std::string repeat_seq = argv[3];
    int beam_width = (argc > 4) ? std::stoi(argv[4]) : 20;
    int max_depth = (argc > 5) ? std::stoi(argv[5]) : 3000;
    
    std::cout << "=== CAS Gene Detection Test (Beam Search) ===" << std::endl;
    std::cout << "Beam width: " << beam_width << std::endl;
    std::cout << "Max depth: " << max_depth << " bp" << std::endl;
    std::cout << std::endl;
    
    SDBG sdbg;
    sdbg.LoadFromFile(graph_prefix.c_str());
    std::cout << "Graph loaded: " << sdbg.size() << " nodes, k=" << sdbg.k() << std::endl;
    
    // Encode repeat sequence
    std::vector<uint8_t> repeat_encoded(repeat_seq.length());
    for (size_t i = 0; i < repeat_seq.length(); ++i) {
        switch (repeat_seq[i]) {
            case 'A': case 'a': repeat_encoded[i] = 1; break;
            case 'C': case 'c': repeat_encoded[i] = 2; break;
            case 'G': case 'g': repeat_encoded[i] = 3; break;
            case 'T': case 't': repeat_encoded[i] = 4; break;
            default: repeat_encoded[i] = 0; break;
        }
    }
    uint64_t repeat_node = sdbg.IndexBinarySearch(repeat_encoded.data());
    
    if (repeat_node == UINT64_MAX) {
        std::cerr << "ERROR: Repeat sequence not found in graph" << std::endl;
        return 1;
    }
    std::cout << "Repeat node: " << repeat_node << std::endl;
    
    // Find ALL start codons in range 150-400
    auto all_starts = FindAllStartCodons(sdbg, repeat_node, 150, 400);
    std::cout << "Found " << all_starts.size() << " start codons in range [150, 400] bp" << std::endl;
    
    if (all_starts.empty()) {
        std::cerr << "ERROR: No start codons found in range" << std::endl;
        return 1;
    }
    
    // Show all start codons found
    std::cout << "\nStart codons:" << std::endl;
    for (size_t i = 0; i < all_starts.size(); ++i) {
        std::cout << "  [" << i << "] Node: " << all_starts[i].node 
                  << ", Distance: " << all_starts[i].distance << " bp"
                  << ", Offset: " << all_starts[i].offset << std::endl;
    }
    std::cout << std::endl;
    
    // For each start codon, run all profiles and find best hit
    std::cout << std::left << std::setw(8) << "START"
              << std::setw(8) << "DIST"
              << std::setw(40) << "BEST_PROFILE" 
              << std::right << std::setw(12) << "BIT_SCORE"
              << std::setw(10) << "HMM_POS"
              << std::setw(8) << "AA_LEN"
              << std::setw(10) << "COMPLETE"
              << std::endl;
    std::cout << std::string(96, '-') << std::endl;
    
    // Track global best
    double global_best_score = -1e9;
    int global_best_start_idx = -1;
    std::string global_best_profile;
    AminoAcidPathInfo global_best_path;
    
    for (size_t si = 0; si < all_starts.size(); ++si) {
        const auto& start = all_starts[si];
        
        double best_score = -1e9;
        std::string best_profile_name;
        AminoAcidPathInfo best_path_info;
        
        for (const auto& ps : HMMProfiles::ALL_PROFILES) {
            Profile hmm;
            std::string profile_path = profiles_dir + "/" + ps.filename;
            
            if (!hmm.LoadFromFile(profile_path)) continue;
            
            CasGeneDetector detector(sdbg, &hmm);
            auto paths = detector.BeamSearchAminoAcids(
                start.node, 
                beam_width, 
                max_depth, 
                start.offset
            );
            
            if (!paths.empty() && paths[0].total_score > best_score) {
                best_score = paths[0].total_score;
                best_profile_name = ps.filename;
                best_path_info = paths[0];
            }
        }
        
        if (best_score > -1e9) {
            std::cout << std::left << std::setw(8) << si
                      << std::setw(8) << start.distance
                      << std::setw(40) << best_profile_name
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(12) << best_score
                      << std::setw(10) << best_path_info.hmm_position
                      << std::setw(8) << best_path_info.amino_acids.size()
                      << std::setw(10) << (best_path_info.is_complete ? "YES" : "NO")
                      << std::endl;
            
            if (best_score > global_best_score) {
                global_best_score = best_score;
                global_best_start_idx = static_cast<int>(si);
                global_best_profile = best_profile_name;
                global_best_path = best_path_info;
            }
        }
    }
    
    std::cout << std::endl;
    std::cout << "=== Global Best Result ===" << std::endl;
    
    if (global_best_score > -1e9) {
        std::cout << "Start codon index: " << global_best_start_idx << std::endl;
        std::cout << "Start codon distance: " << all_starts[global_best_start_idx].distance << " bp" << std::endl;
        std::cout << "Best profile: " << global_best_profile << std::endl;
        std::cout << "Bit score: " << std::fixed << std::setprecision(2) << global_best_score << std::endl;
        std::cout << "HMM position reached: " << global_best_path.hmm_position << std::endl;
        std::cout << "Complete alignment: " << (global_best_path.is_complete ? "YES" : "NO") << std::endl;
        std::cout << "DNA length: " << global_best_path.dna_sequence.length() << " bp" << std::endl;
        std::cout << "AA length: " << global_best_path.amino_acids.size() << std::endl;
        
        // Show first 100 amino acids
        std::cout << "\nFirst 100 AA: ";
        for (size_t i = 0; i < std::min((size_t)100, global_best_path.amino_acids.size()); ++i) {
            std::cout << global_best_path.amino_acids[i];
        }
        std::cout << std::endl;
        
        // Show first 300 bp of DNA
        std::cout << "\nFirst 300 bp DNA: ";
        std::cout << global_best_path.dna_sequence.substr(0, std::min((size_t)300, global_best_path.dna_sequence.length())) << std::endl;
        
        // Also show full AA for HMMER validation
        std::cout << "\nFull AA sequence:" << std::endl;
        for (size_t i = 0; i < global_best_path.amino_acids.size(); ++i) {
            std::cout << global_best_path.amino_acids[i];
        }
        std::cout << std::endl;
    } else {
        std::cout << "No valid paths found with any profile" << std::endl;
    }
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    return 0;
}