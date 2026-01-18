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

StartCodon FindFirstStartCodon(SDBG& sdbg, uint64_t repeat_node, int min_dist, int max_dist) {
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    q.push({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!q.empty()) {
        auto [node, dist] = q.front();
        q.pop();
        
        if (dist <= max_dist) {
            std::string seq = GetNodeSeq(sdbg, node);
            int offset = GetStartCodonOffset(seq);
            if (offset >= 0) {
                int bp_dist = dist + offset;
                if (bp_dist >= min_dist && bp_dist <= max_dist) {
                    return {node, bp_dist, offset};
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
    return {0, -1, -1};
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
    
    // Find first start codon in valid range
    StartCodon start = FindFirstStartCodon(sdbg, repeat_node, 250, 1000);
    if (start.distance < 0) {
        std::cerr << "ERROR: No start codon found in range [600, 1000] bp from repeat" << std::endl;
        return 1;
    }
    
    std::cout << "Start codon found:" << std::endl;
    std::cout << "  Node: " << start.node << std::endl;
    std::cout << "  Distance from repeat: " << start.distance << " bp" << std::endl;
    std::cout << "  Offset in k-mer: " << start.offset << std::endl;
    std::cout << std::endl;
    
    // Test with each profile using beam search
    std::cout << std::left << std::setw(40) << "PROFILE" 
              << std::right << std::setw(12) << "BIT_SCORE"
              << std::setw(10) << "HMM_POS"
              << std::setw(10) << "AA_LEN"
              << std::setw(10) << "DNA_LEN"
              << std::setw(10) << "PATHS"
              << std::setw(12) << "COMPLETE"
              << std::endl;
    std::cout << std::string(104, '-') << std::endl;
    
    for (const auto& ps : HMMProfiles::ALL_PROFILES) {
        Profile hmm;
        std::string profile_path = profiles_dir + "/" + ps.filename;
        
        if (!hmm.LoadFromFile(profile_path)) {
            std::cout << std::left << std::setw(40) << ps.filename 
                      << "  LOAD_FAILED" << std::endl;
            continue;
        }
        
        // Create detector with this profile
        CasGeneDetector detector(sdbg, &hmm);
        
        // Run beam search
        auto paths = detector.BeamSearchAminoAcids(
            start.node, 
            beam_width, 
            max_depth, 
            start.offset
        );
        
        if (paths.empty()) {
            std::cout << std::left << std::setw(40) << ps.filename 
                      << "  NO_PATHS" << std::endl;
            continue;
        }
        
        // Best path is first (sorted by score)
        const auto& best = paths[0];
        
        std::cout << std::left << std::setw(40) << ps.filename
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(12) << best.total_score
                  << std::setw(10) << best.hmm_position
                  << std::setw(10) << best.amino_acids.size()
                  << std::setw(10) << best.dna_sequence.length()
                  << std::setw(10) << paths.size()
                  << std::setw(12) << (best.is_complete ? "YES" : "NO")
                  << std::endl;
    }
    
    std::cout << std::endl;
    
    // Detailed output for best overall profile
    std::cout << "=== Detailed Analysis of Best Profile ===" << std::endl;
    
    double best_score = -1e9;
    std::string best_profile;
    AminoAcidPathInfo best_path;
    
    for (const auto& ps : HMMProfiles::ALL_PROFILES) {
        Profile hmm;
        if (!hmm.LoadFromFile(profiles_dir + "/" + ps.filename)) continue;
        
        CasGeneDetector detector(sdbg, &hmm);
        auto paths = detector.BeamSearchAminoAcids(start.node, beam_width, max_depth, start.offset);
        
        if (!paths.empty() && paths[0].total_score > best_score) {
            best_score = paths[0].total_score;
            best_profile = ps.filename;
            best_path = paths[0];
        }
    }
    
    if (best_score > -1e9) {
        std::cout << "Best profile: " << best_profile << std::endl;
        std::cout << "Bit score: " << best_score << std::endl;
        std::cout << "HMM position reached: " << best_path.hmm_position << std::endl;
        std::cout << "Complete alignment: " << (best_path.is_complete ? "YES" : "NO") << std::endl;
        std::cout << "DNA length: " << best_path.dna_sequence.length() << " bp" << std::endl;
        std::cout << "AA length: " << best_path.amino_acids.size() << std::endl;
        std::cout << "Node path length: " << best_path.node_path.size() << " nodes" << std::endl;
        
        // Show first 100 amino acids
        std::cout << "\nFirst 100 AA: ";
        for (size_t i = 0; i < std::min((size_t)100, best_path.amino_acids.size()); ++i) {
            std::cout << best_path.amino_acids[i];
        }
        std::cout << std::endl;
        
        // Show first 300 bp of DNA
        std::cout << "\nFirst 300 bp DNA: ";
        std::cout << best_path.dna_sequence.substr(0, 300) << std::endl;
        
        // Validate: reconstruct DNA from node path
        if (best_path.node_path.size() > 1) {
            std::cout << "\n=== Node Path Validation ===" << std::endl;
            std::cout << "First 10 nodes: ";
            for (size_t i = 0; i < std::min((size_t)10, best_path.node_path.size()); ++i) {
                std::cout << best_path.node_path[i] << " ";
            }
            std::cout << std::endl;
            
            // Reconstruct DNA from path
            std::string reconstructed = GetNodeSeq(sdbg, best_path.node_path[0]).substr(start.offset);
            for (size_t i = 1; i < best_path.node_path.size(); ++i) {
                reconstructed += GetNodeSeq(sdbg, best_path.node_path[i]).back();
            }
            
            bool match = (reconstructed == best_path.dna_sequence);
            std::cout << "DNA reconstruction from path: " << (match ? "MATCH" : "MISMATCH") << std::endl;
            if (!match) {
                std::cout << "  Stored length: " << best_path.dna_sequence.length() << std::endl;
                std::cout << "  Reconstructed length: " << reconstructed.length() << std::endl;
            }
        }
    } else {
        std::cout << "No valid paths found with any profile" << std::endl;
    }
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    return 0;
}