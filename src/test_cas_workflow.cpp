#include <iostream>
#include <string>
#include <queue>
#include <set>
#include <vector>
#include <cstdint>
#include <climits>
#include <unordered_map>
#include "sdbg/sdbg.h"
#include "profile.h"
#include "buckets.h"

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

char CodonToAA(const std::string& codon) {
    static const std::unordered_map<std::string, char> t = {
        {"TTT",'F'},{"TTC",'F'},{"TTA",'L'},{"TTG",'L'},
        {"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},
        {"TAT",'Y'},{"TAC",'Y'},{"TAA",'*'},{"TAG",'*'},
        {"TGT",'C'},{"TGC",'C'},{"TGA",'*'},{"TGG",'W'},
        {"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
        {"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
        {"CAT",'H'},{"CAC",'H'},{"CAA",'Q'},{"CAG",'Q'},
        {"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},
        {"ATT",'I'},{"ATC",'I'},{"ATA",'I'},{"ATG",'M'},
        {"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},
        {"AAT",'N'},{"AAC",'N'},{"AAA",'K'},{"AAG",'K'},
        {"AGT",'S'},{"AGC",'S'},{"AGA",'R'},{"AGG",'R'},
        {"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
        {"GCT",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},
        {"GAT",'D'},{"GAC",'D'},{"GAA",'E'},{"GAG",'E'},
        {"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'}
    };
    if (codon.length() != 3) return 'X';
    auto it = t.find(codon);
    return (it != t.end()) ? it->second : 'X';
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

// Greedy traversal with backtracking via IncomingEdges
std::string GreedyTraverse(SDBG& sdbg, uint64_t start_node, int offset, int max_bp) {
    std::string dna = GetNodeSeq(sdbg, start_node).substr(offset);
    std::set<uint64_t> visited;
    visited.insert(start_node);
    
    uint64_t current = start_node;
    
    for (int i = 0; i < max_bp; ++i) {
        uint64_t out[4];
        int outdeg = sdbg.OutgoingEdges(current, out);
        
        // Find first unvisited outgoing edge
        uint64_t next = UINT64_MAX;
        for (int j = 0; j < outdeg; ++j) {
            if (sdbg.IsValidEdge(out[j]) && !visited.count(out[j])) {
                next = out[j];
                break;
            }
        }
        
        // If stuck, backtrack
        while (next == UINT64_MAX && current != start_node) {
            // Go back via incoming edges
            uint64_t inc[4];
            int indeg = sdbg.IncomingEdges(current, inc);
            
            uint64_t prev = UINT64_MAX;
            for (int j = 0; j < indeg; ++j) {
                if (sdbg.IsValidEdge(inc[j]) && visited.count(inc[j])) {
                    prev = inc[j];
                    break;
                }
            }
            
            if (prev == UINT64_MAX) break;  // Can't backtrack
            
            // Remove last char from dna
            if (!dna.empty()) dna.pop_back();
            
            current = prev;
            
            // Try outgoing edges from this node
            outdeg = sdbg.OutgoingEdges(current, out);
            for (int j = 0; j < outdeg; ++j) {
                if (sdbg.IsValidEdge(out[j]) && !visited.count(out[j])) {
                    next = out[j];
                    break;
                }
            }
        }
        
        if (next == UINT64_MAX) break;  // Truly stuck
        
        visited.insert(next);
        dna += GetNodeSeq(sdbg, next).back();
        current = next;
    }
    
    return dna;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <graph_prefix> <profiles_dir> <repeat_sequence>" << std::endl;
        return 1;
    }
    
    std::string graph_prefix = argv[1];
    std::string profiles_dir = argv[2];
    std::string repeat_seq = argv[3];
    
    SDBG sdbg;
    sdbg.LoadFromFile(graph_prefix.c_str());
    
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
        std::cerr << "Repeat not found" << std::endl;
        return 1;
    }
    
    StartCodon start = FindFirstStartCodon(sdbg, repeat_node, 600, 1000);
    if (start.distance < 0) {
        std::cerr << "No start codon" << std::endl;
        return 1;
    }
    
    std::cout << "Start codon: node=" << start.node << " dist=" << start.distance << " offset=" << start.offset << std::endl;
    
    // Greedy traversal
    std::string dna = GreedyTraverse(sdbg, start.node, start.offset, 3000);
    std::cout << "DNA length: " << dna.length() << " bp" << std::endl;
    
    // Translate - skip stop codons (HMMER style)
    std::string aa;
    for (size_t i = 0; i + 2 < dna.length(); i += 3) {
        char c = CodonToAA(dna.substr(i, 3));
        if (c == '*') continue;  // Skip stop codon
        aa += c;
    }
    std::cout << "AA length: " << aa.length() << std::endl;
    std::cout << "First 100 AA: " << aa.substr(0, 100) << std::endl;
    
    // Score with ViterbiAlign directly
    std::cout << "\nPROFILE\tBIT_SCORE\tHMM_POS" << std::endl;
    
    std::vector<std::string> aa_vec;
    for (char c : aa) aa_vec.push_back(std::string(1, c));
    
    for (const auto& ps : HMMProfiles::ALL_PROFILES) {
        Profile hmm;
        if (!hmm.LoadFromFile(profiles_dir + "/" + ps.filename)) {
            std::cout << ps.filename << "\tLOAD_FAILED" << std::endl;
            continue;
        }
        
        auto [bit_score, path, hmm_pos] = hmm.ViterbiAlign(aa_vec);
        std::cout << ps.filename << "\t" << bit_score << "\t" << hmm_pos << std::endl;
    }
    
    return 0;
}