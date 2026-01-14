/**
 * @file orf_finder.cpp
 * @brief Implementation of ORF finder for sequence De Bruijn graphs
 */

#include "orf_finder.h"
#include <queue>
#include <unordered_set>
#include <unordered_map>

ORFFinder::ORFFinder(const SDBG& sdbg) : sdbg(sdbg) {}

std::string ORFFinder::NodeToSequence(uint64_t node_id) const {
    const uint32_t k = this->sdbg.k();
    uint8_t seq[k];
    this->sdbg.GetLabel(node_id, seq);
    
    std::string result;
    result.reserve(k);
    for (uint32_t i = 0; i < k; ++i) {
        // seq values: 1=A, 2=C, 3=G, 4=T
        switch(seq[i]) {
            case 1: result += 'A'; break;
            case 2: result += 'C'; break;
            case 3: result += 'G'; break;
            case 4: result += 'T'; break;
            default: result += 'N'; break;
        }
    }
    return result;
}

bool ORFFinder::ContainsStartCodon(const std::string& sequence) const {
    const size_t len = sequence.size();
    if (len < 3) return false;
    
    for (size_t i = 0; i <= len - 3; ++i) {
        const std::string codon = sequence.substr(i, 3);
        if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
            return true;
        }
    }
    return false;
}

bool ORFFinder::ContainsStopCodon(const std::string& sequence) const {
    const size_t len = sequence.size();
    if (len < 3) return false;
    
    for (size_t i = 0; i <= len - 3; ++i) {
        const std::string codon = sequence.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return true;
        }
    }
    return false;
}

bool ORFFinder::ContainsInFrameStopCodon(const std::string& sequence, int frame_offset) const {
    const size_t len = sequence.size();
    if (len < 3) return false;
    
    for (size_t i = frame_offset; i + 2 < len; i += 3) {
        const std::string codon = sequence.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return true;
        }
    }
    return false;
}

// Helper: Check if a specific codon is a stop codon
static bool IsStopCodon(const std::string& codon) {
    return codon == "TAA" || codon == "TAG" || codon == "TGA";
}

std::optional<ORFInfo> ORFFinder::ScanForStopCodon(
    uint64_t start_node,
    int distance_from_repeat,
    int max_orf_length
) const {
    const uint32_t k = this->sdbg.k();
    const int MIN_ORF_LENGTH = 108;  // Minimum ORF length in bp (shortest HMM profile)
    
    // Get the starting sequence and find start codon position
    std::string start_seq = NodeToSequence(start_node);
    int start_codon_pos = -1;
    for (size_t i = 0; i + 2 < start_seq.size(); ++i) {
        std::string codon = start_seq.substr(i, 3);
        if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
            start_codon_pos = i;
            break;
        }
    }
    
    if (start_codon_pos == -1) {
        return std::nullopt;
    }
    
    // BFS state: node, accumulated_sequence, path
    struct BFSState {
        uint64_t node;
        std::string sequence;  // Full accumulated sequence from start
        std::vector<uint64_t> path;
    };
    
    std::queue<BFSState> queue;
    std::unordered_set<uint64_t> visited;
    
    // Initialize: start from the start codon position
    BFSState initial;
    initial.node = start_node;
    initial.sequence = start_seq.substr(start_codon_pos);  // Sequence from start codon onwards
    initial.path.push_back(start_node);
    
    queue.push(initial);
    visited.insert(start_node);
    
    while (!queue.empty()) {
        BFSState current = queue.front();
        queue.pop();
        
        int seq_len = current.sequence.length();
        
        // Check if we exceeded max length
        if (seq_len > max_orf_length) {
            continue;
        }
        
        // Check for in-frame stop codon in the accumulated sequence
        // Reading frame is always 0 because we started from the start codon
        // Check only complete codons
        int num_complete_codons = seq_len / 3;
        
        if (num_complete_codons >= MIN_ORF_LENGTH / 3) {
            // Check the last complete codon
            for (int codon_idx = MIN_ORF_LENGTH / 3 - 1; codon_idx < num_complete_codons; ++codon_idx) {
                std::string codon = current.sequence.substr(codon_idx * 3, 3);
                if (IsStopCodon(codon)) {
                    // Found stop codon!
                    int orf_length = (codon_idx + 1) * 3;  // Length up to and including stop codon
                    return ORFInfo(start_node, current.node, distance_from_repeat, orf_length, current.path);
                }
            }
        }
        
        // Expand to neighbors
        int outdegree = this->sdbg.EdgeOutdegree(current.node);
        if (outdegree > 0) {
            std::vector<uint64_t> outgoings(outdegree);
            if (this->sdbg.OutgoingEdges(current.node, outgoings.data()) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    uint64_t neighbor = outgoings[i];
                    if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        
                        // Get the new nucleotide (last char of neighbor's k-mer)
                        std::string neighbor_seq = NodeToSequence(neighbor);
                        char new_nucleotide = neighbor_seq.back();
                        
                        BFSState next;
                        next.node = neighbor;
                        next.sequence = current.sequence + new_nucleotide;
                        next.path = current.path;
                        next.path.push_back(neighbor);
                        
                        queue.push(next);
                    }
                }
            }
        }
    }
    
    return std::nullopt;
}

std::optional<ORFInfo> ORFFinder::FindFirstORF(
    uint64_t repeat_node,
    int min_distance,
    int max_distance,
    int min_orf_length
) const {
    if (!this->sdbg.IsValidEdge(repeat_node)) {
        return std::nullopt;
    }
    
    const uint32_t k = this->sdbg.k();
    
    // BFS to find nodes with start codons at target distance
    std::vector<std::pair<uint64_t, int>> current_layer;
    std::vector<std::pair<uint64_t, int>> next_layer;
    std::unordered_set<uint64_t> visited;
    std::unordered_map<uint64_t, int> candidate_distances;
    
    current_layer.push_back({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!current_layer.empty()) {
        next_layer.clear();
        
        for (const auto& [current_node, edge_count] : current_layer) {
            int sequence_distance = k + edge_count;
            
            if (sequence_distance > max_distance) {
                continue;
            }
            
            if (sequence_distance >= min_distance && sequence_distance <= max_distance) {
                std::string seq = NodeToSequence(current_node);
                if (ContainsStartCodon(seq)) {
                    candidate_distances[current_node] = sequence_distance;
                }
            }
            
            if (sequence_distance < max_distance) {
                int outdegree = this->sdbg.EdgeOutdegree(current_node);
                if (outdegree > 0) {
                    std::vector<uint64_t> outgoings(outdegree);
                    if (this->sdbg.OutgoingEdges(current_node, outgoings.data()) != -1) {
                        for (int i = 0; i < outdegree; ++i) {
                            uint64_t neighbor = outgoings[i];
                            if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                                visited.insert(neighbor);
                                next_layer.push_back({neighbor, edge_count + 1});
                            }
                        }
                    }
                }
            }
        }
        
        current_layer = std::move(next_layer);
    }
    
    // Scan each candidate for complete ORF
    for (const auto& [start_candidate, dist_from_repeat] : candidate_distances) {
        auto orf_result = ScanForStopCodon(start_candidate, dist_from_repeat, MAX_ORF_LENGTH);
        if (orf_result.has_value()) {
            if (min_orf_length > 0 && orf_result->orf_length < min_orf_length) {
                continue;
            }
            return orf_result;
        }
    }
    
    return std::nullopt;
}

std::vector<ORFInfo> ORFFinder::FindAllORFs(
    uint64_t repeat_node,
    int min_distance,
    int max_distance,
    int min_orf_length,
    int max_traverse
) const {
    std::vector<ORFInfo> all_orfs;
    
    uint32_t k = this->sdbg.k();
    
    // BFS to find all nodes with start codons at target distance
    std::vector<std::pair<uint64_t, int>> current_layer;
    std::vector<std::pair<uint64_t, int>> next_layer;
    std::unordered_set<uint64_t> visited;
    std::unordered_map<uint64_t, int> candidate_distances;
    
    current_layer.push_back({repeat_node, 0});
    visited.insert(repeat_node);
    
    if (min_distance < static_cast<int>(k)) {
        std::string start_seq = NodeToSequence(repeat_node);
        if (ContainsStartCodon(start_seq)) {
            candidate_distances[repeat_node] = k;
        }
    }
    
    while (!current_layer.empty()) {
        next_layer.clear();
        
        for (const auto& [current_node, edge_count] : current_layer) {
            int sequence_distance = k + edge_count;
            
            if (sequence_distance > max_distance) {
                continue;
            }
            
            if (sequence_distance >= min_distance && sequence_distance <= max_distance) {
                std::string seq = NodeToSequence(current_node);
                if (ContainsStartCodon(seq)) {
                    candidate_distances[current_node] = sequence_distance;
                    
                    if (static_cast<int>(candidate_distances.size()) >= max_traverse) {
                        goto done_searching;
                    }
                }
            }
            
            if (sequence_distance < max_distance) {
                int outdegree = this->sdbg.EdgeOutdegree(current_node);
                if (outdegree > 0) {
                    std::vector<uint64_t> outgoings(outdegree);
                    if (this->sdbg.OutgoingEdges(current_node, outgoings.data()) != -1) {
                        for (int i = 0; i < outdegree; ++i) {
                            uint64_t neighbor = outgoings[i];
                            if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                                visited.insert(neighbor);
                                next_layer.push_back({neighbor, edge_count + 1});
                            }
                        }
                    }
                }
            }
        }
        
        current_layer = std::move(next_layer);
    }
    
done_searching:
    
    // IMPORTANT: Each ScanForStopCodon needs its own visited set
    // Otherwise ORFs sharing nodes will fail
    for (const auto& [start_candidate, dist_from_repeat] : candidate_distances) {
        auto orf_result = ScanForStopCodon(start_candidate, dist_from_repeat, MAX_ORF_LENGTH);
        if (orf_result.has_value()) {
            if (min_orf_length > 0 && orf_result->orf_length < min_orf_length) {
                continue;
            }
            all_orfs.push_back(orf_result.value());
        }
    }
    
    std::sort(all_orfs.begin(), all_orfs.end(), 
        [](const ORFInfo& a, const ORFInfo& b) {
            return a.distance_from_repeat < b.distance_from_repeat;
        });
    
    return all_orfs;
}
