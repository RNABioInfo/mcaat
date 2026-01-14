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

// Check if sequence contains an IN-FRAME stop codon relative to start
bool ORFFinder::ContainsInFrameStopCodon(const std::string& sequence, int frame_offset) const {
    const size_t len = sequence.size();
    if (len < 3) return false;
    
    // frame_offset tells us which reading frame we're in (0, 1, or 2)
    for (size_t i = frame_offset; i + 2 < len; i += 3) {
        const std::string codon = sequence.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return true;
        }
    }
    return false;
}

std::optional<ORFInfo> ORFFinder::ScanForStopCodon(
    uint64_t start_node,
    int distance_from_repeat,
    int max_orf_length
) const {
    const uint32_t k = this->sdbg.k();
    const int MIN_ORF_LENGTH = 47;  // Minimum 141bp = 47 codons for a real gene
    
    // First, find the reading frame of the start codon in the start node
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
        return std::nullopt;  // No start codon found (shouldn't happen)
    }
    
    // Use BFS to traverse forward from start_node
    // Track (node, edge_count, frame_offset)
    std::queue<std::tuple<uint64_t, int, int>> queue;
    std::unordered_set<uint64_t> visited;
    
    visited.insert(start_node);
    
    // Get outgoing edges from start node
    int outdegree = this->sdbg.EdgeOutdegree(start_node);
    if (outdegree > 0) {
        std::vector<uint64_t> outgoings(outdegree);
        if (this->sdbg.OutgoingEdges(start_node, outgoings.data()) != -1) {
            for (int i = 0; i < outdegree; ++i) {
                uint64_t neighbor = outgoings[i];
                if (this->sdbg.IsValidEdge(neighbor)) {
                    visited.insert(neighbor);
                    // After moving one edge, frame shifts by 1 (since we add 1 bp in de Bruijn graph)
                    // New frame = (start_codon_pos + k + 1) % 3
                    int new_frame = (start_codon_pos + k) % 3;
                    queue.push({neighbor, 1, new_frame});
                }
            }
        }
    }
    
    while (!queue.empty()) {
        auto [current_node, edge_count, frame_offset] = queue.front();
        queue.pop();
        
        // Actual sequence length = k + edge_count
        int sequence_length = k + edge_count;
        
        if (sequence_length > max_orf_length) {
            continue;
        }
        
        // Get sequence of current node
        std::string seq = NodeToSequence(current_node);
        
        // Check for IN-FRAME stop codon, and only if ORF is long enough
        if (sequence_length >= MIN_ORF_LENGTH && ContainsInFrameStopCodon(seq, frame_offset)) {
            return ORFInfo(start_node, current_node, distance_from_repeat, sequence_length);
        }
        
        // Get outgoing edges
        int outdegree_curr = this->sdbg.EdgeOutdegree(current_node);
        if (outdegree_curr > 0) {
            std::vector<uint64_t> outgoings_curr(outdegree_curr);
            if (this->sdbg.OutgoingEdges(current_node, outgoings_curr.data()) != -1) {
                for (int i = 0; i < outdegree_curr; ++i) {
                    uint64_t neighbor = outgoings_curr[i];
                    if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        // Frame shifts by 1 for each edge
                        int new_frame = (frame_offset + 1) % 3;
                        queue.push({neighbor, edge_count + 1, new_frame});
                    }
                }
            }
        }
    }
    
    return std::nullopt;  // No stop codon found
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
    
    // BFS to find nodes at distance D from repeat_node
    // We store layers to keep only the last layer in memory
    std::vector<std::pair<uint64_t, int>> current_layer;  // (node, distance)
    std::vector<std::pair<uint64_t, int>> next_layer;
    std::unordered_set<uint64_t> visited;
    
    current_layer.push_back({repeat_node, 0});
    visited.insert(repeat_node);
    
    // Vector to store candidate start nodes within the distance range
    std::vector<uint64_t> candidate_starts;
    
    while (!current_layer.empty()) {
        next_layer.clear();
        
        for (const auto& [current_node, distance] : current_layer) {
            // If we're beyond max_distance, skip
            if (distance > max_distance) {
                continue;
            }
            
            // Check if current node is within the target distance range
            if (distance >= min_distance && distance <= max_distance) {
                std::string seq = NodeToSequence(current_node);
                if (ContainsStartCodon(seq)) {
                    candidate_starts.push_back(current_node);
                }
            }
            
            // Expand to neighbors (only if not beyond max_distance)
            if (distance < max_distance) {
                int outdegree = this->sdbg.EdgeOutdegree(current_node);
                if (outdegree > 0) {
                    std::vector<uint64_t> outgoings(outdegree);
                    if (this->sdbg.OutgoingEdges(current_node, outgoings.data()) != -1) {
                        for (int i = 0; i < outdegree; ++i) {
                            uint64_t neighbor = outgoings[i];
                            if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                                visited.insert(neighbor);
                                next_layer.push_back({neighbor, distance + 1});
                            }
                        }
                    }
                }
            }
        }
        
        // Move to next layer (only keep last layer in memory)
        current_layer = std::move(next_layer);
    }
    
    // Now scan each candidate start node for the first complete ORF
    // We need to track which distance each candidate was found at
    std::unordered_map<uint64_t, int> candidate_distances;
    
    // Rebuild the BFS to get distances for candidates
    current_layer.clear();
    next_layer.clear();
    visited.clear();
    
    current_layer.push_back({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!current_layer.empty()) {
        next_layer.clear();
        
        for (const auto& [current_node, edge_count] : current_layer) {
            // Actual sequence distance = k + edge_count
            int sequence_distance = k + edge_count;
            
            if (sequence_distance > max_distance) {
                continue;
            }
            
            // Store distance for candidates
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
    
    // Now scan each candidate with its proper distance
    for (const auto& [start_candidate, dist_from_repeat] : candidate_distances) {
        auto orf_result = ScanForStopCodon(start_candidate, dist_from_repeat, MAX_ORF_LENGTH);
        if (orf_result.has_value()) {
            return orf_result;
        }
    }
    
    return std::nullopt;  // No complete ORF found
}
