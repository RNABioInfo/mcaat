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

std::optional<ORFInfo> ORFFinder::ScanForStopCodon(
    uint64_t start_node,
    int max_orf_length
) const {
    const uint32_t k = this->sdbg.k();
    
    // Use BFS to traverse forward from start_node
    std::queue<std::pair<uint64_t, int>> queue;  // (node, distance)
    std::unordered_set<uint64_t> visited;
    
    queue.push({start_node, 0});
    visited.insert(start_node);
    
    while (!queue.empty()) {
        auto [current_node, distance] = queue.front();
        queue.pop();
        
        // Check if we've exceeded max ORF length
        if (distance > max_orf_length) {
            continue;
        }
        
        // Get sequence of current node
        std::string seq = NodeToSequence(current_node);
        
        // Check for stop codon
        if (ContainsStopCodon(seq)) {
            return ORFInfo(start_node, current_node, distance);
        }
        
        // Get outgoing edges
        int outdegree = this->sdbg.EdgeOutdegree(current_node);
        if (outdegree > 0) {
            std::vector<uint64_t> outgoings(outdegree);
            if (this->sdbg.OutgoingEdges(current_node, outgoings.data()) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    uint64_t neighbor = outgoings[i];
                    if (this->sdbg.IsValidEdge(neighbor) && visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        // Distance increases by 1 (since we move by 1 bp in De Bruijn graph with k-mer overlap)
                        queue.push({neighbor, distance + 1});
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
    int max_distance
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
    for (uint64_t start_candidate : candidate_starts) {
        auto orf_result = ScanForStopCodon(start_candidate, MAX_ORF_LENGTH);
        if (orf_result.has_value()) {
            // Return the first ORF found
            return orf_result;
        }
    }
    
    return std::nullopt;  // No complete ORF found
}
