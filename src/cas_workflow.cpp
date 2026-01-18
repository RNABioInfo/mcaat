#include "cas_workflow.h"
#include "cas_gene_detector.h"
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include <omp.h>

CasWorkflow::CasWorkflow(SDBG& sdbg, const std::string& profiles_dir)
    : sdbg_(sdbg), profiles_dir_(profiles_dir) {}

std::string CasWorkflow::GetNodeSequence(uint64_t node_id) {
    const uint32_t k = sdbg_.k();
    std::vector<uint8_t> seq(k);
    sdbg_.GetLabel(node_id, seq.data());
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

int CasWorkflow::GetStartCodonOffset(uint64_t node_id) {
    std::string seq = GetNodeSequence(node_id);
    if (seq.length() < 3) return -1;
    for (size_t i = 0; i + 2 < seq.length(); ++i) {
        std::string codon = seq.substr(i, 3);
        if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool CasWorkflow::HasStartCodon(uint64_t node_id) {
    return GetStartCodonOffset(node_id) >= 0;
}

std::vector<StartCodonCandidate> CasWorkflow::FindStartCodonCandidates(
    uint64_t repeat_node, int min_dist, int max_dist, int max_candidates) {
    
    std::vector<StartCodonCandidate> candidates;
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    
    q.push({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!q.empty() && static_cast<int>(candidates.size()) < max_candidates) {
        auto [node, dist] = q.front();
        q.pop();
        
        if (dist <= max_dist) {
            int offset = GetStartCodonOffset(node);
            if (offset >= 0) {
                int bp_dist = dist + offset;
                if (bp_dist >= min_dist && bp_dist <= max_dist) {
                    candidates.push_back({node, bp_dist, offset});
                }
            }
        }
        
        if (dist < max_dist) {
            uint64_t out[4];
            int outdeg = sdbg_.OutgoingEdges(node, out);
            for (int i = 0; i < outdeg; ++i) {
                if (sdbg_.IsValidEdge(out[i]) && !visited.count(out[i])) {
                    q.push({out[i], dist + 1});
                    visited.insert(out[i]);
                }
            }
        }
    }
    return candidates;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithProfile(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    const HMMProfiles::ProfileSize& profile_size) {
    
    DetectedCasGene result;
    result.start_node = start_node;
    result.distance_from_repeat = distance_from_repeat;
    result.profile_name = profile_size.filename;
    result.bit_score = -1e9;
    result.normalized_score = -1e9;
    result.is_complete = false;
    
    Profile hmm;
    if (!hmm.LoadFromFile(profiles_dir_ + "/" + profile_size.filename)) {
        return result;
    }
    
    CasGeneDetector detector(sdbg_, &hmm);
    auto paths = detector.BeamSearchAminoAcids(
        start_node, params_.BEAM_WIDTH, profile_size.max_bp, start_codon_offset);
    
    for (const auto& path : paths) {
        if (path.amino_acids.empty()) continue;
        
        int len_bp = static_cast<int>(path.dna_sequence.size());
        if (len_bp < profile_size.min_bp || len_bp > profile_size.max_bp) continue;
        
        double bit = path.total_score;
        int hmm_pos = std::max(1, path.hmm_position);
        double norm = bit / hmm_pos;
        
        if (bit > result.bit_score) {
            result.bit_score = bit;
            result.normalized_score = norm;
            result.gene_length = len_bp;
            result.end_node = path.node_path.back();
            result.node_path = path.node_path;
            result.is_complete = true;
            result.amino_acids.clear();
            for (const auto& aa : path.amino_acids) result.amino_acids += aa;
        }
    }
    return result;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithAllProfiles(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset) {
    
    DetectedCasGene best;
    best.bit_score = -1e9;
    best.is_complete = false;
    
    const auto& profiles = HMMProfiles::ALL_PROFILES;
    std::vector<DetectedCasGene> results(profiles.size());
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < profiles.size(); ++i) {
        std::cout << "START " << start_node << " : " << profiles[i].filename << std::endl;
        results[i] = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, profiles[i]);
    }
    
    for (const auto& r : results) {
        if (r.is_complete && r.bit_score > best.bit_score) {
            best = r;
        }
    }
    return best;
}

std::vector<DetectedCasGene> CasWorkflow::DetectCasGenes(uint64_t repeat_node) {
    std::vector<DetectedCasGene> genes;
    
    // Step 1: Find start codons
    auto candidates = FindStartCodonCandidates(
        repeat_node,
        params_.FIRST_GENE_MIN_DIST,
        params_.FIRST_GENE_MAX_DIST,
        params_.MAX_START_CANDIDATES
    );
    
    if (candidates.empty()) return genes;
    
    // Step 2: Score all candidates, keep absolute best
    DetectedCasGene best;
    best.bit_score = -1e9;
    best.is_complete = false;
    
    for (const auto& c : candidates) {
        auto result = ScoreStartNodeWithAllProfiles(c.node, c.distance, c.offset);
        if (result.is_complete && result.bit_score > best.bit_score) {
            best = result;
        }
    }
    
    // Return best gene regardless of threshold
    if (best.is_complete) {
        genes.push_back(best);
    }
    
    return genes;
}
