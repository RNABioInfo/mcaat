#include "cas_workflow.h"
#include "cas_gene_detector.h"
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include <numeric>
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

bool CasWorkflow::HasStartCodon(uint64_t node_id) {
    std::string seq = GetNodeSequence(node_id);
    if (seq.length() < 3) return false;
    
    std::string first_three = seq.substr(0, 3);
    return (first_three == "ATG" || first_three == "GTG" || first_three == "TTG");
}

bool CasWorkflow::HasStopCodon(const std::string& sequence) {
    if (sequence.length() < 3) return false;
    
    for (size_t i = 0; i + 2 < sequence.length(); i += 3) {
        std::string codon = sequence.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return true;
        }
    }
    return false;
}

std::map<uint64_t, int> CasWorkflow::FindStartCodonCandidates(
    uint64_t repeat_node,
    int min_dist,
    int max_dist,
    int max_candidates) {
    
    std::map<uint64_t, int> candidates;  // node -> distance
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    
    q.push({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!q.empty() && static_cast<int>(candidates.size()) < max_candidates) {
        auto [node, dist] = q.front();
        q.pop();
        
        // Check if in valid range and has start codon
        if (dist >= min_dist && dist <= max_dist) {
            if (HasStartCodon(node)) {
                candidates[node] = dist;
            }
        }
        
        // Continue BFS if not exceeded max distance
        if (dist < max_dist) {
            uint64_t outgoings[4];
            int outdegree = sdbg_.OutgoingEdges(node, outgoings);
            
            for (int i = 0; i < outdegree; ++i) {
                uint64_t next_node = outgoings[i];
                if (sdbg_.IsValidEdge(next_node) && visited.find(next_node) == visited.end()) {
                    q.push({next_node, dist + 1});
                    visited.insert(next_node);
                }
            }
        }
    }
    
    return candidates;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithProfile(
    uint64_t start_node,
    int distance_from_repeat,
    const HMMProfiles::ProfileSize& profile_size) {
    
    DetectedCasGene result;
    result.start_node = start_node;
    result.distance_from_repeat = distance_from_repeat;
    result.profile_name = profile_size.filename;
    result.normalized_score = -1e9;
    result.bit_score = 0.0;
    
    // Load HMM profile
    Profile hmm_profile;
    std::string profile_path = profiles_dir_ + "/" + profile_size.filename;
    if (!hmm_profile.LoadFromFile(profile_path)) {
        std::cerr << "Warning: Failed to load profile " << profile_path << std::endl;
        return result;
    }
    
    // Create detector with this profile
    CasGeneDetector detector(sdbg_, &hmm_profile);
    
    // Run BeamSearchAminoAcids with max_depth from profile
    auto paths = detector.BeamSearchAminoAcids(
        start_node,
        params_.BEAM_WIDTH,
        profile_size.max_bp
    );
    
    // Find best path from all returned paths
    for (const auto& path : paths) {
        if (path.amino_acids.empty()) continue;
        
        // Get bit_score and hmm_end_pos from the path
        double bit_score = path.total_score;
        int hmm_end_pos = std::max(1, path.hmm_position);
        
        // Calculate normalized score
        double normalized_score = bit_score / static_cast<double>(hmm_end_pos);
        
        // Keep best
        if (normalized_score > result.normalized_score) {
            result.normalized_score = normalized_score;
            result.bit_score = bit_score;
            result.end_node = path.node_path.empty() ? start_node : path.node_path.back();
            result.gene_length = static_cast<int>(path.dna_sequence.size());
            result.node_path = path.node_path;
            
            // Concatenate amino acids
            result.amino_acids.clear();
            for (const auto& aa : path.amino_acids) {
                result.amino_acids += aa;
            }
        }
    }
    
    return result;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithAllProfiles(
    uint64_t start_node,
    int distance_from_repeat) {
    
    DetectedCasGene best_gene;
    best_gene.normalized_score = -1e9;
    
    const int num_profiles = HMMProfiles::ALL_PROFILES.size();
    std::vector<DetectedCasGene> thread_results(num_profiles);
    
    // Parallelize profile scoring - each profile independently
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_profiles; ++i) {
        thread_results[i] = ScoreStartNodeWithProfile(
            start_node,
            distance_from_repeat,
            HMMProfiles::ALL_PROFILES[i]
        );
    }
    
    // Find best from all thread results
    for (const auto& gene : thread_results) {
        if (gene.normalized_score > best_gene.normalized_score) {
            best_gene = gene;
        }
    }
    
    return best_gene;
}

std::vector<DetectedCasGene> CasWorkflow::DetectCasGenes(uint64_t repeat_node) {
    std::vector<DetectedCasGene> genes;
    int cumulative_length = 0;
    
    // ========================================================================
    // STEP 0 & 1 & 2: Find and score first gene
    // ========================================================================
    
    auto start_candidates = FindStartCodonCandidates(
        repeat_node,
        params_.FIRST_GENE_MIN_DIST,
        params_.FIRST_GENE_MAX_DIST,
        params_.MAX_START_CANDIDATES
    );
    
    if (start_candidates.empty()) {
        std::cout << "No start codon candidates found for first gene" << std::endl;
        return genes;
    }
    
    std::cout << "Found " << start_candidates.size() << " start codon candidates for first gene" << std::endl;
    
    // Score each candidate with all profiles (parallelized)
    DetectedCasGene best_first_gene;
    best_first_gene.normalized_score = params_.MIN_NORMALIZED_SCORE - 1.0;
    
    // Convert map to vector for OpenMP indexing
    std::vector<std::pair<uint64_t, int>> candidates_vec(start_candidates.begin(), start_candidates.end());
    std::vector<DetectedCasGene> candidate_results(candidates_vec.size());
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < candidates_vec.size(); ++i) {
        candidate_results[i] = ScoreStartNodeWithAllProfiles(
            candidates_vec[i].first,
            candidates_vec[i].second
        );
    }
    
    // Find best candidate
    for (const auto& gene : candidate_results) {
        if (gene.normalized_score > best_first_gene.normalized_score) {
            best_first_gene = gene;
        }
    }
    
    // Check if first gene passes threshold
    if (best_first_gene.normalized_score < params_.MIN_NORMALIZED_SCORE) {
        std::cout << "No gene found with score >= " << params_.MIN_NORMALIZED_SCORE << std::endl;
        return genes;
    }
    
    genes.push_back(best_first_gene);
    cumulative_length = best_first_gene.distance_from_repeat + best_first_gene.gene_length;
    
    std::cout << "First gene found: " << best_first_gene.profile_name 
              << " (score=" << best_first_gene.normalized_score << ")" << std::endl;
    
    // ========================================================================
    // STEP 3: Find subsequent genes
    // ========================================================================
    
    while (cumulative_length < params_.MAX_LOCUS_BP) {
        const auto& previous_gene = genes.back();
        int previous_gene_end = previous_gene.distance_from_repeat + previous_gene.gene_length;
        
        int search_start = previous_gene_end - params_.OVERLAP_ALLOWANCE;
        int search_end = previous_gene_end + params_.INTERGENIC_MAX;
        
        auto next_candidates = FindStartCodonCandidates(
            repeat_node,
            search_start,
            search_end,
            params_.MAX_START_CANDIDATES
        );
        
        if (next_candidates.empty()) {
            std::cout << "No more start codon candidates found" << std::endl;
            break;
        }
        
        // Score each candidate (parallelized)
        DetectedCasGene best_next_gene;
        best_next_gene.normalized_score = params_.MIN_NORMALIZED_SCORE - 1.0;
        
        // Convert map to vector for OpenMP indexing
        std::vector<std::pair<uint64_t, int>> next_candidates_vec(next_candidates.begin(), next_candidates.end());
        std::vector<DetectedCasGene> next_candidate_results(next_candidates_vec.size());
        
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < next_candidates_vec.size(); ++i) {
            next_candidate_results[i] = ScoreStartNodeWithAllProfiles(
                next_candidates_vec[i].first,
                next_candidates_vec[i].second
            );
        }
        
        // Find best candidate
        for (const auto& gene : next_candidate_results) {
            if (gene.normalized_score > best_next_gene.normalized_score) {
                best_next_gene = gene;
            }
        }
        
        // Check threshold
        if (best_next_gene.normalized_score < params_.MIN_NORMALIZED_SCORE) {
            std::cout << "No more genes with sufficient score" << std::endl;
            break;
        }
        
        genes.push_back(best_next_gene);
        cumulative_length = best_next_gene.distance_from_repeat + best_next_gene.gene_length;
        
        std::cout << "Gene " << genes.size() << " found: " << best_next_gene.profile_name 
                  << " (score=" << best_next_gene.normalized_score << ")" << std::endl;
    }
    
    std::cout << "Total genes detected: " << genes.size() << std::endl;
    return genes;
}
