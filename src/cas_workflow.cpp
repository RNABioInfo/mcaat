/**
 * @file cas_workflow.cpp
 * @brief Implementation of CAS gene detection workflow
 */

#include "cas_workflow.h"
#include <algorithm>
#include <limits>
#include <iostream>
#include <string>

CasWorkflow::CasWorkflow(SDBG& sdbg, const std::string& profile_dir)
    : sdbg(sdbg), 
      orf_finder(sdbg),
      profile_directory(profile_dir),
      initial_search_distance_min(0),
      initial_search_distance_max(1000),
      subsequent_search_distance(50),  // CAS genes typically 0-50bp apart, most <20bp
      max_total_length(41591),
      beam_width(10),
      max_search_depth(500) {
}

void CasWorkflow::SetParameters(
    int init_min,
    int init_max,
    int subsequent_dist,
    int max_len,
    int beam,
    int depth
) {
    initial_search_distance_min = init_min;
    initial_search_distance_max = init_max;
    subsequent_search_distance = subsequent_dist;
    max_total_length = max_len;
    beam_width = beam;
    max_search_depth = depth;
}

Profile* CasWorkflow::LoadProfile(const std::string& profile_filename) {
    std::string full_path = profile_directory + "/" + profile_filename;
    
    Profile* profile = new Profile();
    if (!profile->LoadFromFile(full_path)) {
        delete profile;
        return nullptr;
    }
    
    return profile;
}

ProfileScoringResult CasWorkflow::ScoreORFWithProfile(
    const ORFInfo& orf,
    Profile* profile
) {
    ProfileScoringResult result;
    result.profile_name = profile->GetName();
    
    // Create CasGeneDetector with the profile
    CasGeneDetector detector(sdbg, profile);
    
    // Calculate max depth based on ORF length
    // Each amino acid needs 3 bp, and we traverse node by node (k-mer by k-mer)
    uint32_t k = sdbg.k();
    int estimated_depth = (orf.orf_length / 3) * 3 / static_cast<int>(k) + 100;  // Add buffer
    int search_depth = std::min(estimated_depth, max_search_depth);
    
    // Run beam search from ORF start node
    std::vector<AminoAcidPathInfo> paths = detector.BeamSearchAminoAcids(
        orf.start_node,
        beam_width,
        search_depth
    );
    
    // Find the best scoring path
    if (paths.empty()) {
        result.score = -std::numeric_limits<double>::infinity();
        return result;
    }
    
    // Select path with highest score
    auto best_path = std::max_element(
        paths.begin(), 
        paths.end(),
        [](const AminoAcidPathInfo& a, const AminoAcidPathInfo& b) {
            return a.total_score < b.total_score;
        }
    );
    
    result.score = best_path->total_score;
    result.amino_acids = best_path->amino_acids;
    result.node_path = best_path->node_path;
    result.hmm_end_position = best_path->hmm_position;
    
    return result;
}

ProfileScoringResult CasWorkflow::ScoreORFWithBestProfile(
    const ORFInfo& orf,
    const std::vector<std::string>& profile_filenames
) {
    ProfileScoringResult best_result;
    best_result.score = -std::numeric_limits<double>::infinity();
    
    std::cout << "Scoring ORF (length=" << orf.orf_length 
              << "bp) with " << profile_filenames.size() 
              << " candidate profiles..." << std::endl;
    
    for (const auto& profile_filename : profile_filenames) {
        // Load profile
        Profile* profile = LoadProfile(profile_filename);
        if (!profile) {
            std::cerr << "Failed to load profile: " << profile_filename << std::endl;
            continue;
        }
        
        // Score ORF with this profile
        ProfileScoringResult result = ScoreORFWithProfile(orf, profile);
        result.profile_name = profile_filename;  // Store filename for later use
        
        std::cout << "  Profile " << profile_filename 
                  << " score: " << result.score << std::endl;
        
        // Update best if this is better
        if (result.score > best_result.score) {
            best_result = result;
        }
        
        delete profile;
    }
    
    if (best_result.score > -std::numeric_limits<double>::infinity()) {
        std::cout << "Best profile: " << best_result.profile_name 
                  << " (score=" << best_result.score << ")" << std::endl;
    }
    
    return best_result;
}

CasOperonResult CasWorkflow::DetectCasOperon(uint64_t repeat_node) {
    CasOperonResult operon_result;
    operon_result.repeat_node = repeat_node;
    
    std::cout << "\n=== Starting CAS Operon Detection ===" << std::endl;
    std::cout << "Repeat node: " << repeat_node << std::endl;
    std::cout << "Max total length: " << max_total_length << " bp" << std::endl;
    
    // Current search position (node from which to search for next ORF)
    uint64_t current_search_node = repeat_node;
    int current_distance_from_repeat = 0;
    
    // For the first ORF, use initial search parameters
    int search_min = initial_search_distance_min;
    int search_max = initial_search_distance_max;
    bool is_first_orf = true;
    
    int iteration = 0;
    
    while (current_distance_from_repeat < max_total_length) {
        iteration++;
        std::cout << "\n--- Iteration " << iteration << " ---" << std::endl;
        std::cout << "Current distance from repeat: " 
                  << current_distance_from_repeat << " bp" << std::endl;
        
        // Step 1: Find ORF
        std::cout << "Step 1: Finding ORF at distance [" 
                  << search_min << ", " << search_max << "] from current position..." 
                  << std::endl;
        
        // Get minimum ORF length from HMM profiles (shortest profile is 108bp)
        int min_orf_len = HMMProfiles::HMMPicker::getMinimumORFLength();
        
        std::optional<ORFInfo> orf_opt = orf_finder.FindFirstORF(
            current_search_node,
            search_min,
            search_max,
            min_orf_len  // Skip ORFs shorter than shortest HMM profile
        );
        
        if (!orf_opt.has_value()) {
            std::cout << "No ORF found. Stopping." << std::endl;
            break;
        }
        
        ORFInfo orf = orf_opt.value();
        std::cout << "Found ORF: start_node=" << orf.start_node 
                  << ", length=" << orf.orf_length << " bp"
                  << ", distance_from_search_node=" << orf.distance_from_repeat << " bp" 
                  << std::endl;
        
        // Update absolute distance from original repeat
        int absolute_distance = current_distance_from_repeat + orf.distance_from_repeat;
        
        // Step 2: Select profiles based on ORF length
        std::cout << "Step 2: Selecting HMM profiles for ORF length " 
                  << orf.orf_length << " bp..." << std::endl;
        
        std::vector<HMMProfiles::ProfileSize> matching_profiles = 
            HMMProfiles::HMMPicker::filterByLength(orf.orf_length);
        
        if (matching_profiles.empty()) {
            std::cout << "No matching profiles for this ORF length. Stopping." << std::endl;
            break;
        }
        
        std::cout << "Found " << matching_profiles.size() 
                  << " matching profiles" << std::endl;
        
        // Extract profile filenames
        std::vector<std::string> profile_filenames;
        for (const auto& ps : matching_profiles) {
            profile_filenames.push_back(ps.filename);
        }
        
        // Step 3 & 4: Score with all profiles and select best
        std::cout << "Step 3-4: Scoring ORF with all candidate profiles..." << std::endl;
        
        ProfileScoringResult best_scoring = ScoreORFWithBestProfile(
            orf,
            profile_filenames
        );
        
        if (best_scoring.score <= -std::numeric_limits<double>::infinity()) {
            std::cout << "All profiles failed to score. Stopping." << std::endl;
            break;
        }
        
        // Store detected gene
        DetectedCasGene gene;
        gene.gene_name = best_scoring.profile_name;
        gene.start_node = orf.start_node;
        gene.end_node = orf.end_node;
        gene.distance_from_repeat = absolute_distance;
        gene.gene_length = orf.orf_length;
        gene.score = best_scoring.score;
        gene.amino_acids = best_scoring.amino_acids;
        gene.node_path = best_scoring.node_path;
        
        operon_result.genes.push_back(gene);
        
        std::cout << "Detected CAS gene #" << operon_result.genes.size() 
                  << ": " << gene.gene_name << std::endl;
        
        // Update total length (distance to end of this ORF)
        operon_result.total_length = absolute_distance + orf.orf_length;
        
        std::cout << "Total operon length so far: " 
                  << operon_result.total_length << " bp" << std::endl;
        
        // Check if we've reached max length
        if (operon_result.total_length >= max_total_length) {
            std::cout << "Reached maximum total length (" 
                      << max_total_length << " bp). Stopping." << std::endl;
            break;
        }
        
        // Step 5: Prepare for next iteration
        // To handle 1-4bp overlaps where gene2 start codon appears in tail of gene1:
        // Search from the LAST FEW NODES of the current gene's path, not just end_node
        // This allows finding start codons that overlap with the stop codon region
        
        // Get candidate search nodes from the tail of the gene path
        std::vector<uint64_t> search_candidate_nodes;
        uint32_t k = sdbg.k();
        int overlap_region_nodes = std::min(10, static_cast<int>(gene.node_path.size()));  // Last ~10 nodes to cover overlap region
        
        for (int i = gene.node_path.size() - overlap_region_nodes; i < gene.node_path.size(); ++i) {
            if (i >= 0) {
                search_candidate_nodes.push_back(gene.node_path[i]);
            }
        }
        
        // Find the FIRST ORF from any of these candidate positions
        std::optional<ORFInfo> next_orf_candidate;
        int min_distance_found = std::numeric_limits<int>::max();
        uint64_t best_search_node = orf.end_node;
        int best_node_distance = 0;  // Distance from original repeat to the search node we used
        
        if (is_first_orf) {
            is_first_orf = false;
        }
        
        // Calculate distance from original repeat to each candidate node in the path
        for (size_t idx = 0; idx < search_candidate_nodes.size(); ++idx) {
            uint64_t candidate_node = search_candidate_nodes[idx];
            
            // Approximate distance: absolute_distance + remaining length from this node to end
            int nodes_from_here_to_end = search_candidate_nodes.size() - idx - 1;
            int approx_distance_to_candidate = absolute_distance + orf.orf_length - (nodes_from_here_to_end * static_cast<int>(k));
            
            // Search for ORF from this candidate with range [0, subsequent_search_distance]
            int min_orf_len = HMMProfiles::HMMPicker::getMinimumORFLength();
            std::optional<ORFInfo> candidate_orf = orf_finder.FindFirstORF(
                candidate_node,
                0,
                subsequent_search_distance,
                min_orf_len  // Skip ORFs shorter than shortest HMM profile
            );
            
            if (candidate_orf.has_value()) {
                // Calculate actual distance from original repeat to this new ORF start
                int total_distance = approx_distance_to_candidate + candidate_orf->distance_from_repeat;
                
                if (total_distance < min_distance_found) {
                    min_distance_found = total_distance;
                    next_orf_candidate = candidate_orf;
                    best_search_node = candidate_node;
                    best_node_distance = approx_distance_to_candidate;
                }
            }
        }
        
        // If we found an ORF from the tail region, use it for next iteration
        if (next_orf_candidate.has_value()) {
            current_search_node = best_search_node;
            current_distance_from_repeat = best_node_distance;
            std::cout << "Next search will start from node in tail region (distance=" 
                      << current_distance_from_repeat << "bp from repeat)" << std::endl;
        } else {
            // No ORF found in overlap region, just continue from end_node
            current_search_node = orf.end_node;
            current_distance_from_repeat = absolute_distance + orf.orf_length;
            std::cout << "Next search will start from end_node (distance=" 
                      << current_distance_from_repeat << "bp from repeat)" << std::endl;
        }
        
        search_min = 0;
        search_max = subsequent_search_distance;
    }
    
    std::cout << "\n=== CAS Operon Detection Complete ===" << std::endl;
    std::cout << "Total genes detected: " << operon_result.genes.size() << std::endl;
    std::cout << "Total operon length: " << operon_result.total_length << " bp" << std::endl;
    
    return operon_result;
}
