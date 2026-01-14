/**
 * @file cas_workflow.cpp
 * @brief Implementation of CAS gene detection workflow
 */

#include "cas_workflow.h"
#include <algorithm>
#include <limits>
#include <iostream>
#include <string>
#include <unordered_set>

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
        
        // Step 1: Find ALL possible ORFs in the search range
        std::cout << "Step 1: Finding ALL ORFs at distance [" 
                  << search_min << ", " << search_max << "] from current position..." 
                  << std::endl;
        
        // Get minimum ORF length from HMM profiles (shortest profile is 108bp)
        int min_orf_len = HMMProfiles::HMMPicker::getMinimumORFLength();
        
        std::vector<ORFInfo> all_orfs;
        
        if (is_first_orf) {
            // First gene: search from repeat node at 50-5000bp
            all_orfs = orf_finder.FindAllORFs(
                current_search_node,
                search_min,
                search_max,
                min_orf_len,
                5000  // Max 5000 candidates
            );
        } else {
            // Subsequent genes: search from END of previous gene
            // Simply search forward 0-100bp to catch overlaps and gaps
            all_orfs = orf_finder.FindAllORFs(
                current_search_node,  // This is already set to prev gene's end_node
                search_min,  // 0
                search_max,  // 100
                min_orf_len,
                5000  // Max 5000 candidates
            );
        }
        
        if (all_orfs.empty()) {
            std::cout << "No ORFs found. Stopping." << std::endl;
            break;
        }
        
        std::cout << "Found " << all_orfs.size() << " candidate ORF(s)" << std::endl;
        
        // Step 2-4: Score ALL ORFs with their matching profiles and select the best
        std::cout << "Step 2-4: Scoring all candidate ORFs..." << std::endl;
        
        ProfileScoringResult best_overall;
        best_overall.score = -std::numeric_limits<double>::infinity();
        ORFInfo best_orf;
        int best_absolute_distance = 0;
        
        for (const auto& orf : all_orfs) {
            // orf.distance_from_repeat is distance from current_search_node
            // current_distance_from_repeat is distance from repeat to current_search_node
            // So absolute distance from repeat = current_distance + orf distance
            int absolute_distance_to_orf_start = current_distance_from_repeat + orf.distance_from_repeat;
            
            std::cout << "  Evaluating ORF: length=" << orf.orf_length 
                      << "bp, start at " << absolute_distance_to_orf_start << "bp from repeat" << std::endl;
            
            // Select profiles based on ORF length
            std::vector<HMMProfiles::ProfileSize> matching_profiles = 
                HMMProfiles::HMMPicker::filterByLength(orf.orf_length);
            
            if (matching_profiles.empty()) {
                std::cout << "    No matching profiles, skipping." << std::endl;
                continue;
            }
            
            // Extract profile filenames
            std::vector<std::string> profile_filenames;
            for (const auto& ps : matching_profiles) {
                profile_filenames.push_back(ps.filename);
            }
            
            // Score with all matching profiles
            ProfileScoringResult scoring_result = ScoreORFWithBestProfile(
                orf,
                profile_filenames
            );
            
            // Update best if this is better
            if (scoring_result.score > best_overall.score) {
                best_overall = scoring_result;
                best_orf = orf;
                best_absolute_distance = absolute_distance_to_orf_start;
                std::cout << "    *** New best score: " << scoring_result.score 
                          << " (profile: " << scoring_result.profile_name << ")" << std::endl;
            }
        }
        
        if (best_overall.score <= -std::numeric_limits<double>::infinity()) {
            std::cout << "All ORF candidates failed to score. Stopping." << std::endl;
            break;
        }
        
        std::cout << "Best ORF selected with score: " << best_overall.score << std::endl;
        
        // Store detected gene
        DetectedCasGene gene;
        gene.gene_name = best_overall.profile_name;
        gene.start_node = best_orf.start_node;
        gene.end_node = best_orf.end_node;
        gene.distance_from_repeat = best_absolute_distance;
        gene.gene_length = best_orf.orf_length;
        gene.score = best_overall.score;
        gene.amino_acids = best_overall.amino_acids;
        gene.node_path = best_overall.node_path;
        
        operon_result.genes.push_back(gene);
        
        std::cout << "Detected CAS gene #" << operon_result.genes.size() 
                  << ": " << gene.gene_name << std::endl;
        
        // Update total length (distance to end of this ORF)
        operon_result.total_length = best_absolute_distance + best_orf.orf_length;
        
        std::cout << "Total operon length so far: " 
                  << operon_result.total_length << " bp" << std::endl;
        
        // Check if we've reached max length
        if (operon_result.total_length >= max_total_length) {
            std::cout << "Reached maximum total length (" 
                      << max_total_length << " bp). Stopping." << std::endl;
            break;
        }
        
        // Step 5: Prepare for next iteration
        if (is_first_orf) {
            is_first_orf = false;
        }
        
        // Update search parameters for subsequent genes
        search_min = 0;
        search_max = subsequent_search_distance;  // 100bp max
        
        // Update position: next search starts from END of current gene
        // current_distance_from_repeat should be distance to END of current gene
        current_search_node = best_orf.end_node;
        current_distance_from_repeat = best_absolute_distance + best_orf.orf_length;
    }
    
    std::cout << "\n=== CAS Operon Detection Complete ===" << std::endl;
    std::cout << "Total genes detected: " << operon_result.genes.size() << std::endl;
    std::cout << "Total operon length: " << operon_result.total_length << " bp" << std::endl;
    
    return operon_result;
}
