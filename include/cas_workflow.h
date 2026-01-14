/**
 * @file cas_workflow.h
 * @brief CAS gene detection workflow implementation
 * 
 * This module implements the iterative CAS gene detection workflow:
 * 1. Find ORF at distance D from repeat node
 * 2. Select HMM profiles based on ORF length
 * 3. Score all candidate profiles using Viterbi beam search
 * 4. Select best scoring profile
 * 5. Find next ORF and repeat until total length reaches threshold
 */

#ifndef INCLUDE_CAS_WORKFLOW_H_
#define INCLUDE_CAS_WORKFLOW_H_

#include <vector>
#include <string>
#include <cstdint>
#include <optional>
#include "sdbg/sdbg.h"
#include "orf_finder.h"
#include "buckets.h"
#include "profile.h"
#include "cas_gene_detector.h"

/**
 * @brief Result of HMM profile scoring
 */
struct ProfileScoringResult {
    std::string profile_name;
    double score;
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    int hmm_end_position;
    
    ProfileScoringResult() : profile_name(""), score(0.0), hmm_end_position(0) {}
};

/**
 * @brief Detected CAS gene with complete information
 */
struct DetectedCasGene {
    std::string gene_name;
    uint64_t start_node;
    uint64_t end_node;
    int distance_from_repeat;
    int gene_length;
    double score;
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
    
    DetectedCasGene() : gene_name(""), start_node(SDBG::kNullID), 
                        end_node(SDBG::kNullID), distance_from_repeat(0), 
                        gene_length(0), score(0.0) {}
};

/**
 * @brief Complete CAS operon detection result
 */
struct CasOperonResult {
    std::vector<DetectedCasGene> genes;
    int total_length;  // Total distance from repeat to last gene end
    uint64_t repeat_node;
    
    CasOperonResult() : total_length(0), repeat_node(SDBG::kNullID) {}
};

/**
 * @brief CAS gene detection workflow
 */
class CasWorkflow {
private:
    SDBG& sdbg;
    ORFFinder orf_finder;
    std::string profile_directory;
    
    // Workflow parameters
    int initial_search_distance_min;
    int initial_search_distance_max;
    int subsequent_search_distance;  // Distance for finding next ORFs (default 50bp)
                                     // Note: Searches from last ~10 nodes of previous gene to handle 1-4bp overlaps
    int max_total_length;  // Maximum total length to search (default 41591bp)
    int beam_width;  // Beam width for beam search
    int max_search_depth;  // Maximum depth for beam search
    
    /**
     * @brief Load HMM profile from file
     * @param profile_filename Name of the profile file
     * @return Loaded profile or nullptr if failed
     */
    Profile* LoadProfile(const std::string& profile_filename);
    
    /**
     * @brief Score an ORF with a specific HMM profile
     * @param orf The ORF to score
     * @param profile The HMM profile to use
     * @return Scoring result with score and alignment information
     */
    ProfileScoringResult ScoreORFWithProfile(
        const ORFInfo& orf,
        Profile* profile
    );
    
    /**
     * @brief Score an ORF with all candidate profiles and select the best
     * @param orf The ORF to score
     * @param profile_filenames List of candidate profile filenames
     * @return Best scoring result
     */
    ProfileScoringResult ScoreORFWithBestProfile(
        const ORFInfo& orf,
        const std::vector<std::string>& profile_filenames
    );

public:
    /**
     * @brief Constructor
     * @param sdbg Reference to the sequence De Bruijn graph
     * @param profile_dir Directory containing HMM profile files
     */
    CasWorkflow(SDBG& sdbg, const std::string& profile_dir);
    
    /**
     * @brief Set workflow parameters
     * @param init_min Minimum initial search distance (default 50bp)
     * @param init_max Maximum initial search distance (default 5000bp) 
     * @param subsequent_dist Max distance for subsequent ORF searches (default 100bp)
     * @param max_len Maximum total operon length (default 41591)
     * @param beam Beam width for search (default 10)
     * @param depth Maximum search depth (default 500)
     */
    void SetParameters(
        int init_min = 50,
        int init_max = 5000,
        int subsequent_dist = 100,
        int max_len = 41591,
        int beam = 10,
        int depth = 500
    );
    
    /**
     * @brief Run the complete CAS detection workflow
     * 
     * Iteratively finds and scores ORFs:
     * 1. Find first ORF at distance D from repeat
     * 2. Filter HMM profiles by ORF length
     * 3. Score with all candidate profiles
     * 4. Select best profile
     * 5. Find next ORF and repeat
     * 
     * Stops when total length reaches max_total_length
     * 
     * @param repeat_node The repeat node to start from
     * @return Complete operon detection result
     */
    CasOperonResult DetectCasOperon(uint64_t repeat_node);
};

#endif  // INCLUDE_CAS_WORKFLOW_H_
