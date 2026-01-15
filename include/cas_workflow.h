#ifndef CAS_WORKFLOW_H
#define CAS_WORKFLOW_H

#include <vector>
#include <string>
#include <cstdint>
#include <map>
#include "sdbg/sdbg.h"
#include "profile.h"
#include "buckets.h"

// DetectedCasGene structure from requirements
struct DetectedCasGene {
    std::string profile_name;
    uint64_t start_node;
    uint64_t end_node;
    int distance_from_repeat;
    int gene_length;
    double normalized_score;
    double bit_score;
    std::string amino_acids;
    std::vector<uint64_t> node_path;
};

// Parameters from cas_workflow.txt
struct CasWorkflowParams {
    int FIRST_GENE_MIN_DIST = 50;
    int FIRST_GENE_MAX_DIST = 1000;
    int MAX_START_CANDIDATES = 5000;
    int MAX_LOCUS_BP = 41591;
    int OVERLAP_ALLOWANCE = 15;
    int INTERGENIC_MAX = 100;
    int BEAM_WIDTH = 10;
    double MIN_NORMALIZED_SCORE = 0.0;
};

class CasWorkflow {
public:
    CasWorkflow(SDBG& sdbg, const std::string& profiles_dir);
    
    // Main workflow entry point - returns list of detected Cas genes for given repeat node
    std::vector<DetectedCasGene> DetectCasGenes(uint64_t repeat_node);
    
    // Set custom parameters (optional)
    void SetParams(const CasWorkflowParams& params) { params_ = params; }
    
private:
    SDBG& sdbg_;
    std::string profiles_dir_;
    CasWorkflowParams params_;
    
    // Helper to get node sequence
    std::string GetNodeSequence(uint64_t node_id);
    
    // Check if node contains start codon (ATG, GTG, TTG)
    bool HasStartCodon(uint64_t node_id);
    
    // Check if sequence contains stop codon (TAA, TAG, TGA)
    bool HasStopCodon(const std::string& sequence);
    
    // STEP 0: Find first gene start codon candidates via BFS
    std::map<uint64_t, int> FindStartCodonCandidates(
        uint64_t repeat_node,
        int min_dist,
        int max_dist,
        int max_candidates
    );
    
    // STEP 1: Score a start node with all HMM profiles, return best
    DetectedCasGene ScoreStartNodeWithAllProfiles(
        uint64_t start_node,
        int distance_from_repeat
    );
    
    // Score a start node with one specific profile
    DetectedCasGene ScoreStartNodeWithProfile(
        uint64_t start_node,
        int distance_from_repeat,
        const HMMProfiles::ProfileSize& profile_size
    );
};

#endif // CAS_WORKFLOW_H
