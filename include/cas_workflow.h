#ifndef CAS_WORKFLOW_H
#define CAS_WORKFLOW_H

#include <vector>
#include <string>
#include <map>
#include <cstdint>
#include "sdbg/sdbg.h"
#include "profile.h"
#include "buckets.h"

/**
 * @brief Workflow parameters for CAS gene detection
 */
struct CasWorkflowParams {
    int FIRST_GENE_MIN_DIST = 50;      // Min distance from repeat to first gene
    int FIRST_GENE_MAX_DIST = 1000;    // Max distance from repeat to first gene
    int MAX_START_CANDIDATES = 5000;   // Max start codon candidates to evaluate
    int MAX_LOCUS_BP = 41591;          // Max total locus length
    int OVERLAP_ALLOWANCE = 15;        // Allow overlap between genes
    int INTERGENIC_MAX = 100;          // Max gap between genes
    int BEAM_WIDTH = 10;               // Beam search width (reduced to save RAM)
    double MIN_NORMALIZED_SCORE = 0.5; // Minimum normalized score
};

/**
 * @brief Detected CAS gene information
 */
struct DetectedCasGene {
    std::string profile_name;          // HMM profile that matched
    uint64_t start_node;               // Starting node in graph
    uint64_t end_node;                 // Ending node in graph
    int distance_from_repeat;          // Distance from repeat to start
    int gene_length;                   // Gene length in bp
    double normalized_score;           // bit_score / hmm_positions
    double bit_score;                  // Raw Viterbi bit score
    std::string amino_acids;           // Amino acid sequence
    std::vector<uint64_t> node_path;   // Path through graph
    bool is_complete;                  // True if has stop codon or HMM fully aligned
    
    DetectedCasGene() 
        : start_node(0), end_node(0), distance_from_repeat(0), 
          gene_length(0), normalized_score(-1e9), bit_score(0.0),
          is_complete(false) {}
};

/**
 * @brief Start codon candidate with position info
 */
struct StartCodonCandidate {
    uint64_t node;
    int distance;
    int offset;  // Position of start codon within k-mer
};

/**
 * @brief CAS gene detection workflow
 */
class CasWorkflow {
public:
    CasWorkflow(SDBG& sdbg, const std::string& profiles_dir);
    
    /**
     * @brief Detect CAS genes starting from a repeat node
     * @param repeat_node The CRISPR repeat node to start from
     * @return Vector of detected CAS genes
     */
    std::vector<DetectedCasGene> DetectCasGenes(uint64_t repeat_node);
    
    /**
     * @brief Set workflow parameters
     */
    void SetParams(const CasWorkflowParams& params) { params_ = params; }
    
    /**
     * @brief Get current parameters
     */
    const CasWorkflowParams& GetParams() const { return params_; }
    
private:
    SDBG& sdbg_;
    std::string profiles_dir_;
    CasWorkflowParams params_;
    
    /**
     * @brief Get nucleotide sequence for a node
     */
    std::string GetNodeSequence(uint64_t node_id);
    
    /**
     * @brief Check if node contains a start codon anywhere
     */
    bool HasStartCodon(uint64_t node_id);
    
    /**
     * @brief Get position of first start codon in node (-1 if none)
     */
    int GetStartCodonOffset(uint64_t node_id);
    
    /**
     * @brief Check if sequence contains an in-frame stop codon
     */
    bool HasStopCodon(const std::string& sequence);
    
    /**
     * @brief Find start codon candidates via BFS
     */
    std::vector<StartCodonCandidate> FindStartCodonCandidates(
        uint64_t repeat_node,
        int min_dist,
        int max_dist,
        int max_candidates
    );
    
    /**
     * @brief Score a start node with a specific HMM profile
     */
    DetectedCasGene ScoreStartNodeWithProfile(
        uint64_t start_node,
        int distance_from_repeat,
        int start_codon_offset,
        const HMMProfiles::ProfileSize& profile_size
    );
    
    /**
     * @brief Score a start node with all HMM profiles
     */
    DetectedCasGene ScoreStartNodeWithAllProfiles(
        uint64_t start_node,
        int distance_from_repeat,
        int start_codon_offset
    );
};

#endif // CAS_WORKFLOW_H
