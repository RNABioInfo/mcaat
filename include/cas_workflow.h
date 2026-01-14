#ifndef CAS_WORKFLOW_H
#define CAS_WORKFLOW_H

#include <string>
#include <vector>
#include <unordered_set>
#include "sdbg/sdbg.h"
#include "orf_finder.h"

struct DetectedCasGene {
    std::string gene_name;
    uint64_t start_node;
    uint64_t end_node;
    int distance_from_repeat;  // bp (approx)
    int gene_length;  // bp
    double score;
    std::vector<std::string> amino_acids;
    std::vector<uint64_t> node_path;
};

struct CasOperonResult {
    uint64_t repeat_node;
    std::vector<DetectedCasGene> genes;
    int total_length = 0; // bp
};

class CasWorkflow {
public:
    CasWorkflow(SDBG& sdbg, const std::string& profiles_dir);

    // Set parameters
    void SetParameters(int init_min = 50,
                       int init_max = 5000,
                       int subsequent_dist = 100,
                       int max_len = 41591,
                       int beam = 10,
                       int depth = 500);

    CasOperonResult DetectCasOperon(uint64_t repeat_node);

private:
    SDBG& sdbg;
    ORFFinder orf_finder;
    std::string profiles_dir;

    int initial_search_distance_min;
    int initial_search_distance_max;
    int subsequent_search_distance;
    int max_total_length;
    int beam_width;
    int search_depth;

    // Helpers
    struct ProfileScoringResult {
        double score;
        std::string profile_name;
        std::vector<std::string> amino_acids;
        std::vector<uint64_t> node_path;
        int hmm_end_position;
    };

    ProfileScoringResult ScoreORFWithBestProfile(const ORFInfo& orf, const std::vector<std::string>& profile_filenames);
    ProfileScoringResult ScoreORFWithProfile(const ORFInfo& orf, Profile* profile);
    Profile* LoadProfile(const std::string& filename);
};

#endif // CAS_WORKFLOW_H
