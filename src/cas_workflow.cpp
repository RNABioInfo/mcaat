#include "cas_workflow.h"
#include "profile.h"
#include "cas_gene_detector.h"
#include "buckets.h"
#include <iostream>
#include <limits>
#include <algorithm>
#include <filesystem>

using namespace HMMProfiles;

CasWorkflow::CasWorkflow(SDBG& sdbg_, const std::string& profiles_dir_)
    : sdbg(sdbg_), orf_finder(sdbg_), profiles_dir(profiles_dir_),
      initial_search_distance_min(50), initial_search_distance_max(5000),
      subsequent_search_distance(100), max_total_length(41591),
      beam_width(100), search_depth(500) {}

void CasWorkflow::SetParameters(int init_min, int init_max, int subsequent_dist, int max_len, int beam, int depth) {
    initial_search_distance_min = init_min;
    initial_search_distance_max = init_max;
    subsequent_search_distance = subsequent_dist;
    max_total_length = max_len;
    beam_width = beam;
    search_depth = depth;
}

Profile* CasWorkflow::LoadProfile(const std::string& filename) {
    std::string path = profiles_dir + "/" + filename;
    Profile* p = new Profile();
    if (!p->LoadFromFile(path)) {
        delete p;
        return nullptr;
    }
    return p;
}

CasWorkflow::ProfileScoringResult CasWorkflow::ScoreORFWithProfile(const ORFInfo& orf, const ProfileSize& profile_meta) {
    ProfileScoringResult res;
    res.score = -std::numeric_limits<double>::infinity();
    res.profile_name = profile_meta.filename;

    // Load profile file
    Profile* profile = LoadProfile(profile_meta.filename);
    if (!profile) {
        std::cerr << "    Failed to load profile: " << profile_meta.filename << std::endl;
        return res;
    }

    // Use CasGeneDetector with this profile to perform beam search and score
    int max_depth = profile_meta.max_bp; // use profile's max_bp as dynamic depth (in bp)
    // Convert max_depth (bp) to node depth estimate and cap by workflow search_depth
    int depth_nodes = std::min(std::max(1, max_depth), search_depth);

    CasGeneDetector detector(sdbg, profile);
    auto paths = detector.BeamSearchAminoAcids(orf.start_node, beam_width, depth_nodes);
    if (paths.empty()) {
        delete profile;
        return res;
    }

    // pick best path by total_score
    auto best_it = std::max_element(paths.begin(), paths.end(), [](const AminoAcidPathInfo&a, const AminoAcidPathInfo&b){ return a.total_score < b.total_score; });
    res.score = best_it->total_score;
    res.hmm_node_path = best_it->node_path; // HMM alignment path (for reporting)
    res.node_path = best_it->node_path; // keep for compatibility (but ORF biological path is stored in ORFInfo)
    res.amino_acids = best_it->amino_acids;
    res.hmm_end_position = best_it->hmm_position;

    delete profile;
    return res;
}

CasWorkflow::ProfileScoringResult CasWorkflow::ScoreORFWithBestProfile(const ORFInfo& orf, const std::vector<ProfileSize>& profile_sizes) {
    ProfileScoringResult best;
    best.score = -std::numeric_limits<double>::infinity();
    best.profile_name = "";

    std::cout << "Scoring ORF (length=" << orf.orf_length << "bp) with " << profile_sizes.size() << " candidate profiles..." << std::endl;

    for (const auto& ps : profile_sizes) {
        std::cout << "  Candidate profile: " << ps.filename << " (min_bp=" << ps.min_bp << " max_bp=" << ps.max_bp << ")" << std::endl;
        auto r = ScoreORFWithProfile(orf, ps);
        r.profile_name = ps.filename;
        std::cout << "    Profile " << ps.filename << " score: " << r.score << std::endl;
        if (r.score > best.score) best = r; // r.node_path is HMM alignment path; ORF path is in 'orf' variable when selecting best ORF
    }
    if (best.score > -std::numeric_limits<double>::infinity()) {
        std::cout << "Best profile: " << best.profile_name << " (score=" << best.score << ")" << std::endl;
    }
    return best;
}

CasOperonResult CasWorkflow::DetectCasOperon(uint64_t repeat_node) {
    CasOperonResult res;
    res.repeat_node = repeat_node;

    std::cout << "\n=== Starting CAS Operon Detection ===" << std::endl;
    std::cout << "Repeat node: " << repeat_node << std::endl;
    std::cout << "Max total length: " << max_total_length << " bp" << std::endl;

    uint64_t current_search_node = repeat_node;
    int current_distance_from_repeat = 0;
    int search_min = initial_search_distance_min;
    int search_max = initial_search_distance_max;
    bool is_first_orf = true;

    // track used start+end keys to prevent cycles and duplicates
    std::unordered_set<std::string> used_start_end_keys; // key = start:end

    while (current_distance_from_repeat < max_total_length) {
        std::cout << "\n--- Iteration " << (res.genes.size() + 1) << " ---" << std::endl;
        std::cout << "Current distance from repeat: " << current_distance_from_repeat << " bp" << std::endl;
        std::cout << "Step 1: Finding ALL ORFs at distance [" << search_min << ", " << search_max << "] from current position..." << std::endl;

        int min_orf_len = HMMPicker::getMinimumORFLength();

        std::vector<ORFInfo> all_orfs;

        if (is_first_orf) {
            // find all ORFs in big window (cap candidates to MAX_CANDIDATES)
            const int MAX_CANDIDATES = 5000;
            all_orfs = orf_finder.FindAllORFs(current_search_node, search_min, search_max, min_orf_len, MAX_CANDIDATES);
        } else {
            // search from last ~15 tail nodes
            const auto& prev_gene = res.genes.back();
            std::vector<uint64_t> tail_nodes;
            int num_tail = std::min(15, static_cast<int>(prev_gene.orf_node_path.size()));
            int start_idx = static_cast<int>(prev_gene.orf_node_path.size()) - num_tail;
            for (int i = start_idx; i < static_cast<int>(prev_gene.orf_node_path.size()); ++i) {
                if (i >= 0) tail_nodes.push_back(prev_gene.orf_node_path[i]);
            }
            if (tail_nodes.empty()) tail_nodes.push_back(current_search_node);

            std::unordered_set<std::string> seen_local;
            for (size_t idx = 0; idx < tail_nodes.size(); ++idx) {
                uint64_t tn = tail_nodes[idx];
                int node_index_in_path = start_idx + static_cast<int>(idx);
                // approximate dist from repeat to this tail node
                int node_dist_from_repeat = prev_gene.distance_from_repeat + node_index_in_path;
                std::cout << "  Tail node " << idx << ": id=" << tn << ", dist_from_repeat=" << node_dist_from_repeat << std::endl;

                const int MAX_CANDIDATES = 5000;
                auto node_orfs = orf_finder.FindAllORFs(tn, search_min, search_max, min_orf_len, MAX_CANDIDATES);
                if (!node_orfs.empty()) std::cout << "    Found " << node_orfs.size() << " ORF(s) from tail node " << idx << std::endl;

                for (auto orf : node_orfs) {
                    int absolute_start = node_dist_from_repeat + orf.distance_from_repeat; // relative->absolute
                    // store absolute distance in ORFInfo
                    orf.distance_from_repeat = absolute_start;

                    // build uniqueness key using start:end
                    std::string key = std::to_string(orf.start_node) + ":" + std::to_string(orf.end_node);

                    // debug print for each ORF found from tail node
                    std::cout << "      ORF found: start=" << orf.start_node << " end=" << orf.end_node
                              << " dist=" << orf.distance_from_repeat << " len=" << orf.orf_length << std::endl;

                    if (seen_local.find(key) == seen_local.end()) {
                        if (used_start_end_keys.find(key) == used_start_end_keys.end()) {
                            all_orfs.push_back(orf);
                            seen_local.insert(key);
                            std::cout << "        Added ORF (unique)" << std::endl;
                        } else {
                            std::cout << "        Skipping ORF (already used in previous iterations): " << key << std::endl;
                        }
                    } else {
                        std::cout << "        Duplicate ORF within tail nodes, skipping: " << key << std::endl;
                    }
                }
            }
            std::cout << "  Total unique ORFs found from all tail nodes: " << all_orfs.size() << std::endl;
        }

        if (all_orfs.empty()) {
            std::cout << "No ORFs found. Stopping." << std::endl;
            break;
        }

        std::cout << "Found " << all_orfs.size() << " candidate ORF(s)" << std::endl;

        // Score all ORFs with profile sets
        ProfileScoringResult best_overall;
        best_overall.score = -std::numeric_limits<double>::infinity();
        ORFInfo best_orf;
        int best_absolute_distance = 0;

        for (const auto& orf : all_orfs) {
            // ORF distances in all_orfs are stored as ABSOLUTE distances from repeat
            int absolute_start = orf.distance_from_repeat;
            std::cout << "  Evaluating ORF: length=" << orf.orf_length << "bp, start at " << absolute_start << "bp from repeat" << std::endl;

            // get matching profiles (ProfileSize includes max_bp for dynamic depth)
            auto matching_profiles = HMMPicker::filterByLength(orf.orf_length);
            std::cout << "    Matching profiles count: " << matching_profiles.size() << std::endl;

            if (matching_profiles.empty()) {
                std::cout << "    No matching profiles, skipping." << std::endl;
                continue;
            }

            // Score ORF with all matching profiles (use their max_bp for beam depth)
            auto scoring_result = ScoreORFWithBestProfile(orf, matching_profiles);
            if (scoring_result.score > best_overall.score) {
                best_overall = scoring_result;
                best_orf = orf;
                best_absolute_distance = absolute_start;
                std::cout << "    *** New best score: " << scoring_result.score << " (profile: " << scoring_result.profile_name << ")" << std::endl;
            }
        }

        if (best_overall.score <= -std::numeric_limits<double>::infinity()) {
            std::cout << "All ORF candidates failed to score. Stopping." << std::endl;
            break;
        }

        std::cout << "Best ORF selected with score: " << best_overall.score << std::endl;

        // mark used start+end key to prevent future duplicate detection
        {
            std::string used_key = std::to_string(best_orf.start_node) + ":" + std::to_string(best_orf.end_node);
            used_start_end_keys.insert(used_key);
        }

        // store detected gene
        DetectedCasGene gene;
        gene.gene_name = best_overall.profile_name;
        gene.start_node = best_orf.start_node;
        gene.end_node = best_orf.end_node;
        gene.distance_from_repeat = best_absolute_distance;
        gene.gene_length = best_orf.orf_length;
        gene.score = best_overall.score;
        gene.amino_acids = best_overall.amino_acids;
        // Biological ORF path (start->stop) must come from ORFInfo
        gene.orf_node_path = best_orf.node_path;
        // HMM alignment node path from scoring
        gene.hmm_node_path = best_overall.hmm_node_path;

        res.genes.push_back(gene);

        res.total_length = best_absolute_distance + best_orf.orf_length;
        std::cout << "Total operon length so far: " << res.total_length << " bp" << std::endl;

        if (res.total_length >= max_total_length) {
            std::cout << "Reached maximum total length (" << max_total_length << " bp). Stopping." << std::endl;
            break;
        }

        // prepare next iteration
        if (is_first_orf) is_first_orf = false;
        search_min = 0;
        search_max = subsequent_search_distance;

        // update current_search_node and distance (search from END node)
        current_search_node = best_orf.end_node;
        current_distance_from_repeat = best_absolute_distance + best_orf.orf_length;
    }

    std::cout << "\n=== CAS Operon Detection Complete ===" << std::endl;
    std::cout << "Total genes detected: " << res.genes.size() << std::endl;
    std::cout << "Total operon length: " << res.total_length << " bp" << std::endl;

    return res;
}
