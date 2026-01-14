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
      beam_width(10), search_depth(500) {}

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

CasWorkflow::ProfileScoringResult CasWorkflow::ScoreORFWithProfile(const ORFInfo& orf, Profile* profile) {
    ProfileScoringResult res;
    res.score = -std::numeric_limits<double>::infinity();
    res.profile_name = profile->GetName();

    if (!profile) return res;

    // Use CasGeneDetector with this profile to perform beam search and score
    CasGeneDetector detector(sdbg, profile);
    auto paths = detector.BeamSearchAminoAcids(orf.start_node, beam_width, search_depth);
    if (paths.empty()) {
        return res;
    }

    // pick best path by total_score
    auto best_it = std::max_element(paths.begin(), paths.end(), [](const AminoAcidPathInfo&a, const AminoAcidPathInfo&b){ return a.total_score < b.total_score; });
    res.score = best_it->total_score;
    res.node_path = best_it->node_path;
    res.amino_acids = best_it->amino_acids;
    res.hmm_end_position = best_it->hmm_position;

    return res;
}

CasWorkflow::ProfileScoringResult CasWorkflow::ScoreORFWithBestProfile(const ORFInfo& orf, const std::vector<std::string>& profile_filenames) {
    ProfileScoringResult best;
    best.score = -std::numeric_limits<double>::infinity();
    best.profile_name = "";

    std::cout << "Scoring ORF (length=" << orf.orf_length << "bp) with " << profile_filenames.size() << " candidate profiles..." << std::endl;

    for (const auto& pf : profile_filenames) {
        std::cout << "  Loading profile: " << pf << std::endl;
        Profile* p = LoadProfile(pf);
        if (!p) {
            std::cerr << "    Failed to load profile: " << pf << std::endl;
            continue;
        }
        auto r = ScoreORFWithProfile(orf, p);
        r.profile_name = pf;
        std::cout << "    Profile " << pf << " score: " << r.score << std::endl;
        if (r.score > best.score) best = r;
        delete p;
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

    // track used start nodes to prevent cycles
    std::unordered_set<uint64_t> used_start_nodes;

    while (current_distance_from_repeat < max_total_length) {
        std::cout << "\n--- Iteration " << (res.genes.size() + 1) << " ---" << std::endl;
        std::cout << "Current distance from repeat: " << current_distance_from_repeat << " bp" << std::endl;
        std::cout << "Step 1: Finding ALL ORFs at distance [" << search_min << ", " << search_max << "] from current position..." << std::endl;

        int min_orf_len = HMMPicker::getMinimumORFLength();

        std::vector<ORFInfo> all_orfs;

        if (is_first_orf) {
            // find all ORFs in big window
            all_orfs = orf_finder.FindAllORFs(current_search_node, search_min, search_max, min_orf_len, std::numeric_limits<int>::max());
        } else {
            // search from last ~15 tail nodes
            const auto& prev_gene = res.genes.back();
            std::vector<uint64_t> tail_nodes;
            int num_tail = std::min(15, static_cast<int>(prev_gene.node_path.size()));
            int start_idx = static_cast<int>(prev_gene.node_path.size()) - num_tail;
            for (int i = start_idx; i < static_cast<int>(prev_gene.node_path.size()); ++i) {
                if (i >= 0) tail_nodes.push_back(prev_gene.node_path[i]);
            }
            if (tail_nodes.empty()) tail_nodes.push_back(current_search_node);

            std::unordered_set<uint64_t> seen_local;
            for (size_t idx = 0; idx < tail_nodes.size(); ++idx) {
                uint64_t tn = tail_nodes[idx];
                int node_index_in_path = start_idx + static_cast<int>(idx);
                // approximate dist from repeat to this tail node
                int node_dist_from_repeat = prev_gene.distance_from_repeat + node_index_in_path;
                std::cout << "  Tail node " << idx << ": id=" << tn << ", dist_from_repeat=" << node_dist_from_repeat << std::endl;

                auto node_orfs = orf_finder.FindAllORFs(tn, search_min, search_max, min_orf_len, std::numeric_limits<int>::max());
                if (!node_orfs.empty()) std::cout << "    Found " << node_orfs.size() << " ORF(s) from tail node " << idx << std::endl;

                for (auto orf : node_orfs) {
                    int absolute_start = node_dist_from_repeat + orf.distance_from_repeat; // debug: see logs
                    orf.distance_from_repeat = absolute_start;
                    if (seen_local.find(orf.start_node) == seen_local.end() && used_start_nodes.find(orf.start_node) == used_start_nodes.end()) {
                        all_orfs.push_back(orf);
                        seen_local.insert(orf.start_node);
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
            // compute absolute start distance
            int absolute_start = current_distance_from_repeat + orf.distance_from_repeat; // for first gene current_distance_from_repeat=0
            std::cout << "  Evaluating ORF: length=" << orf.orf_length << "bp, start at " << absolute_start << "bp from repeat" << std::endl;

            // get matching profiles
            auto matching_profiles = HMMPicker::filterByLength(orf.orf_length);
            std::cout << "    Matching profiles count: " << matching_profiles.size() << std::endl;
            std::vector<std::string> filenames;
            for (auto &ps : matching_profiles) filenames.push_back(ps.filename);

            if (filenames.empty()) {
                std::cout << "    No matching profiles, skipping." << std::endl;
                continue;
            }

            auto scoring_result = ScoreORFWithBestProfile(orf, filenames);
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

        // mark used start node
        used_start_nodes.insert(best_orf.start_node);

        // store detected gene
        DetectedCasGene gene;
        gene.gene_name = best_overall.profile_name;
        gene.start_node = best_orf.start_node;
        gene.end_node = best_orf.end_node;
        gene.distance_from_repeat = best_absolute_distance;
        gene.gene_length = best_orf.orf_length;
        gene.score = best_overall.score;
        gene.amino_acids = best_overall.amino_acids;
        gene.node_path = best_overall.node_path;

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
