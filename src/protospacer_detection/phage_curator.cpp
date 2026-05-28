#include "protospacer_detection/phage_curator.h"
#include <fstream>
#include <deque>
#include <cstdio> // for std::remove
#include <parallel_hashmap/phmap.h>
#include <functional>
#include <set>
#include <algorithm>
#include <limits>
#include <cassert>
#include <mutex>

// ---------- tiny helpers (do not change logic) ----------
namespace {
inline char safe_base_from_code(uint8_t code) {
    // Original logic used "ACGT"[code-1]; out-of-range now maps to 'N'
    static const char map[] = "NACGT"; // index 0..4
    if (code <= 4) return map[code];
    return 'N';
}
}

// constructor
PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {
    const bool validation = RevalidateAllNodesButSingleton();
    if (validation) {
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }
}

PhageCurator::PhageCurator(
    SDBG& sdbg,
    const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths,
    const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), grouped_paths(grouped_paths), cycles(cycles) {

    const bool validation = RevalidateAllNodesButSingleton();
    if (validation) {
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }

    float sum_mult = 0.0f;
    for (const auto& [id, cycle] : cycles) {
        for (const auto& path : cycle) {
            for (uint64_t node : path) {
                cycle_nodes.insert(node);
            }
        }
        // Compute avg_spacers for this cycle (preserve original accumulation semantics)
        phmap::flat_hash_set<uint64_t> least_frequent_node_multiplicities;
        for (const auto& path : cycle) {
            for (uint64_t node : path) {
                least_frequent_node_multiplicities.insert(node);
            }
        }
        for (uint64_t node : least_frequent_node_multiplicities) {
            sum_mult += sdbg.EdgeMultiplicity(node);
        }
        const size_t denom = least_frequent_node_multiplicities.size();
        avg_spacers[id] = denom ? (sum_mult / static_cast<float>(denom)) : 0.0f;
    }
}

std::string PhageCurator::_ReconstructPath(const std::vector<uint64_t>& path) {
    if (path.empty()) return std::string();
    std::string result = _FetchFirstNode(path.front());
    for (size_t i = 1; i < path.size(); ++i) {
        result += _FetchNodeLastBase(path[i]);
    }
    return result;
}

std::string PhageCurator::_FetchFirstNode(size_t node) {
    std::string label;
    const size_t k = sdbg.k();
    if (k == 0) return label;

    std::vector<uint8_t> seq(k, 0);
    sdbg.GetLabel(node, seq.data());

    // Original wrote: for (int i = k-1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    // Keep the same order and reverse afterward, but with safe indexing.
    for (size_t t = 0; t < k; ++t) {
        const size_t i = k - 1 - t;
        label.push_back(safe_base_from_code(seq[i]));
    }
    std::reverse(label.begin(), label.end());
    return label;
}

std::string PhageCurator::_FetchNodeLastBase(size_t node) {
    const size_t k = sdbg.k();
    if (k == 0) return std::string();

    std::vector<uint8_t> seq(k, 0);
    sdbg.GetLabel(node, seq.data());
    char base = safe_base_from_code(seq[0]);
    return std::string(1, base);
}

void PhageCurator::ReconstructPaths(std::vector<std::vector<uint64_t>> paths) {
    reconstructed_sequences.clear();
    reconstructed_sequences.reserve(paths.size());
    for (const auto& path : paths) {
        if (path.empty()) continue; // guard
        std::string result_path = _FetchFirstNode(path.front());
        for (size_t i = 1; i < path.size(); ++i) {
            result_path += _FetchNodeLastBase(path[i]);
        }
        reconstructed_sequences.push_back(std::move(result_path));
    }
}

void PhageCurator::WriteSequencesToFasta(const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    int count = 1;
    for (const auto& seq : reconstructed_sequences) {
        out << ">sequence_" << count << "\n" << seq << "\n";
        ++count;
    }
}

bool PhageCurator::RevalidateAllNodesButSingleton() {
    // If SDBG is not thread-safe, uncontrolled parallel SetValidEdge can crash.
    // Preserve parallel intent but guard calls with a critical section.
    #pragma omp parallel for
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        // We don’t rely on IsValidEdge() value to avoid TOCTTOU; SetValidEdge is assumed idempotent.
        #pragma omp critical(phage_curator_set_valid_edge)
        {
            if (!sdbg.IsValidEdge(node)) {
                sdbg.SetValidEdge(node);
            }
        }
    }
    return true;
}

std::vector<std::vector<uint64_t>> PhageCurator::BeamSearchPathsAvoiding(
    uint64_t start,
    int lower,
    int higher,
    const std::set<uint64_t>& forbidden,
    int beam_width,
    double min_mult,
    double max_mult,
    std::function<void(const std::vector<uint64_t>&)> path_callback) {

    std::vector<std::vector<uint64_t>> all_paths;
    if (!path_callback) {
        all_paths.reserve(static_cast<size_t>(std::max(0, beam_width)));
    }

    // thread-local pools
    static thread_local std::vector<std::vector<uint64_t>> path_pool;
    static thread_local std::vector<double> score_pool;
    static thread_local std::vector<uint64_t> current_pool;

    path_pool.clear();
    score_pool.clear();
    current_pool.clear();

    const size_t reserve_hint = (higher > 0) ? static_cast<size_t>(higher + 64) : 128;
    if (path_pool.capacity() == 0) path_pool.reserve(reserve_hint);
    if (score_pool.capacity() == 0) score_pool.reserve(reserve_hint);
    if (current_pool.capacity() == 0) current_pool.reserve(reserve_hint);

    auto comp = std::greater<std::pair<double, size_t>>();
    std::set<std::pair<double, size_t>, decltype(comp)> beam_set(comp);

    size_t unique_id = 0;

    // Initial multiplicity check
    const double initial_mult = sdbg.EdgeMultiplicity(start);
    if (initial_mult <= 1.0 || initial_mult < min_mult || initial_mult > max_mult) {
        return all_paths;
    }

    path_pool.emplace_back(std::vector<uint64_t>{start});
    path_pool.back().reserve((higher > 0) ? static_cast<size_t>(higher + 1) : path_pool.back().capacity());
    score_pool.push_back(initial_mult);
    current_pool.push_back(start);
    beam_set.insert({initial_mult, unique_id++});

    while (!beam_set.empty()) {
        // Pop best
        const auto it = beam_set.begin();
        const double score = it->first;
        const size_t id = it->second;
        beam_set.erase(it);

        // Do NOT bind a reference to an element of a growing vector; use index each time.
        const std::vector<uint64_t>& path_ref = path_pool[id];
        const uint64_t v = current_pool[id];
        const int current_depth = static_cast<int>(path_ref.size()) - 1;

        if (current_depth >= lower && current_depth <= higher) {
            if (path_callback) {
                path_callback(path_ref);
            } else {
                all_paths.push_back(path_ref);
            }
            // Preserve original behavior: do not expand after collection
            continue;
        }

        if (sdbg.EdgeOutdegreeZero(v)) {
            continue;
        }

        const int outdegree = sdbg.EdgeOutdegree(v);
        if (outdegree <= 0 || !sdbg.IsValidEdge(v)) {
            continue;
        }

        // SAFETY: cap iteration by buffer capacity
        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT] = {0, 0, 0, 0};
        const int flag = sdbg.OutgoingEdges(v, neighbors);
        if (flag == -1) {
            continue;
        }
        const int n_iter = std::min(outdegree, MAX_EDGE_COUNT);

        for (int i = 0; i < n_iter; ++i) {
            const uint64_t neighbor = neighbors[i];

            // Cycle check (on current path)
            if (std::find(path_ref.begin(), path_ref.end(), neighbor) != path_ref.end()) {
                continue;
            }

            // Forbidden check
            if (forbidden.find(neighbor) != forbidden.end() && neighbor != start) {
                continue;
            }

            // Multiplicity guard
            const double mult = sdbg.EdgeMultiplicity(neighbor);
            if (mult <= 1.0 || mult < min_mult || mult > max_mult) {
                continue;
            }

            // Copy base path to avoid dangling references upon vector growth
            std::vector<uint64_t> new_path = path_ref;
            new_path.push_back(neighbor);

            double new_score = (score * std::max(0, current_depth) + mult) /
                               static_cast<double>(std::max(1, current_depth + 1));
            uint64_t current = neighbor;

            // Extend along simple path (preserving original loop’s semantics)
            while (true) {
                const bool depth_limit = (static_cast<int>(new_path.size()) - 1) >= higher;
                if (depth_limit) break;

                const uint64_t next = sdbg.NextSimplePathEdge(current);
                if (next == SDBG::kNullID) break;

                // Cycle check
                if (std::find(new_path.begin(), new_path.end(), next) != new_path.end()) {
                    break;
                }
                // Forbidden check
                if (forbidden.find(next) != forbidden.end() && next != start) {
                    break;
                }
                // Multiplicity check
                const double next_mult = sdbg.EdgeMultiplicity(next);
                if (next_mult <= 1.0 || next_mult < min_mult || next_mult > max_mult) {
                    break;
                }

                new_path.push_back(next);
                const int new_depth = static_cast<int>(new_path.size()) - 1;
                new_score = (new_score * (new_depth - 1) + next_mult) / static_cast<double>(new_depth);
                current = next;
            }

            // Append new entries atomically to the pools
            const size_t new_id = unique_id++;
            path_pool.emplace_back(std::move(new_path));
            current_pool.push_back(current);
            score_pool.push_back(new_score);

            beam_set.insert({new_score, new_id});

            // Prune the lowest score if over beam width
            if (beam_width > 0 && beam_set.size() > static_cast<size_t>(beam_width)) {
                auto prune_it = beam_set.end();
                --prune_it;
                beam_set.erase(prune_it);
            }
        }
    }

    return all_paths;
}

// Backward beam search: follows IncomingEdges instead of OutgoingEdges.
// Returned paths grow leftward: path[0]=start, path[1]=prev1, ...
// Caller must reverse path[1:] and prepend to the protospacer to get final order.
std::vector<std::vector<uint64_t>> PhageCurator::BeamSearchBackwardAvoiding(
    uint64_t start,
    int lower,
    int higher,
    const std::set<uint64_t>& forbidden,
    int beam_width,
    double min_mult,
    double max_mult) {

    std::vector<std::vector<uint64_t>> all_paths;

    static thread_local std::vector<std::vector<uint64_t>> path_pool;
    static thread_local std::vector<double> score_pool;
    static thread_local std::vector<uint64_t> current_pool;
    path_pool.clear(); score_pool.clear(); current_pool.clear();

    auto comp = std::greater<std::pair<double, size_t>>();
    std::set<std::pair<double, size_t>, decltype(comp)> beam_set(comp);
    size_t unique_id = 0;

    const double initial_mult = sdbg.EdgeMultiplicity(start);
    if (initial_mult <= 1.0 || initial_mult < min_mult || initial_mult > max_mult)
        return all_paths;

    path_pool.push_back({start});
    score_pool.push_back(initial_mult);
    current_pool.push_back(start);
    beam_set.insert({initial_mult, unique_id++});

    while (!beam_set.empty()) {
        const auto it = beam_set.begin();
        const double score = it->first;
        const size_t id = it->second;
        beam_set.erase(it);

        const std::vector<uint64_t>& path_ref = path_pool[id];
        const uint64_t v = current_pool[id];
        const int current_depth = static_cast<int>(path_ref.size()) - 1;

        if (current_depth >= lower && current_depth <= higher) {
            all_paths.push_back(path_ref);
            continue;
        }

        const int indegree = sdbg.EdgeIndegree(v);
        if (indegree <= 0 || !sdbg.IsValidEdge(v)) continue;

        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT] = {0,0,0,0};
        const int flag = sdbg.IncomingEdges(v, neighbors);
        if (flag == -1) {
            continue;
        }
        const int n_iter = std::min(flag, MAX_EDGE_COUNT);

        for (int i = 0; i < n_iter; ++i) {
            const uint64_t neighbor = neighbors[i];

            if (std::find(path_ref.begin(), path_ref.end(), neighbor) != path_ref.end()) continue;
            if (forbidden.find(neighbor) != forbidden.end()) continue;

            const double mult = sdbg.EdgeMultiplicity(neighbor);
            if (mult <= 1.0 || mult < min_mult || mult > max_mult) continue;

            std::vector<uint64_t> new_path = path_ref;
            new_path.push_back(neighbor);
            double new_score = (score * std::max(0, current_depth) + mult) /
                               static_cast<double>(std::max(1, current_depth + 1));
            uint64_t current = neighbor;

            // Extend along simple backward path
            while (true) {
                if (static_cast<int>(new_path.size()) - 1 >= higher) break;
                const uint64_t prev = sdbg.PrevSimplePathEdge(current);
                if (prev == SDBG::kNullID) break;
                if (std::find(new_path.begin(), new_path.end(), prev) != new_path.end()) break;
                if (forbidden.find(prev) != forbidden.end()) break;
                const double pm = sdbg.EdgeMultiplicity(prev);
                if (pm <= 1.0 || pm < min_mult || pm > max_mult) break;
                new_path.push_back(prev);
                const int nd = static_cast<int>(new_path.size()) - 1;
                new_score = (new_score * (nd - 1) + pm) / static_cast<double>(nd);
                current = prev;
            }

            const size_t new_id = unique_id++;
            path_pool.push_back(std::move(new_path));
            score_pool.push_back(new_score);
            current_pool.push_back(current);
            beam_set.insert({new_score, new_id});

            if (beam_width > 0 && beam_set.size() > static_cast<size_t>(beam_width)) {
                auto prune = beam_set.end(); --prune;
                beam_set.erase(prune);
            }
        }
    }
    return all_paths;
}

// Find paths of form: incoming_path --> spacer_nodes --> outgoing_path
// where total sequence length <= max_total_bp.
std::map<std::string, std::vector<std::string>>
PhageCurator::FindContiguousPathsFromGroupedPaths(int max_total_bp, const std::string& filename, int beam_width) {
    std::ofstream out(filename, std::ios::app);
    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }

    std::map<std::string, std::vector<std::string>> consensus_map;
    int path_count = 0;
    const int k = static_cast<int>(sdbg.k());

    for (const auto& [group_id, cycle_map] : grouped_paths) {
        for (const auto& [cycle_id, paths_vec] : cycle_map) {
            for (const auto& path : paths_vec) {
                if (path.empty()) continue;

                const int spacer_bp = k - 1 + static_cast<int>(path.size());
                const int remaining_bp = max_total_bp - spacer_bp;
                if (remaining_bp <= 0) continue;

                auto make_mult_range = [&](uint64_t node, double& mn, double& mx) {
                    double m = sdbg.EdgeMultiplicity(node);
                    mn = std::max(1.0, 0.1 * m);
                    mx = 5.0 * m;
                    if (mx < mn) mx = mn * 50.0;
                };

                double bk_min, bk_max, fw_min, fw_max;
                make_mult_range(path.front(), bk_min, bk_max);
                make_mult_range(path.back(),  fw_min, fw_max);

                const int left_budget = remaining_bp / 2;
                const int right_budget = remaining_bp - left_budget;

                auto back_paths = BeamSearchBackwardAvoiding(
                    path.front(), left_budget, left_budget,
                    cycle_nodes, beam_width, bk_min, bk_max);

                auto fwd_paths = BeamSearchPathsAvoiding(
                    path.back(), right_budget, right_budget,
                    cycle_nodes, beam_width, fw_min, fw_max, nullptr);

                if (back_paths.empty() || fwd_paths.empty()) continue;

                const std::vector<uint64_t> empty_stub;
                const auto& bp = back_paths.empty() ? empty_stub : back_paths[0];
                const auto& fp = fwd_paths.empty()  ? empty_stub : fwd_paths[0];

                std::vector<uint64_t> combined;
                if (bp.size() > 1) {
                    for (int i = static_cast<int>(bp.size()) - 1; i >= 1; --i)
                        combined.push_back(bp[i]);
                }
                for (uint64_t n : path) combined.push_back(n);
                if (fp.size() > 1) {
                    for (size_t i = 1; i < fp.size(); ++i)
                        combined.push_back(fp[i]);
                }

                if (combined.empty()) continue;

                std::string seq = _FetchFirstNode(combined.front());
                for (size_t i = 1; i < combined.size(); ++i)
                    seq += _FetchNodeLastBase(combined[i]);

                if (static_cast<int>(seq.size()) > max_total_bp) continue;

                out << ">contiguous_path_" << path_count << "\n" << seq << "\n";
                ++path_count;
            }
        }
    }

    std::cout << "\nSaved " << path_count << " contiguous paths in " << filename << std::endl;
    out.close();
    return consensus_map;
}

// New function using beam search, mirroring FindQualityPathsDLSFromGroupedPaths
std::map<std::string, std::vector<std::string>>
PhageCurator::FindQualityPathsBeamSearchFromGroupedPaths(int min_length, int max_length, const std::string& filename, int beam_width) {
    std::ofstream out(filename, std::ios::app);  // Append mode to add to file
    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }

    std::map<std::string, std::vector<std::string>> consensus_map;
    int path_count = 0;

    struct counts {
        int first;
        int second;
        int third;
    };
    counts counts_by_nodes = {0, 0, 0};

    for (const auto& [group_id, cycle_map] : grouped_paths) {
        std::vector<std::string> quality_paths;
        counts_by_nodes.first += 1;
        std::cout << "Processed groups: " << counts_by_nodes.first << "\r";

        for (const auto& [cycle_id, paths_vec] : cycle_map) {
            (void)cycle_id; // unused in current logic
            counts_by_nodes.second += 1;

            for (const auto& path : paths_vec) {
                counts_by_nodes.third += 1;

                if (path.empty()) continue; // guard

                // Start from last node as in original code
                const uint64_t start = path.back();

                // min/max multiplicity window around start
                double min_mult = 0.1 * sdbg.EdgeMultiplicity(start);
                if (min_mult < 1.0) min_mult = 1.0;
                double max_mult = 5.0 * sdbg.EdgeMultiplicity(start);
                if (max_mult < min_mult) max_mult = min_mult * 50.0;

                auto extended_paths = BeamSearchPathsAvoiding(
                    start, min_length, max_length, cycle_nodes, beam_width, min_mult, max_mult, nullptr);

                if (extended_paths.empty()) continue;

                // NOTE: keep original parameter order/types to avoid changing logic
                std::vector<std::vector<uint64_t>> best_paths =
                    GetTopPathsFromBeamPaths(extended_paths, /*max*/ static_cast<int>(min_mult), /*min*/ static_cast<int>(max_mult), /*top_n*/ 2);

                if (best_paths.empty()) continue;

                for (const auto& ext_path : best_paths) {
                    if (ext_path.empty()) continue;

                    std::string result_path = _FetchFirstNode(ext_path.front());
                    for (size_t i = 1; i < ext_path.size(); ++i) {
                        result_path += _FetchNodeLastBase(ext_path[i]);
                    }

                    // Keep original output side-effect and counter behavior
                    out << ">quality_path_" << path_count << "\n" << result_path << "\n";
                    ++path_count;
                    // quality_paths.push_back(result_path); // commented in original; keep as-is
                }
            }
        }

        std::cout << "Processed counts:" << counts_by_nodes.first << ","
                  << counts_by_nodes.second << "," << counts_by_nodes.third << ";";

        std::string group_id_str = _FetchFirstNode(group_id);
        consensus_map[group_id_str] = quality_paths;
    }

    std::cout << "\nSaved in " << filename << std::endl;
    out.close();
    return consensus_map;
}

std::string PhageCurator::ComputeConsensusForCurrentGroup(std::vector<std::string> sequences) {
    if (sequences.empty()) return std::string();

    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }
    auto consensus = graph.GenerateConsensus();
    return consensus;
}

std::vector<std::vector<uint64_t>> PhageCurator::GetTopPathsFromBeamPaths(
    const std::vector<std::vector<uint64_t>>& beam_paths,
    int max,
    int min,
    size_t top_n) {

    struct PathScore {
        std::vector<uint64_t> path;
        double avg_mult;
        double delta;
    };

    std::vector<PathScore> scored_paths;
    scored_paths.reserve(beam_paths.size());

    for (const auto& path : beam_paths) {
        if (path.empty()) continue;

        double total = 0.0;
        for (uint64_t node : path) {
            total += sdbg.EdgeMultiplicity(node);
        }
        const double avg = (path.size() ? (total / static_cast<double>(path.size())) : 0.0);

        double delta = 0.0;
        if (avg < static_cast<double>(min))
            delta = static_cast<double>(min) - avg;
        else if (avg > static_cast<double>(max))
            delta = avg - static_cast<double>(max);
        else
            delta = 0.0;

        scored_paths.push_back({path, avg, delta});
    }

    std::sort(scored_paths.begin(), scored_paths.end(),
              [](const PathScore& a, const PathScore& b) {
                  return a.delta < b.delta;
              });

    const size_t keep = std::min(top_n, scored_paths.size());
    std::vector<std::vector<uint64_t>> top_paths;
    top_paths.reserve(keep);
    for (size_t i = 0; i < keep; ++i) {
        top_paths.push_back(scored_paths[i].path);
    }
    return top_paths;
}
