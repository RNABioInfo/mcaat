#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <functional>
#include <chrono>
#include <iomanip>
#include <ctime>
#include "spoa/spoa.hpp"
#include "core/settings.h"

namespace fs = std::filesystem;

/**
 * PostProcessor: Deduplicate cycles and output CRISPR arrays.
 *
 * Pipeline:
 * 1-5: ORIGINAL — glue, group by 23-mer, unanimous extend, extract spacers, dedup
 * 6:   Cluster similar repeats by edit distance, SPOA consensus (front repeat)
 * 7:   Per array: detect tail repeat in spacer suffixes, validate it matches
 *       a prefix of front consensus, SPOA consensus on tails, strip + prepend
 * 8:   Sanity filter on repeat/spacer dimensions
 * 9:   Output
 */
class PostProcessor {
private:
    Settings& settings;

    static constexpr int GROUP_KMER_SIZE = 23;
    static constexpr int DEDUP_KMER_SIZE = 23;
    static constexpr int TAIL_KMER_SIZE = 22;
    static constexpr int MAX_LINES_PER_FILE = 10000;

    // SPOA merge parameters
    static constexpr int MAX_EDIT_DIST = 4;

    // Tail detection parameters
    static constexpr int MAX_TAIL_SCAN = 40;
    static constexpr double TAIL_AGREE_THRESH = 0.6;
    static constexpr int TAIL_SUSTAIN_WINDOW = 3;
    static constexpr int MIN_TAIL_LEN = 5;
    static constexpr int MIN_SPACER_AFTER_TRIM = 15;

    // Tail validation: max edit distance between tail_consensus and front_consensus prefix
    static constexpr int MAX_TAIL_PREFIX_EDIT_DIST = 2;

    // Sanity filters
    static constexpr int MAX_REPEAT_LEN = 55;
    static constexpr int MIN_REPEAT_LEN = 20;
    static constexpr int MIN_SPACER_LEN = 20;
    static constexpr double MIN_MEDIAN_SPACER_REPEAT_RATIO = 0.5;
    // Maximum allowed similarity (0–1) between a spacer and the consensus repeat.
    // Real CRISPR spacers are foreign DNA — completely different from the repeat.
    // Tandem repeat false positives have spacers that look like the repeat.
    static constexpr double MAX_SPACER_REPEAT_SIMILARITY = 0.5;

    struct ConsensusArray {
        std::string consensus;
        std::vector<std::pair<std::string, std::string>> entries;
    };
    std::vector<ConsensusArray> consensus_arrays;

    std::map<std::string, std::vector<std::string>> crispr_arrays;

    struct SpacerData {
        std::string spacer;
        std::string repeat;
        std::string cycle;
    };

    std::vector<int> parent;

    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int x, int y) {
        int px = find(x), py = find(y);
        if (px != py) parent[px] = py;
    }

    std::unordered_set<std::string> get_kmers(const std::string& seq, int k) {
        std::unordered_set<std::string> kmers;
        if ((int)seq.size() >= k) {
            for (size_t i = 0; i <= seq.size() - k; ++i)
                kmers.insert(seq.substr(i, k));
        }
        return kmers;
    }

    std::string glue_kmers(const std::string& line) {
        std::istringstream iss(line);
        std::vector<std::string> kmers;
        std::string kmer;
        while (iss >> kmer) kmers.push_back(kmer);
        if (kmers.empty()) return "";

        std::string result = kmers[0];
        for (size_t i = 1; i < kmers.size(); ++i) {
            const std::string& next = kmers[i];
            bool found = false;
            for (int k = (int)result.size(); k > 0; --k) {
                if ((int)next.size() >= k &&
                    result.substr(result.size() - k) == next.substr(0, k)) {
                    result += next.substr(k);
                    found = true;
                    break;
                }
            }
            if (!found) result += next;
        }
        return result;
    }

    int edit_distance(const std::string& a, const std::string& b, int max_dist) {
        int m = (int)a.size(), n = (int)b.size();
        if (std::abs(m - n) > max_dist) return max_dist + 1;

        std::vector<int> prev(n + 1), curr(n + 1);
        for (int j = 0; j <= n; ++j) prev[j] = j;

        for (int i = 1; i <= m; ++i) {
            curr[0] = i;
            int row_min = curr[0];
            for (int j = 1; j <= n; ++j) {
                int cost = (a[i-1] == b[j-1]) ? 0 : 1;
                curr[j] = std::min({prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost});
                row_min = std::min(row_min, curr[j]);
            }
            if (row_min > max_dist) return max_dist + 1;
            std::swap(prev, curr);
        }
        return prev[n];
    }

    static std::string revcomp(const std::string& s) {
        static const auto ct = [](){
            std::array<char,256> t{};
            t['A']='T'; t['T']='A'; t['C']='G'; t['G']='C';
            t['a']='T'; t['t']='A'; t['c']='G'; t['g']='C';
            t['N']='N'; t['n']='N';
            return t;
        }();
        std::string rc(s.size(), 'N');
        for (size_t i = 0; i < s.size(); ++i)
            rc[s.size()-1-i] = ct[(unsigned char)s[i]];
        return rc;
    }

    // Canonical form of a sequence: min(seq, RC(seq)).
    // Ensures that a sequence and its reverse complement map to the same key,
    // so the graph finding both strands does not double-count arrays.
    static std::string canonical_seq(const std::string& s) {
        std::string rc = revcomp(s);
        return (s <= rc) ? s : rc;
    }

    std::string spoa_consensus(const std::vector<std::pair<std::string, int>>& seq_weights) {
        if (seq_weights.empty()) return "";
        if (seq_weights.size() == 1) return seq_weights[0].first;

        auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 5, -4, -8, -6);

        spoa::Graph graph;

        for (const auto& [seq, weight] : seq_weights) {
            for (int w = 0; w < weight; ++w) {
                auto alignment = alignment_engine->Align(seq, graph);
                graph.AddAlignment(alignment, seq);
            }
        }

        return graph.GenerateConsensus();
    }

    /**
     * Detect tail repeat length by reverse majority voting on spacer suffixes.
     */
    int detect_tail_length(const std::vector<std::string>& spacers) {
        if (spacers.size() < 2) return 0;

        int min_len = (int)spacers[0].size();
        for (const auto& s : spacers)
            min_len = std::min(min_len, (int)s.size());

        int scan_limit = std::min(MAX_TAIL_SCAN, min_len - MIN_SPACER_AFTER_TRIM);
        if (scan_limit < MIN_TAIL_LEN) return 0;

        int tail_len = 0;
        int low_streak = 0;

        for (int rpos = 0; rpos < scan_limit; ++rpos) {
            int counts[4] = {0, 0, 0, 0};
            int total = 0;
            for (const auto& s : spacers) {
                int pos = (int)s.size() - 1 - rpos;
                switch (s[pos]) {
                    case 'A': counts[0]++; total++; break;
                    case 'C': counts[1]++; total++; break;
                    case 'G': counts[2]++; total++; break;
                    case 'T': counts[3]++; total++; break;
                }
            }
            if (total == 0) break;

            int best = *std::max_element(counts, counts + 4);
            double agreement = (double)best / total;

            if (agreement >= TAIL_AGREE_THRESH) {
                low_streak = 0;
                tail_len = rpos + 1;
            } else {
                low_streak++;
                if (low_streak >= TAIL_SUSTAIN_WINDOW) break;
            }
        }

        return (tail_len >= MIN_TAIL_LEN) ? tail_len : 0;
    }

    /**
     * Validate that tail_consensus approximately matches a prefix of front_consensus.
     * This is the key biological constraint: the tail is the beginning of the next
     * repeat occurrence in the cycle, so it must look like the front repeat.
     */
    bool validate_tail_matches_front(const std::string& tail, const std::string& front) {
        if (tail.empty() || front.empty()) return false;

        int check_len = std::min((int)tail.size(), (int)front.size());
        std::string tail_portion = tail.substr(tail.size() - check_len);
        std::string front_prefix = front.substr(0, check_len);

        int dist = edit_distance(tail_portion, front_prefix, MAX_TAIL_PREFIX_EDIT_DIST + 1);
        return dist <= MAX_TAIL_PREFIX_EDIT_DIST;
    }

    /**
     * Compute median of a vector of ints.
     */
    double median(std::vector<int>& vals) {
        if (vals.empty()) return 0;
        std::sort(vals.begin(), vals.end());
        int n = (int)vals.size();
        if (n % 2 == 0) return (vals[n/2 - 1] + vals[n/2]) / 2.0;
        return vals[n/2];
    }

    /**
     * Rotate s to start at the position of its highest-frequency GROUP_KMER_SIZE-mer
     * as measured across all cycles (supplied via freq_map).
     */
    std::string freq_rotation(const std::string& s,
                              const std::unordered_map<std::string,int>& freq_map) const {
        int n = (int)s.size();
        if (n < GROUP_KMER_SIZE) return s;
        int best_pos = 0, best_freq = 0;
        for (int i = 0; i + GROUP_KMER_SIZE <= n; ++i) {
            auto it = freq_map.find(s.substr(i, GROUP_KMER_SIZE));
            int f = (it != freq_map.end()) ? it->second : 0;
            if (f > best_freq) { best_freq = f; best_pos = i; }
        }
        if (best_pos == 0) return s;
        return s.substr(best_pos) + s.substr(0, best_pos);
    }

    /**
     * Returns true if seq is a tandem repeat of a shorter unit (period ≤ len/2)
     * with at most 15% mismatch. CRISPR repeats are not tandem repeats.
     */
    bool is_tandem_repeat(const std::string& seq) const {
        int n = (int)seq.size();
        for (int p = 2; p <= n / 2; ++p) {
            int mismatches = 0;
            for (int i = p; i < n; ++i)
                if (seq[i] != seq[i % p]) ++mismatches;
            if ((double)mismatches / n < 0.15) return true;
        }
        return false;
    }

public:
    explicit PostProcessor(Settings& s) : settings(s) {}

    std::map<std::string, std::vector<std::string>> getSystems() const {
        return crispr_arrays;
    }

    void run_analysis() {
        std::string input_path = settings.output_folder + "/cycles.txt";
        std::string output_dir = settings.output_folder;
        fs::create_directories(output_dir);

        // ==========================================================
        // ORIGINAL Step 1: Read all cycles and glue k-mers
        // ==========================================================
        std::vector<std::string> cycles;
        {
            // First pass: read and glue k-mers into raw sequences
            std::vector<std::string> raw_cycles;
            {
                std::ifstream in(input_path);
                if (!in.is_open()) {
                    std::cerr << "Error: Cannot open cycles file: " << input_path << std::endl;
                    return;
                }
                std::string line;
                while (std::getline(in, line)) {
                    if (!line.empty()) {
                        line.erase(0, line.find_first_not_of(" \t\r\n"));
                        line.erase(line.find_last_not_of(" \t\r\n") + 1);
                        if (!line.empty()) {
                            std::string glued = glue_kmers(line);
                            if (!glued.empty()) raw_cycles.push_back(glued);
                        }
                    }
                }
            }

            // Build global k-mer frequency map — repeat k-mers are high-frequency.
            std::unordered_map<std::string,int> kmer_freq;
            kmer_freq.reserve(raw_cycles.size() * 4);
            for (const auto& seq : raw_cycles)
                for (int i = 0; i + GROUP_KMER_SIZE <= (int)seq.size(); ++i)
                    kmer_freq[seq.substr(i, GROUP_KMER_SIZE)]++;

            cycles.reserve(raw_cycles.size());
            for (const auto& seq : raw_cycles)
                cycles.push_back(freq_rotation(seq, kmer_freq));
        }
        std::cout << "  ▸ Loaded " << cycles.size() << " cycles" << std::endl;

        // ==========================================================
        // ORIGINAL Step 2: Group by first 23-mer
        // ==========================================================
        std::unordered_map<std::string, std::vector<std::string>> groups;
        for (const auto& seq : cycles) {
            if ((int)seq.size() >= GROUP_KMER_SIZE)
                groups[seq.substr(0, GROUP_KMER_SIZE)].push_back(seq);
        }
        std::cout << "  ▸ Grouped into " << groups.size() << " repeat groups" << std::endl;

        // ==========================================================
        // Step 3 & 4: Unanimous extend, extract spacers
        // ==========================================================
        std::vector<SpacerData> spacer_data;

        for (const auto& [kmer, group_cycles] : groups) {
            int t = 0;
            while (true) {
                std::unordered_set<char> chars;
                bool valid = true;
                for (const auto& c : group_cycles) {
                    int pos = GROUP_KMER_SIZE + t;
                    if (pos >= (int)c.size() - TAIL_KMER_SIZE) {
                        valid = false;
                        break;
                    }
                    chars.insert(c[pos]);
                }
                if (!valid || chars.size() != 1) break;
                ++t;
            }

            int repeat_len = GROUP_KMER_SIZE + t;

            for (const auto& c : group_cycles) {
                if ((int)c.size() > repeat_len + TAIL_KMER_SIZE) {
                    std::string repeat = c.substr(0, repeat_len);
                    std::string spacer = c.substr(repeat_len, c.size() - repeat_len - TAIL_KMER_SIZE);
                    if ((int)spacer.size() >= MIN_SPACER_LEN) {
                        spacer_data.push_back({spacer, repeat, c});
                    }
                }
            }
        }

        std::cout << "  ▸ Extracted " << spacer_data.size() << " spacers" << std::endl;
        if (spacer_data.empty()) {
            std::cout << "  No spacers found." << std::endl;
            return;
        }

        // ==========================================================
        // ORIGINAL Step 5: Deduplicate spacers via k-mer union-find
        // ==========================================================
        std::unordered_map<std::string, std::vector<int>> kmer_index;
        for (size_t idx = 0; idx < spacer_data.size(); ++idx) {
            auto kmers = get_kmers(spacer_data[idx].spacer, DEDUP_KMER_SIZE);
            for (const auto& km : kmers)
                kmer_index[km].push_back(static_cast<int>(idx));
        }
        std::cout << "  ▸ Index: " << kmer_index.size() << " unique " << DEDUP_KMER_SIZE << "-mers" << std::endl;

        parent.resize(spacer_data.size());
        std::iota(parent.begin(), parent.end(), 0);

        for (const auto& [km, indices] : kmer_index) {
            for (size_t i = 1; i < indices.size(); ++i)
                unite(indices[0], indices[i]);
        }

        std::unordered_map<int, std::vector<int>> clusters;
        for (size_t idx = 0; idx < spacer_data.size(); ++idx)
            clusters[find(static_cast<int>(idx))].push_back(static_cast<int>(idx));

        std::cout << "  ▸ Spacer clusters: " << clusters.size() << std::endl;

        // Pick longest spacer per cluster, group by repeat
        std::map<std::string, std::vector<std::string>> repeat_to_spacers;
        for (const auto& [root, indices] : clusters) {
            int best_idx = indices[0];
            for (int idx : indices) {
                if (spacer_data[idx].spacer.size() > spacer_data[best_idx].spacer.size())
                    best_idx = idx;
            }
            const auto& sd = spacer_data[best_idx];
            repeat_to_spacers[sd.repeat].push_back(sd.spacer);
        }

        // Remove singletons
        for (auto it = repeat_to_spacers.begin(); it != repeat_to_spacers.end(); ) {
            if (it->second.size() < 2)
                it = repeat_to_spacers.erase(it);
            else
                ++it;
        }

        std::cout << "  ▸ Arrays before merge: " << repeat_to_spacers.size() << std::endl;

        // ==========================================================
        // Step 6: Cluster similar repeats by edit distance,
        //         SPOA consensus on front repeats
        // ==========================================================
        std::vector<std::string> all_repeats;
        for (const auto& [rep, _] : repeat_to_spacers)
            all_repeats.push_back(rep);

        int nr = (int)all_repeats.size();

        std::vector<int> rep_parent(nr);
        std::iota(rep_parent.begin(), rep_parent.end(), 0);

        std::function<int(int)> rep_find = [&](int x) -> int {
            if (rep_parent[x] != x) rep_parent[x] = rep_find(rep_parent[x]);
            return rep_parent[x];
        };
        auto rep_unite = [&](int x, int y) {
            int px = rep_find(x), py = rep_find(y);
            if (px != py) rep_parent[px] = py;
        };

        for (int i = 0; i < nr; ++i) {
            for (int j = i + 1; j < nr; ++j) {
                if (edit_distance(all_repeats[i], all_repeats[j], MAX_EDIT_DIST) <= MAX_EDIT_DIST)
                    rep_unite(i, j);
            }
        }

        std::unordered_map<int, std::vector<int>> repeat_clusters;
        for (int i = 0; i < nr; ++i)
            repeat_clusters[rep_find(i)].push_back(i);

        std::cout << "  ▸ Repeat groups (post-SPOA merge): " << repeat_clusters.size() << std::endl;

        // ==========================================================
        // Step 7: Per consensus group — front SPOA + validated tail SPOA
        // ==========================================================
        consensus_arrays.clear();
        crispr_arrays.clear();

        for (const auto& [croot, rep_indices] : repeat_clusters) {
            std::vector<std::pair<std::string, int>> weighted_repeats;
            std::vector<std::pair<std::string, std::string>> merged_entries;

            for (int ri : rep_indices) {
                const std::string& rep = all_repeats[ri];
                const auto& spacers = repeat_to_spacers[rep];
                weighted_repeats.emplace_back(rep, (int)spacers.size());
                for (const auto& sp : spacers)
                    merged_entries.emplace_back(rep, sp);
            }

            if (merged_entries.size() < 2) continue;

            // --- Front consensus via SPOA ---
            std::string front_consensus = spoa_consensus(weighted_repeats);

            // --- Tail detection + validation ---
            std::vector<std::string> all_spacers;
            for (const auto& [rv, sp] : merged_entries)
                all_spacers.push_back(sp);

            int tail_len = detect_tail_length(all_spacers);
            std::string tail_consensus;
            bool tail_valid = false;

            if (tail_len > 0) {
                // Extract tail portions, SPOA consensus
                std::vector<std::pair<std::string, int>> tail_seqs;
                for (const auto& sp : all_spacers) {
                    std::string tail = sp.substr(sp.size() - tail_len);
                    tail_seqs.emplace_back(tail, 1);
                }
                tail_consensus = spoa_consensus(tail_seqs);

                // VALIDATE: tail must match a prefix of front_consensus
                tail_valid = validate_tail_matches_front(tail_consensus, front_consensus);

                if (tail_valid) {
                    // Trim tail from spacers
                    std::vector<std::pair<std::string, std::string>> trimmed_entries;
                    for (const auto& [rv, sp] : merged_entries) {
                        int trim = std::min(tail_len, (int)sp.size() - MIN_SPACER_AFTER_TRIM);
                        if (trim > 0) {
                            std::string trimmed = sp.substr(0, sp.size() - trim);
                            if ((int)trimmed.size() >= MIN_SPACER_AFTER_TRIM) {
                                trimmed_entries.emplace_back(rv, trimmed);
                                continue;
                            }
                        }
                        trimmed_entries.push_back({rv, sp});
                    }
                    merged_entries = std::move(trimmed_entries);
                } else {
                    tail_consensus.clear();
                }
            }

            if (merged_entries.size() < 2) continue;

            // Full consensus = tail_consensus + front_consensus
            std::string full_consensus = tail_consensus + front_consensus;

            // --- Sanity filter: repeat length ---
            if ((int)full_consensus.size() > MAX_REPEAT_LEN ||
                (int)full_consensus.size() < MIN_REPEAT_LEN) {
                std::cout << "  [filtered] repeat len=" << full_consensus.size()
                          << "bp (allowed " << MIN_REPEAT_LEN << "–" << MAX_REPEAT_LEN << "): "
                          << full_consensus.substr(0, 30) << "\n";
                continue;
            }
            // --- Sanity filter: median spacer / repeat ratio ---
            std::vector<int> spacer_lens;
            for (const auto& [rv, sp] : merged_entries)
                spacer_lens.push_back((int)sp.size());

            double med_spacer = median(spacer_lens);
            double ratio = med_spacer / (double)full_consensus.size();
            if (ratio < MIN_MEDIAN_SPACER_REPEAT_RATIO) {
                std::cout << "  [filtered] spacer/repeat ratio="
                          << std::defaultfloat << std::setprecision(4) << ratio
                          << " (min " << MIN_MEDIAN_SPACER_REPEAT_RATIO << "): "
                          << full_consensus.substr(0, 30) << "\n";
                continue;
            }

            consensus_arrays.push_back({full_consensus, merged_entries});

            std::vector<std::string> spacers;
            for (const auto& [rv, sp] : merged_entries)
                spacers.push_back(sp);
            crispr_arrays[full_consensus] = spacers;
        }

        // ==========================================================
        // Step 8: Output
        // ==========================================================
        // Generate timestamp string for file header
        auto now = std::chrono::system_clock::now();
        auto now_t = std::chrono::system_clock::to_time_t(now);
        std::ostringstream ts;
        ts << std::put_time(std::localtime(&now_t), "%Y-%m-%d %H:%M:%S");
        std::string timestamp_str = ts.str();

        int total_spacers_final = 0;
        for (const auto& arr : consensus_arrays)
            total_spacers_final += (int)arr.entries.size();

        int file_index = 1;
        int current_lines = 0;
        std::ofstream out;

        auto open_new_file = [&]() {
            if (out.is_open()) out.close();
            std::string filename = output_dir + "/CRISPR_Arrays_" +
                                   std::to_string(file_index) + ".txt";
            out.open(filename);
            ++file_index;
            current_lines = 0;
            // Write file header
            out << "# MCAAT — CRISPR Array Output\n";
            out << "# Generated : " << timestamp_str << "\n";
            out << "# Arrays    : " << consensus_arrays.size() << "\n";
            out << "# Spacers   : " << total_spacers_final << "\n\n";
            std::cout << "  ▸ Writing " << filename << std::endl;
        };

        open_new_file();

        int array_num = 0;
        for (const auto& arr : consensus_arrays) {
            ++array_num;
            int lines_needed = 2 + (int)arr.entries.size() + 1;
            if (current_lines + lines_needed > MAX_LINES_PER_FILE && current_lines > 0)
                open_new_file();

            out << ">Array_" << array_num << "  spacers=" << arr.entries.size() << "\n";
            out << arr.consensus << "\n";
            current_lines += 2;

            for (const auto& [repeat_var, spacer] : arr.entries) {
                const std::string& rv = (repeat_var == arr.consensus)
                    ? std::string(arr.consensus.size(), '-')
                    : repeat_var;
                out << "        " << rv << "\t" << spacer << "\n";
                ++current_lines;
            }

            out << "\n";
            ++current_lines;
        }

        if (out.is_open()) out.close();

        int total_spacers = 0;
        for (const auto& arr : consensus_arrays)
            total_spacers += (int)arr.entries.size();

        std::cout << "  ▸ Arrays: " << consensus_arrays.size()
                  << "   Spacers: " << total_spacers
                  << "   Files: " << (file_index - 1) << std::endl;
    }
};

#endif // POST_PROCESSING_H
