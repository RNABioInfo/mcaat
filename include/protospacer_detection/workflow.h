#ifndef MCAAT_PROTOSPACER_WORKFLOW_H
#define MCAAT_PROTOSPACER_WORKFLOW_H

#include "sdbg/sdbg.h"
#include <array>
#include <csignal>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace mcaat::proto {

struct WorkflowOptions {
    std::filesystem::path arrays, graph, output, protospacers;
    size_t threads = 1, beam_width = 64, max_flank = 750, min_flank = 200;
    size_t max_work = 2000000, placement_work = 400000, max_placements = 128;
    size_t max_visits = 2, context_bases = 24, context_width = 256, read_length = 150, min_external = 200;
    double identity = .90, coverage = .85, repeat_identity = .85;
    double boundary_support = .85, branch_ratio = 1.5;
    bool correct_boundaries = true, save_search = false, filter_artifacts = true;
    bool source_array_crossing = true;
};

struct SpacerSource {
    std::string id, array, original_spacer, original_repeat, spacer, repeat;
    std::string source_left, source_right;
    size_t trim_left = 0, trim_right = 0;
};

struct GraphArc {
    uint64_t node = 0, edge = 0;
    double count = 0;
};

struct ExtensionWalk {
    std::vector<uint64_t> nodes;
    std::string bases, stop = "dead_end";
    std::vector<double> counts, capacities, all_counts;
    std::vector<size_t> sample_positions;
    size_t profile_samples = 0;
    double deviation = 0;
    double transition_energy = 0, branch_loss = 0;
    size_t transitions = 0, branch_observations = 0;
    std::vector<double> branch_costs;
    size_t clear_at = 0, work = 0;
    size_t source_until = 0, source_last = 0, anchor_padding = 0;
    double flow = 0;
    bool external = false, limited = false;
    bool crosses_source = false, ambiguous = false, pre_context_fork = false;
    size_t competing_routes = 0, first_ambiguous_base = 0, first_context_fork = 0;
    double runner_up_score = -1.;
};

struct TargetCandidate {
    size_t query = 0, placement = 0, begin = 0, end = 0, mismatches = 0;
    size_t anchor_prefix = 0;
    std::vector<uint64_t> anchor;
    std::string target, sequence, status = "rejected", reason, target_id, extension_id;
    ExtensionWalk left, right;
    std::string context_left, context_right, context_left_stop, context_right_stop;
    size_t clear_left = 0, clear_right = 0;
    double identity = 0, coverage = 0;
    bool detected = false, extended = false, limited = false;
    bool identified = false, allow_crossing = false;
};

// One owner for graph access, detection, extension, filtering and output.
// All graph operations are read-only. Per-query state is local to its worker.
class ProtospacerWorkflow {
  public:
    static int main(int argc, char **argv);
    explicit ProtospacerWorkflow(WorkflowOptions options);
    int run();

  private:
    WorkflowOptions options_;
    std::unique_ptr<SDBG> graph_;
    std::vector<SpacerSource> sources_;
    std::vector<std::pair<std::string, std::string>> supplied_;
    std::map<std::string, std::set<uint64_t>> array_nodes_, array_repeat_nodes_;
    std::vector<std::set<uint64_t>> query_source_nodes_, query_repeat_nodes_;
    std::vector<std::pair<std::string, std::vector<size_t>>> queries_;
    std::vector<std::vector<TargetCandidate>> results_;
    std::vector<std::string> query_status_;
    std::vector<size_t> query_work_, mapping_work_;
    std::vector<uint8_t> query_limited_;
    std::vector<std::array<double, 3>> query_timings_;
    std::map<std::string, double> timings_;
    static volatile std::sig_atomic_t interrupted_;
    static void interrupt(int signal);

    static WorkflowOptions parse_options(int argc, char **argv);
    static void print_help();
    static void self_test();
    static std::string reverse_complement(std::string sequence);
    static std::string canonical(const std::string &sequence);
    static std::string json(const std::string &value);
    static std::string numbered(const std::string &prefix, size_t value);
    static bool dna(const std::string &sequence);
    static double triplet_entropy(const std::string &sequence);
    static bool informative(const std::string &sequence);
    static size_t adapter_position(const std::string &sequence, bool backwards = false);
    static int repeat_context(const std::string &outward, const std::string &remainder,
                              const std::string &repeat, double identity);

    std::string label(uint64_t node) const;
    uint64_t find_node(const std::string &word) const;
    std::vector<GraphArc> neighbors(uint64_t node, bool backwards) const;
    char added_base(const GraphArc &arc, bool backwards) const;
    std::string spell(const std::vector<uint64_t> &path) const;
    std::vector<std::pair<uint64_t, size_t>> find_starts(const std::string &word, size_t errors, size_t &work,
                                                         bool &limited) const;
    void search_range(const std::string &word, size_t errors, uint64_t lo, uint64_t hi, size_t depth,
                      size_t used, size_t &work, bool &limited, std::map<uint64_t, size_t> &found) const;

    void verify_index() const;
    void load_sources();
    void prepare_source_nodes();
    void load_supplied_targets();
    void correct_boundaries();
    void group_queries();
    std::vector<TargetCandidate> map_spacer(size_t query, size_t &work, bool &limited) const;
    void update_profile(ExtensionWalk &walk, bool finish = false) const;
    double walk_score(const ExtensionWalk &walk) const;
    void rank_walks(std::vector<ExtensionWalk> &walks, size_t limit, bool diverse) const;
    bool retain_completed(std::vector<ExtensionWalk> &walks) const;
    void save_walks(const TargetCandidate &candidate, bool backwards, const std::vector<ExtensionWalk> &walks,
                    size_t &serial) const;
    int source_context(size_t query, const TargetCandidate &candidate, const std::string &outward,
                       bool backwards) const;
    size_t source_clearance(size_t query, const TargetCandidate &candidate, bool backwards) const;
    ExtensionWalk extend_flank(size_t query, const TargetCandidate &candidate, bool backwards, size_t &work,
                               size_t work_limit) const;
    bool advance_walk(ExtensionWalk &walk, const GraphArc &arc, size_t query,
                      const TargetCandidate &candidate, bool backwards, double branch_cost = 0.) const;
    bool extension_barrier(ExtensionWalk &walk, const TargetCandidate &candidate, bool backwards) const;
    ExtensionWalk choose_walk(std::vector<ExtensionWalk> walks) const;
    bool same_route(const ExtensionWalk &a, const ExtensionWalk &b) const;
    void trim_walk(ExtensionWalk &walk, size_t length, const std::string &reason) const;
    void select_candidate(std::vector<TargetCandidate> &candidates) const;
    void process_query(size_t query);
    void assemble_candidate(TargetCandidate &candidate) const;
    void consolidate();
    void write_outputs();
    void write_summary(const std::string &status, const std::string &error = "") const;
};
} // namespace mcaat::proto
#endif
