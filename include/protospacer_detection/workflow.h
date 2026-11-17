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
    std::filesystem::path arrays, graph, output, protospacers, extension_seeds, candidate_pool;
    std::string coverage_method = "negbin";
    double nb_dispersion = 10.;
    size_t min_contig_length = 500;
    size_t threads = 1, beam_width = 64, max_flank = 750, min_flank = 200;
    size_t max_work = 2000000, placement_work = 400000, max_placements = 128;
    size_t max_visits = 2, context_bases = 24, context_width = 256, read_length = 150, min_external = 200;
    double identity = .90, coverage = .85, repeat_identity = .85;
    double boundary_support = .85, branch_ratio = 1.5;
    double max_coverage_fold = 2.;
    bool correct_boundaries = true, save_search = false, filter_artifacts = true;
    bool source_array_crossing = true, expand_candidates = false;
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

struct FlowDecision {
    size_t position = 0;
    uint64_t node = 0;
    double incoming = 0., remaining = 0.;
    std::vector<double> outgoing;
    int chosen = -1;
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
    double coverage_baseline = 0., coverage_loss = 0.;
    size_t coverage_seen = 0, local_tips = 0, local_bubbles = 0, local_weak_links = 0;
    double seed_flow = 0., remaining_flow = 0.;
    std::vector<FlowDecision> flow_decisions;
    std::vector<uint64_t> branch_family;
    bool uses_lookahead = false;
    double selection_metric = 0., reverse_branch_loss = 0.;
};
struct LocalSegment {
    std::vector<GraphArc> arcs;
    std::string bases;
    double depth = 0.;
    size_t depth_samples = 0, terminal_indegree = 0;
    bool dead_end = false, complete = false;
};
struct SavedExtension {
    std::string sequence, origin;
    size_t begin = 0, end = 0;
};
struct CandidateEvidence {
    size_t query = 0, begin = 0, end = 0, left_source = 0, right_source = 0;
    std::string sequence, origin;
    std::array<double, 4> features{}; // branch cost, negative coding, negative array/target support
    double rank = 0.;
    bool limited = false;
};
struct CoverageBlock {
    std::vector<double> visits;
    double observed = 0., noise = 1., weight = 1.;
    size_t edges = 0;
};
struct CoverageFit {
    std::vector<double> abundance, predicted;
    double objective = 0., error = 0.;
};
struct TargetCandidate {
    size_t query = 0, placement = 0, begin = 0, end = 0, mismatches = 0;
    size_t anchor_prefix = 0, evidence_rank = size_t(-1);
    double evidence_score = -1., array_support = -1., target_support = -1.;
    std::vector<uint64_t> anchor;
    std::string target, sequence, status = "rejected", reason, target_id, extension_id;
    ExtensionWalk left, right, withheld;
    std::string withheld_sequence;
    bool withheld_left = false;
    size_t coverage_blocks = 0;
    double coverage_min = 0., coverage_max = 0.;
    std::string context_left, context_right, context_left_stop, context_right_stop;
    size_t clear_left = 0, clear_right = 0;
    double identity = 0, coverage = 0;
    bool detected = false, extended = false, limited = false;
    bool identified = false, allow_crossing = false;
    bool uses_mixture = false;
    double mixture_stability = 0., mixture_abundance = 0., mixture_error = 0.;
    size_t mixture_candidates = 0, mixture_backgrounds = 0;
    std::string mixture_origin;
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
    void filter_selected_extension(TargetCandidate &candidate) const;
    ExtensionWalk choose_walk(std::vector<ExtensionWalk> walks) const;
    bool same_route(const ExtensionWalk &a, const ExtensionWalk &b) const;
    void trim_walk(ExtensionWalk &walk, size_t length, const std::string &reason) const;
    void select_candidate(std::vector<TargetCandidate> &candidates) const;
    void process_query(size_t query);
    void assemble_candidate(TargetCandidate &candidate) const;
    void consolidate();
    void write_outputs();
    void write_summary(const std::string &status, const std::string &error = "") const;
    static double count_loss(double count, double mean, const std::string &method, double dispersion);
    static double sequence_similarity(const std::string &a, const std::string &b, double minimum);
    static bool same_graph_region(const std::string &a, const std::string &b, size_t k);
    LocalSegment probe_segment(uint64_t from, GraphArc arc, size_t query, bool backwards,
                               size_t maximum_edges, size_t &work, size_t limit) const;
    void score_seeded_walk(ExtensionWalk &walk) const;
    double seeded_walk_score(const ExtensionWalk &walk) const;
    ExtensionWalk extend_seeded(size_t query, const TargetCandidate &candidate, bool backwards, size_t &work,
                                size_t limit, std::vector<ExtensionWalk> *alternatives = nullptr,
                                int lookahead_override = -1, size_t beam_override = 0) const;
    std::vector<double> future_depths(size_t query, uint64_t from, const GraphArc &first, bool backwards,
                                      size_t distance, size_t &work, size_t limit) const;
    static void reserve_branch_families(std::vector<ExtensionWalk> &walks, size_t width);
    bool compatible_flanks(const TargetCandidate &seed, const ExtensionWalk &left,
                           const ExtensionWalk &right) const;
    static std::pair<int, double> flow_choice(double incoming, const std::vector<double> &outgoing);
    static double support_cost(const TargetCandidate &candidate);
    static bool supported_representative(const TargetCandidate &a, const TargetCandidate &b);
    ExtensionWalk extend_flow(size_t query, const TargetCandidate &candidate, bool backwards, size_t &work,
                              size_t limit) const;
    std::vector<TargetCandidate> load_extension_seeds() const;
    std::map<size_t, std::vector<SavedExtension>> load_candidate_pool(const std::vector<TargetCandidate> &seeds) const;
    TargetCandidate replay_candidate(const TargetCandidate &seed, const SavedExtension &saved, bool branch_evidence = false) const;
    void add_evidence_windows(std::vector<TargetCandidate> &candidates) const;
    std::vector<TargetCandidate> evidence_candidates(const TargetCandidate &seed, size_t &work) const;
    std::vector<TargetCandidate> mixture_candidates(const TargetCandidate &seed, size_t &work) const;
    TargetCandidate select_evidence(const TargetCandidate &seed, std::vector<TargetCandidate> candidates,
                                    std::string &audit, size_t &work, std::vector<CandidateEvidence> *evidence = nullptr) const;
    void select_evidence_regions(const std::vector<TargetCandidate> &seeds,
                                 std::vector<std::vector<CandidateEvidence>> &evidence,
                                 std::vector<TargetCandidate> &predictions) const;
    static double coding_evidence(const std::string &sequence);
    static CoverageFit fit_coverage(const std::vector<CoverageBlock> &blocks, size_t columns,
                                   const std::vector<double> &observations);
    TargetCandidate select_mixture(const TargetCandidate &seed, std::vector<TargetCandidate> candidates,
                                   std::string &audit, size_t &work) const;
    int run_seeded();
};
} // namespace mcaat::proto
#endif
