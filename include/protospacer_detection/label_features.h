// Training data for a learned branch decision, built from reference-confirmed
// protospacer regions.
//
// A label is one place in a reference catalogue where a spacer detected in a
// sample matches at full length with no mismatch, together with 1250 bases of
// reference sequence on each side. The reference region therefore spells the
// path the extension should have taken through that sample's graph.
//
// Walking that path outward from the spacer, every node with more than one
// neighbour is a decision the beam search has to make. The neighbour the
// reference continues into is a positive; its siblings are negatives. The
// reference decides which is which and is then discarded: every feature
// recorded for a neighbour comes from graph topology and edge multiplicity
// alone, so a model fitted on them consults no reference when it runs.
//
//   LabelReader       reads labels_by_dataset.json and hands out one dataset's
//                     entries.
//   TopologyFeatures  maps an entry's region into the graph and collects the
//                     features of every neighbour at every decision on it.
//   FeatureWriter     writes those rows as a table, with a summary beside it.
#ifndef MCAAT_PROTOSPACER_LABEL_FEATURES_H
#define MCAAT_PROTOSPACER_LABEL_FEATURES_H

#include "sdbg/sdbg.h"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace mcaat::proto {

// One reference-confirmed protospacer region.
struct LabelEntry {
    std::string id;          // label identifier, e.g. L00001
    std::string spacer_key;  // spacer identifier, e.g. S000116
    std::string left;        // reference sequence before the spacer
    std::string spacer;      // the spacer, as detected in the sample
    std::string right;       // reference sequence after the spacer
    std::string record;      // reference record the region was taken from
    std::string strand;      // plus or minus, relative to that record
    std::vector<std::string> arrays; // array:index positions in this sample
    size_t region_begin = 0, region_end = 0; // coordinates in the record

    std::string region() const { return left + spacer + right; }
    size_t spacer_begin() const { return left.size(); }
    size_t spacer_end() const { return left.size() + spacer.size(); }
};

// Class 1. Reads the label file and serves the entries of one dataset.
class LabelReader {
  public:
    explicit LabelReader(const std::filesystem::path &file);

    bool has(const std::string &dataset) const;
    const std::vector<LabelEntry> &entries(const std::string &dataset) const;
    std::vector<std::string> datasets() const;
    size_t size() const;  // entries over all datasets
    const std::map<std::string, std::string> &meta() const { return meta_; }

  private:
    std::map<std::string, std::vector<LabelEntry>> by_dataset_;
    std::map<std::string, std::string> meta_;
    static const std::vector<LabelEntry> empty_;
};

// One neighbour at one decision, with its label and its features. Every field
// below the identity block is read from the graph; none of them sees a base.
struct BranchFeature {
    // identity, not for fitting
    std::string dataset, label, spacer_key;
    size_t position = 0;     // bases walked from the spacer when the decision is met
    bool backwards = false;  // true walking left of the spacer, false walking right
    uint64_t node = 0, edge = 0;
    int target = 0;          // 1 the reference continues here, 0 it does not

    // the decision
    size_t siblings = 0;
    double branch_entropy = 0;      // over the siblings' multiplicity shares
    double parent_multiplicity = 0; // multiplicity the walk arrived on
    double flow_residual = 0;       // parent minus the sum of the siblings
    double flow_residual_share = 0;

    // this neighbour
    double multiplicity = 0;
    double multiplicity_share = 0;         // against the sum of the siblings
    double ratio_best_sibling = 0;         // against the largest sibling
    double ratio_running_median = 0;       // against the walk's running median
    double ratio_origin = 0;               // against the multiplicity at the spacer boundary
    int in_window = 0;                     // inside [0.1 m, 5 m] of the running median
    size_t out_degree = 0, in_degree = 0;

    // what lies ahead of this neighbour, bounded by the lookahead
    size_t run_ahead = 0;        // bases before the next decision
    int run_ends_tip = 0;        // that run stops with nothing to continue into
    size_t forks_ahead = 0;
    size_t reachable_nodes = 0;
    double run_depth_mean = 0, run_depth_min = 0, run_depth_cv = 0;
    int rejoins_sibling = 0;     // meets a sibling again, so the pair is a bubble
    size_t bubble_length = 0;

    // the walk so far
    size_t walked = 0, forks_passed = 0;
    double running_median = 0, running_cv = 0;
};

// Class 2. Finds the nodes of a region and collects the features.
class TopologyFeatures {
  public:
    struct Options {
        size_t lookahead = 300;    // bases probed ahead of a neighbour
        size_t probe_nodes = 2000; // nodes a probe may visit
        double window_low = .1, window_high = 5.; // the beam search admission window
        bool require_both_sides = false; // drop a region that maps on one side only
    };
    struct Stats {
        size_t entries = 0, mapped = 0, unmapped = 0;
        size_t anchors = 0, reanchors = 0, steps = 0;
        size_t decisions = 0, positives = 0, negatives = 0;
        size_t decisions_skipped = 0; // the reference base had no neighbour
        size_t sides_mapped = 0, sides_empty = 0;
    };

    explicit TopologyFeatures(const SDBG &graph);
    TopologyFeatures(const SDBG &graph, Options options);

    // Every neighbour at every decision on this entry's reference path.
    std::vector<BranchFeature> collect(const std::string &dataset, const LabelEntry &entry);
    const Stats &stats() const { return stats_; }
    size_t k() const { return k_; }

  private:
    struct Probe {
        size_t run = 0, forks = 0, nodes = 0;
        bool tip = false;
        double mean = 0, minimum = 0, cv = 0;
        std::vector<uint64_t> seen;
    };

    const SDBG &graph_;
    Options options_;
    size_t k_ = 0;
    Stats stats_{};

    // Graph access, kept identical to ProtospacerWorkflow so that a feature
    // measured here is the feature the extension sees.
    std::string label(uint64_t node) const;
    uint64_t find_node(const std::string &word) const;
    struct Arc {
        uint64_t node = 0, edge = 0;
        double count = 0;
    };
    std::vector<Arc> neighbors(uint64_t node, bool backwards) const;
    char added_base(const Arc &arc, bool backwards) const;

    Probe probe(uint64_t from, const Arc &arc, bool backwards) const;
    void walk(const std::string &dataset, const LabelEntry &entry, const std::string &region,
              size_t start, bool backwards, std::vector<BranchFeature> &out);
    static double median(std::vector<double> values);
    static double variation(const std::vector<double> &values);
};

// Class 3. Writes the rows out.
class FeatureWriter {
  public:
    explicit FeatureWriter(const std::filesystem::path &file);
    ~FeatureWriter();

    void write(const std::vector<BranchFeature> &rows);
    void write(const BranchFeature &row);
    // A summary beside the table: <file>.summary.json
    void finish(const TopologyFeatures::Stats &stats, const std::string &graph,
                const std::string &labels, const std::vector<std::string> &datasets);
    size_t rows() const { return rows_; }
    static std::vector<std::string> columns();

  private:
    std::filesystem::path file_;
    std::ofstream out_;
    size_t rows_ = 0;
};

} // namespace mcaat::proto
#endif
