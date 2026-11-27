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

// One neighbour at one decision, with its label and its features.
//
// Nine features, and the reason each one is here. Two conditions decide the
// set. A feature has to differ between the candidates at a single decision,
// because a value they all share cannot say which of them to take. And it has
// to be a quantity whose relevance is plain, not a number that happens to
// correlate.
//
// Multiplicity enters only as a ratio. An absolute multiplicity means
// different things in a metagenome sequenced to one depth than in another, so
// a model fitted on raw counts does not carry across samples. The raw count is
// kept beside the identity fields for inspection, not for fitting.
struct BranchFeature {
    // identity and diagnostics, not for fitting
    std::string dataset, label, spacer_key;
    bool backwards = false;  // true walking left of the spacer, false walking right
    size_t walked = 0;       // bases walked from the spacer when the decision is met
    uint64_t node = 0, edge = 0;
    double multiplicity = 0; // the raw count behind the ratios below
    int target = 0;          // 1 the reference continues here, 0 it does not

    // This candidate's multiplicity over the sum at the decision. The direct
    // statement that this branch carries the coverage and its rivals do not.
    double multiplicity_share = 0;
    // Bases that can be walked from this candidate before the next decision.
    // A long unbranching stretch is a real contig; a branch that forks again
    // at once is inside a tangle.
    size_t run_ahead = 0;
    // That stretch ends with nothing to continue into. A dead end is not a
    // genome.
    int run_ends_tip = 0;
    // The lowest multiplicity on that stretch, over the running median. The
    // bottleneck, which is where a chimeric join shows itself and which is
    // already one of the seven scoring criteria.
    double run_depth_min_ratio = 0;
    // How much the multiplicity varies along that stretch. Coverage holding
    // level is what locally_smooth asks of a whole reconstruction, measured
    // here for one branch.
    double run_depth_cv = 0;
    // Predecessors of the candidate against the direction of travel. A node
    // several paths enter is a repeat entrance, which is where a walk crosses
    // from one genome into another.
    size_t in_degree = 0;
};

// Class 2. Finds the nodes of a region and collects the features.
class TopologyFeatures {
  public:
    struct Options {
        size_t lookahead = 300;    // bases probed ahead of a neighbour
        size_t probe_nodes = 2000; // nodes a probe may visit
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
        size_t run = 0;
        bool tip = false;
        double minimum = 0, cv = 0;
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
