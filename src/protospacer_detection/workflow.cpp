#include "protospacer_detection/workflow.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <zlib.h>

namespace mcaat::proto {
namespace fs = std::filesystem;
static constexpr uint64_t missing = std::numeric_limits<uint64_t>::max();

volatile std::sig_atomic_t ProtospacerWorkflow::interrupted_ = 0;
void ProtospacerWorkflow::interrupt(int signal) { interrupted_ = signal; }

ProtospacerWorkflow::ProtospacerWorkflow(WorkflowOptions options) : options_(std::move(options)) {}

void ProtospacerWorkflow::print_help() {
    std::cout << "MCAAT protospacer workflow v9\n"
                 "mcaat phage-curate --arrays DIR --graph PREFIX --output "
                 "NEW_DIR [options]\n"
                 "Final protospacers must have an accepted extension; all detected "
                 "candidates are saved separately.\n"
                 "One workflow: spacer alignment -> source-context check -> "
                 "extension -> filtering -> output.\n"
                 "--threads 1 --beam-width 64 --max-flank-length 750 "
                 "--min-extension-flank 200\n"
                 "--min-identity 0.90 --min-query-coverage 0.85 "
                 "--array-repeat-identity 0.85\n"
                 "--max-query-states 2000000 --max-placement-states 400000 "
                 "--max-placements 128\n"
                 "--branch-ratio 1.5 --max-node-visits 2\n"
                 "--max-context-paths 256 --read-length 150 --min-external-bases 200\n"
                 "--max-coverage-fold 2: withhold the selected extension if read-span "
                 "medians exceed this ratio; 0 disables.\n"
                 "--correct-spacer-boundaries true --boundary-support 0.85\n"
                 "--filter-extension-artifacts true: stop at adapter cores or "
                 "low-complexity windows.\n"
                 "--save-search false: optionally stream all finished/pruned "
                 "flank paths to gzip TSVs.\n"
                 "The flank limit is a maximum. One best actual route is saved; "
                 "close alternatives and search limits are reported.\n"
                 "Source-array traversal never establishes a target. Graph files "
                 "are never modified.\n"
                 "--protospacers FASTA: extend supplied targets; preserve every "
                 "supplied input separately.\n"
                 "Supplied targets require an exact full-sequence graph embedding.\n"
                 "--source-array-crossing true: one passage for identified "
                 "targets only, never discovery.\n"
                 "--extension-seeds DIR: extend ALL saved candidate starts from this "
                 "workflow directory.\n"
                 "--coverage-method poisson|negbin|huber|flow|lookahead|all --nb-dispersion 10 "
                 "--min-contig-length 500\n"
                 "Seeded mode uses local MEGAHIT-derived graph rules and keeps up to "
                 "max-flank-length on EACH side.\n"
                 "No reference-guided path selection.\n"
                 "Output must be new or empty. All counts and search limits are "
                 "recorded in run_summary.json (seeded mode: summary.json).\n"
                 "--self-test runs arithmetic/context checks without a graph or "
                 "biological dataset.\n";
}
WorkflowOptions ProtospacerWorkflow::parse_options(int argc, char **argv) {
    WorkflowOptions o;
    bool seeded_options = false;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (++i >= argc)
            throw std::invalid_argument("missing value for " + flag);
        std::string value = argv[i];
        seeded_options |= flag == "--coverage-method" || flag == "--nb-dispersion" || flag == "--min-contig-length";
        if (flag == "--arrays")
            o.arrays = value;
        else if (flag == "--protospacers")
            o.protospacers = value;
        else if (flag == "--graph")
            o.graph = value;
        else if (flag == "--output")
            o.output = value;
        else if (flag == "--extension-seeds")
            o.extension_seeds = value;
        else if (flag == "--coverage-method") {
            if (value != "poisson" && value != "negbin" && value != "huber" && value != "flow" && value != "lookahead" && value != "all")
                throw std::invalid_argument("coverage method must be poisson, negbin, huber, flow, lookahead or all");
            o.coverage_method = value;
        } else if (flag == "--correct-spacer-boundaries" || flag == "--save-search" ||
                   flag == "--filter-extension-artifacts" || flag == "--source-array-crossing") {
            if (value != "true" && value != "false")
                throw std::invalid_argument(flag + " requires true or false");
            if (flag == "--save-search")
                o.save_search = value == "true";
            else if (flag == "--source-array-crossing")
                o.source_array_crossing = value == "true";
            else if (flag == "--filter-extension-artifacts")
                o.filter_artifacts = value == "true";
            else
                o.correct_boundaries = value == "true";
        } else if (flag == "--threads" || flag == "--beam-width" || flag == "--max-flank-length" ||
                   flag == "--min-extension-flank" || flag == "--max-query-states" ||
                   flag == "--max-placement-states" || flag == "--max-placements" ||
                   flag == "--max-node-visits" || flag == "--max-context-paths" || flag == "--read-length" ||
                   flag == "--min-external-bases" || flag == "--min-contig-length") {
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument(flag + " requires a positive integer");
            size_t n = std::stoull(value);
            if (!n || n > 10000000)
                throw std::invalid_argument(flag + " is out of range");
            if (flag == "--threads")
                o.threads = n;
            if (flag == "--beam-width")
                o.beam_width = n;
            if (flag == "--max-flank-length")
                o.max_flank = n;
            if (flag == "--min-extension-flank")
                o.min_flank = n;
            if (flag == "--max-query-states")
                o.max_work = n;
            if (flag == "--max-placement-states")
                o.placement_work = n;
            if (flag == "--max-placements")
                o.max_placements = n;
            if (flag == "--max-node-visits")
                o.max_visits = n;
            if (flag == "--max-context-paths")
                o.context_width = n;
            if (flag == "--read-length")
                o.read_length = n;
            if (flag == "--min-external-bases")
                o.min_external = n;
            if (flag == "--min-contig-length")
                o.min_contig_length = n;
        } else if (flag == "--min-identity" || flag == "--min-query-coverage" ||
                   flag == "--array-repeat-identity" || flag == "--boundary-support" ||
                   flag == "--branch-ratio" || flag == "--max-coverage-fold" || flag == "--nb-dispersion") {
            size_t used = 0;
            double n = std::stod(value, &used);
            if (used != value.size() || !std::isfinite(n))
                throw std::invalid_argument("invalid " + flag);
            if (flag == "--nb-dispersion") {
                if (n <= 0.)
                    throw std::invalid_argument("negative-binomial dispersion must be positive");
                o.nb_dispersion = n;
            } else if (flag == "--max-coverage-fold") {
                if (n != 0. && (n <= 1. || n > 1000000.))
                    throw std::invalid_argument("coverage fold must exceed one, or be zero to disable");
                o.max_coverage_fold = n;
            } else if (flag == "--branch-ratio") {
                if (n <= 1 || n > 1000000)
                    throw std::invalid_argument("branch ratio must exceed one");
                o.branch_ratio = n;
            } else {
                if (n <= .5 || n > 1)
                    throw std::invalid_argument(flag + " must be in (0.5, 1]");
                if (flag == "--min-identity")
                    o.identity = n;
                if (flag == "--min-query-coverage")
                    o.coverage = n;
                if (flag == "--array-repeat-identity")
                    o.repeat_identity = n;
                if (flag == "--boundary-support")
                    o.boundary_support = n;
            }
        } else
            throw std::invalid_argument("unsupported workflow option: " + flag + "; see phage-curate --help");
    }
    if (o.arrays.empty() || o.graph.empty() || o.output.empty())
        throw std::invalid_argument("--arrays, --graph and --output are required");
    if (seeded_options && o.extension_seeds.empty())
        throw std::invalid_argument("coverage-model options require --extension-seeds DIR; see phage-curate --help");
    if (!fs::is_directory(o.arrays))
        throw std::invalid_argument("array directory does not exist");
    if (fs::is_directory(o.graph))
        o.graph /= "graph";
    if (!fs::is_regular_file(o.graph.string() + ".sdbg_info"))
        throw std::invalid_argument("graph metadata does not exist");
    if (o.threads > 1024 || o.beam_width > 256 || o.max_placements > 256 || o.max_visits > 4 ||
        o.max_flank > 10000 || o.min_external > o.max_flank || o.context_width > 4096 || o.read_length < 25 ||
        o.read_length > 10000 || o.min_flank > o.max_flank || o.placement_work > o.max_work)
        throw std::invalid_argument("inconsistent or excessive workflow limits");
    if (fs::exists(o.output) && (!fs::is_directory(o.output) || !fs::is_empty(o.output)))
        throw std::invalid_argument("output must be new or empty; existing results are preserved");
    o.arrays = fs::absolute(o.arrays);
    o.graph = fs::absolute(o.graph);
    o.output = fs::absolute(o.output);
    if (!o.extension_seeds.empty()) {
        if (!fs::is_directory(o.extension_seeds) || !o.protospacers.empty())
            throw std::invalid_argument("extension seeds require a workflow directory and cannot be combined "
                                        "with --protospacers");
        o.extension_seeds = fs::absolute(o.extension_seeds);
    }
    if (!o.protospacers.empty()) {
        if (!fs::is_regular_file(o.protospacers))
            throw std::invalid_argument("supplied protospacer FASTA does not exist");
        o.protospacers = fs::absolute(o.protospacers);
        // Already identified graph intervals must not be substituted or clipped.
        o.identity = o.coverage = 1.;
    }
    return o;
}
int ProtospacerWorkflow::main(int argc, char **argv) {
    if (argc == 1 || (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))) {
        print_help();
        return argc == 1 ? 1 : 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        self_test();
        return 0;
    }
    ProtospacerWorkflow workflow(parse_options(argc, argv));
    fs::create_directories(workflow.options_.output);
    std::signal(SIGINT, interrupt);
    std::signal(SIGTERM, interrupt);
    try {
        return workflow.run();
    } catch (const std::exception &e) {
        workflow.write_summary("failed", e.what());
        throw;
    }
}

std::string ProtospacerWorkflow::reverse_complement(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char &c : s) {
        switch (c) {
        case 'A':
            c = 'T';
            break;
        case 'C':
            c = 'G';
            break;
        case 'G':
            c = 'C';
            break;
        case 'T':
            c = 'A';
            break;
        default:
            throw std::invalid_argument("non-ACGT sequence");
        }
    }
    return s;
}

std::string ProtospacerWorkflow::canonical(const std::string &s) {
    return std::min(s, reverse_complement(s));
}

bool ProtospacerWorkflow::dna(const std::string &s) {
    return !s.empty() && s.find_first_not_of("ACGT") == std::string::npos;
}

std::string ProtospacerWorkflow::json(const std::string &s) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\')
            out << '\\' << c;
        else if (c < 32)
            out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
        else
            out << c;
    }
    out << '"';
    return out.str();
}

std::string ProtospacerWorkflow::numbered(const std::string &prefix, size_t value) {
    std::ostringstream s;
    s << prefix << std::setw(6) << std::setfill('0') << value;
    return s.str();
}

double ProtospacerWorkflow::triplet_entropy(const std::string &s) {
    if (s.size() < 3)
        return 0;
    std::array<size_t, 64> counts{};
    size_t word = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        word = ((word << 2) | std::string("ACGT").find(s[i])) & 63;
        if (i >= 2)
            ++counts[word];
    }
    double entropy = 0;
    for (auto n : counts)
        if (n) {
            double p = double(n) / (s.size() - 2);
            entropy -= p * std::log2(p);
        }
    return entropy;
}

bool ProtospacerWorkflow::informative(const std::string &sequence) {
    if (triplet_entropy(sequence) < 2.)
        return false;
    size_t run = 0, longest = 0;
    char previous = 0;
    for (char c : sequence) {
        run = c == previous ? run + 1 : 1;
        previous = c;
        longest = std::max(longest, run);
    }
    return longest * 2 < sequence.size();
}

size_t ProtospacerWorkflow::adapter_position(const std::string &sequence, bool backwards) {
    // The exact shared TruSeq prefix catches library-specific/graph variants
    // beyond that prefix. Full cores additionally tolerate one substitution.
    // Never relax the short prefix: its specificity depends on an exact match.
    static const std::vector<std::string> motifs = [] {
        std::vector<std::string> m{"AGATCGGAAGAGC", "AGATCGGAAGAGCACACGT", "AGATCGGAAGAGCGTCGT"};
        m.reserve(12);
        for (size_t i = 0; i < 3; ++i)
            m.push_back(reverse_complement(m[i]));
        for (size_t i = 0; i < 6; ++i)
            m.emplace_back(m[i].rbegin(), m[i].rend());
        return m;
    }();
    size_t found = std::string::npos;
    for (size_t i = backwards ? 6 : 0; i < (backwards ? 12 : 6); ++i) {
        const auto &motif = motifs[i];
        const size_t allowed = motif.size() == 13 ? 0 : 1;
        for (size_t position = 0; position + motif.size() <= sequence.size() && position < found; ++position) {
            size_t changes = 0;
            for (size_t j = 0; j < motif.size() && changes <= allowed; ++j)
                changes += sequence[position + j] != motif[j];
            if (changes <= allowed) found = position;
        }
    }
    return found;
}

int ProtospacerWorkflow::repeat_context(const std::string &outward, const std::string &remainder,
                                        const std::string &repeat, double identity) {
    // 0=pending, 1=source repeat, 2=incompatible with this source context.
    if (outward.size() < remainder.size())
        return 0;
    const size_t edits = size_t(std::floor((1 - identity) * repeat.size() + 1e-9));
    const size_t available = outward.size() - remainder.size();
    std::vector<size_t> row(repeat.size() + 1), next(row.size());
    std::iota(row.begin(), row.end(), 0);
    for (size_t i = 1; i <= std::min(available, repeat.size() + edits); ++i) {
        next[0] = i;
        for (size_t j = 1; j < row.size(); ++j)
            next[j] = std::min({row[j] + 1, next[j - 1] + 1,
                                row[j - 1] + (outward[remainder.size() + i - 1] != repeat[j - 1])});
        row.swap(next);
        if (i + edits >= repeat.size() && row.back() <= edits)
            return 1;
    }
    return available < repeat.size() + edits && *std::min_element(row.begin(), row.end()) <= edits ? 0 : 2;
}

std::string ProtospacerWorkflow::label(uint64_t node) const {
    if (node >= graph_->size())
        return {};
    std::vector<uint8_t> encoded(graph_->k());
    graph_->GetLabel(node, encoded.data());
    std::string s;
    for (auto c : encoded) {
        if (c < 1 || c > 4)
            return {};
        s += "ACGT"[c - 1];
    }
    return s;
}

uint64_t ProtospacerWorkflow::find_node(const std::string &word) const {
    if (word.size() != size_t(graph_->k()) || !dna(word))
        return missing;
    std::vector<uint8_t> encoded;
    for (char c : word)
        encoded.push_back(std::string("ACGT").find(c) + 1);
    uint64_t n = graph_->IndexBinarySearch(encoded.data());
    return n != missing && label(n) == word ? n : missing;
}

std::vector<GraphArc> ProtospacerWorkflow::neighbors(uint64_t node, bool backwards) const {
    std::vector<GraphArc> arcs;
    if (node >= graph_->size())
        return arcs;
    auto add = [&](uint64_t edge, uint64_t to) {
        auto c = graph_->GetW(edge);
        if (c > 4)
            c -= 4;
        if (c >= 1 && c <= 4 && graph_->IsValidEdge(edge))
            arcs.push_back({to, edge, double(graph_->EdgeMultiplicity(edge))});
    };
    if (!backwards) {
        uint64_t e = node;
        do {
            auto c = graph_->GetW(e);
            if (c > 4)
                c -= 4;
            if (c >= 1 && c <= 4 && graph_->IsValidEdge(e))
                add(e, graph_->Forward(e));
            --e;
        } while (e != missing && !graph_->IsLastOrTip(e));
    } else {
        uint64_t first = graph_->Backward(node);
        if (first >= graph_->size())
            return arcs;
        auto c = graph_->GetW(first);
        unsigned ends = graph_->IsLastOrTip(first);
        add(first, graph_->GetLastIndex(first));
        for (uint64_t e = first + 1; ends < 5 && e < graph_->size(); ++e) {
            ends += graph_->IsLastOrTip(e);
            auto w = graph_->GetW(e);
            if (w == c)
                break;
            if (w == c + 4)
                add(e, graph_->GetLastIndex(e));
        }
    }
    std::sort(arcs.begin(), arcs.end(), [](const auto &a, const auto &b) { return a.node < b.node; });
    return arcs;
}

char ProtospacerWorkflow::added_base(const GraphArc &arc, bool backwards) const {
    if (backwards) {
        auto s = label(arc.node);
        if (s.empty())
            throw std::runtime_error("invalid graph label");
        return s.front();
    }
    auto c = graph_->GetW(arc.edge);
    if (c > 4)
        c -= 4;
    if (c < 1 || c > 4)
        throw std::runtime_error("invalid graph base");
    return "ACGT"[c - 1];
}

std::string ProtospacerWorkflow::spell(const std::vector<uint64_t> &path) const {
    if (path.empty())
        return {};
    std::string s = label(path.front());
    for (size_t i = 1; i < path.size(); ++i) {
        bool connected = false;
        for (const auto &a : neighbors(path[i - 1], false))
            if (a.node == path[i]) {
                s += added_base(a, false);
                connected = true;
                break;
            }
        if (!connected)
            throw std::runtime_error("disconnected graph path");
    }
    return s;
}
void ProtospacerWorkflow::search_range(const std::string &word, size_t errors, uint64_t lo, uint64_t hi,
                                       size_t depth, size_t used, size_t &work, bool &limited,
                                       std::map<uint64_t, size_t> &found) const {
    if (lo >= hi || limited || interrupted_)
        return;
    if (work++ >= options_.placement_work) {
        limited = true;
        return;
    }
    if (used == errors && depth < word.size() && hi - lo > 8) {
        // Once the mismatch budget is spent, locate the remaining exact suffix
        // by two binary searches instead of branching character by character.
        auto exact_bound = [&](bool upper) {
            uint64_t a = lo, b = hi;
            std::vector<uint8_t> encoded(graph_->k());
            while (a < b) {
                if (work++ >= options_.placement_work) {
                    limited = true;
                    return a;
                }
                uint64_t middle = a + (b - a) / 2;
                graph_->GetLabel(middle, encoded.data());
                int comparison = 0;
                for (int j = int(word.size() - depth) - 1; j >= 0; --j) {
                    unsigned actual = encoded[graph_->k() - word.size() + j];
                    unsigned expected = std::string("ACGT").find(word[j]) + 1;
                    if (actual != expected) {
                        comparison = actual < expected ? -1 : 1;
                        break;
                    }
                }
                if (comparison < 0 || (upper && comparison == 0))
                    a = middle + 1;
                else
                    b = middle;
            }
            return a;
        };
        auto first = exact_bound(false), last = exact_bound(true);
        search_range(word, errors, first, last, word.size(), used, work, limited, found);
        return;
    }
    if (depth == word.size() || hi - lo <= 8) {
        for (uint64_t index = lo; index < hi;) {
            if (work++ >= options_.placement_work || found.size() >= 4096) {
                limited = true;
                return;
            }
            uint64_t node = graph_->GetLastIndex(index);
            std::string sequence = label(node);
            if (sequence.size() >= word.size()) {
                size_t errors_here = 0, offset = sequence.size() - word.size();
                for (size_t j = 0; j < word.size(); ++j)
                    errors_here += sequence[offset + j] != word[j];
                if (errors_here <= errors)
                    found.emplace(node, errors_here);
            }
            if (node < index)
                throw std::runtime_error("invalid succinct-index interval");
            index = node + 1;
        }
        return;
    }
    // The index is colexicographic: suffixes form contiguous ranges. Explore
    // the matching character first, then mismatches; limits remain explicit.
    size_t position = graph_->k() - 1 - depth;
    auto bound = [&](unsigned base) {
        uint64_t a = lo, b = hi;
        std::vector<uint8_t> encoded(graph_->k());
        while (a < b) {
            if (work++ >= options_.placement_work) {
                limited = true;
                return a;
            }
            auto middle = a + (b - a) / 2;
            graph_->GetLabel(middle, encoded.data());
            if (encoded[position] < base)
                a = middle + 1;
            else
                b = middle;
        }
        return a;
    };
    unsigned expected = std::string("ACGT").find(word[word.size() - 1 - depth]) + 1;
    for (unsigned turn = 0; turn < 5 && !limited; ++turn) {
        unsigned base = turn ? turn : expected;
        if (turn && base == expected)
            continue;
        size_t changes = used + (base != expected);
        if (changes > errors)
            continue;
        uint64_t first = bound(base), last = bound(base + 1);
        search_range(word, errors, first, last, depth + 1, changes, work, limited, found);
    }
}

std::vector<std::pair<uint64_t, size_t>>
ProtospacerWorkflow::find_starts(const std::string &word, size_t errors, size_t &work, bool &limited) const {
    std::map<uint64_t, size_t> found;
    if (word.size() > size_t(graph_->k()))
        throw std::runtime_error("seed exceeds graph k");
    search_range(word, errors, 0, graph_->size(), 0, 0, work, limited, found);
    std::vector<std::pair<uint64_t, size_t>> result(found.begin(), found.end());
    std::sort(result.begin(), result.end(), [](const auto &a, const auto &b) {
        return std::tie(a.second, a.first) < std::tie(b.second, b.first);
    });
    return result;
}

void ProtospacerWorkflow::verify_index() const {
    size_t checked = 0;
    for (size_t trial = 1; trial <= 64; ++trial) {
        uint64_t node = graph_->GetLastIndex((graph_->size() / 65) * trial);
        auto sequence = label(node);
        if (sequence.empty())
            continue;
        if (find_node(sequence) != node)
            throw std::runtime_error("exact index does not recover graph label");
        for (size_t errors : {size_t(0), size_t(1), size_t(2)}) {
            auto seed = sequence.substr(sequence.size() - std::min(size_t(16), sequence.size()));
            if (errors)
                seed[seed.size() / 2] = seed[seed.size() / 2] == 'A' ? 'C' : 'A';
            if (errors > 1)
                seed[seed.size() / 3] = seed[seed.size() / 3] == 'A' ? 'C' : 'A';
            size_t work = 0;
            bool limited = false;
            auto starts = find_starts(seed, errors, work, limited);
            if (limited)
                continue;
            if (std::none_of(starts.begin(), starts.end(),
                             [&](const auto &item) { return item.first == node; }))
                throw std::runtime_error("suffix seed index missed a known graph node");
            for (auto [hit, cost] : starts) {
                auto word = label(hit);
                size_t actual = 0;
                for (size_t j = 0; j < seed.size(); ++j)
                    actual += word[word.size() - seed.size() + j] != seed[j];
                if (actual != cost || actual > errors)
                    throw std::runtime_error("seed index returned a wrong alignment");
            }
            ++checked;
        }
    }
    std::cout << "Verified " << checked << " exact/mutated seed searches against real graph labels\n"
              << std::flush;
}

void ProtospacerWorkflow::load_sources() {
    auto add = [&](SpacerSource source) {
        if (!dna(source.spacer) || !dna(source.repeat))
            throw std::runtime_error("non-ACGT source observation: " + source.id);
        source.original_spacer = source.spacer;
        source.original_repeat = source.repeat;
        sources_.push_back(std::move(source));
    };
    fs::path table = options_.arrays / "CRISPR_graph_paths.tsv";
    if (fs::is_regular_file(table)) {
        std::ifstream input(table);
        std::string line;
        auto split = [](const std::string &s) {
            std::vector<std::string> v;
            std::stringstream stream(s);
            std::string x;
            while (std::getline(stream, x, '\t'))
                v.push_back(x);
            if (!s.empty() && s.back() == '\t')
                v.emplace_back();
            return v;
        };
        if (!std::getline(input, line))
            throw std::runtime_error("empty graph-path table");
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        auto header = split(line);
        std::map<std::string, size_t> col;
        for (size_t i = 0; i < header.size(); ++i)
            col[header[i]] = i;
        for (const auto &name : {"array_id", "observation_id", "repeat", "spacer"})
            if (!col.count(name))
                throw std::runtime_error("missing graph-path column: " + std::string(name));
        while (std::getline(input, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            if (line.empty())
                continue;
            auto row = split(line);
            if (row.size() != header.size())
                throw std::runtime_error("malformed graph-path row");
            SpacerSource s;
            s.array = row[col["array_id"]];
            s.id = s.array + ":" + row[col["observation_id"]];
            s.repeat = row[col["repeat"]];
            s.spacer = row[col["spacer"]];
            add(std::move(s));
        }
        if (input.bad())
            throw std::runtime_error("failed to read graph-path table");
    } else {
        std::vector<fs::path> files;
        for (const auto &f : fs::directory_iterator(options_.arrays)) {
            auto name = f.path().filename().string();
            if (f.is_regular_file() && name.rfind("CRISPR_Arrays_", 0) == 0 && f.path().extension() == ".txt")
                files.push_back(f.path());
        }
        if (files.empty() && fs::is_regular_file(options_.arrays / "CRISPR_Arrays.txt"))
            files.push_back(options_.arrays / "CRISPR_Arrays.txt");
        std::sort(files.begin(), files.end());
        for (const auto &file : files) {
            std::ifstream input(file);
            std::string line, array, repeat;
            size_t observation = 0;
            if (!input)
                throw std::runtime_error("cannot read " + file.string());
            while (std::getline(input, line)) {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                if (line.empty() || line[0] == '#')
                    continue;
                if (line[0] == '>') {
                    std::istringstream h(line.substr(1));
                    h >> array;
                    repeat.clear();
                    observation = 0;
                    continue;
                }
                if (array.empty())
                    throw std::runtime_error("array sequence before header");
                std::istringstream fields(line);
                std::string first, spacer, extra;
                fields >> first;
                if (!(fields >> spacer)) {
                    if (!dna(first))
                        throw std::runtime_error("invalid array repeat");
                    repeat = first;
                    continue;
                }
                if (fields >> extra || repeat.empty())
                    throw std::runtime_error("malformed spacer observation");
                if (first.find('-') != std::string::npos) {
                    if (first.size() != repeat.size())
                        throw std::runtime_error("repeat variant has incorrect length");
                    for (size_t i = 0; i < first.size(); ++i)
                        if (first[i] == '-')
                            first[i] = repeat[i];
                }
                SpacerSource s;
                s.array = file.filename().string() + ":" + array;
                s.id = s.array + ":" + std::to_string(++observation);
                s.spacer = spacer;
                s.repeat = first;
                add(std::move(s));
            }
            if (input.bad())
                throw std::runtime_error("array-report read failed");
        }
    }
    std::set<std::string> ids;
    for (const auto &s : sources_)
        if (!ids.insert(s.id).second)
            throw std::runtime_error("duplicate source ID: " + s.id);
}

void ProtospacerWorkflow::correct_boundaries() {
    if (!options_.correct_boundaries)
        return;
    std::map<std::string, std::set<std::string>> groups;
    for (const auto &s : sources_)
        groups[s.array].insert(s.spacer);
    std::map<std::string, std::pair<std::string, std::string>> motifs;
    for (const auto &[array, spacers] : groups) {
        if (spacers.size() < 5)
            continue;
        size_t shortest = std::numeric_limits<size_t>::max();
        for (const auto &s : spacers)
            shortest = std::min(shortest, s.size());
        for (bool right : {false, true}) {
            std::string conserved;
            for (size_t n = 1; n < shortest; ++n) {
                std::map<std::string, size_t> counts;
                for (const auto &s : spacers)
                    ++counts[right ? s.substr(s.size() - n) : s.substr(0, n)];
                auto top = std::max_element(counts.begin(), counts.end(),
                                            [](const auto &a, const auto &b) { return a.second < b.second; });
                if (double(top->second) + 1e-12 < options_.boundary_support * spacers.size())
                    break;
                conserved = top->first;
            }
            if (conserved.size() >= 6) {
                if (right)
                    motifs[array].second = conserved;
                else
                    motifs[array].first = conserved;
            }
        }
    }
    std::vector<std::pair<size_t, size_t>> proposed;
    std::map<std::string, std::set<std::pair<size_t, size_t>>> consistent;
    for (const auto &s : sources_) {
        const auto &m = motifs[s.array];
        size_t begin = 0, end = s.spacer.size();
        if (!m.first.empty() && s.spacer.compare(0, m.first.size(), m.first) == 0)
            begin = m.first.size();
        if (!m.second.empty() &&
            s.spacer.compare(s.spacer.size() - m.second.size(), m.second.size(), m.second) == 0)
            end -= m.second.size();
        if (end < begin || end - begin < size_t(graph_->k())) {
            begin = 0;
            end = s.spacer.size();
        }
        proposed.push_back({begin, end});
        auto key = canonical(s.spacer);
        consistent[key].insert(key == s.spacer
                                   ? std::make_pair(begin, end)
                                   : std::make_pair(s.spacer.size() - end, s.spacer.size() - begin));
    }
    for (size_t i = 0; i < sources_.size(); ++i) {
        auto &s = sources_[i];
        auto [begin, end] = proposed[i];
        if (consistent[canonical(s.spacer)].size() != 1)
            continue;
        s.trim_left = begin;
        s.trim_right = s.spacer.size() - end;
        s.repeat = s.spacer.substr(end) + s.repeat + s.spacer.substr(0, begin);
        s.spacer = s.spacer.substr(begin, end - begin);
    }
}

void ProtospacerWorkflow::prepare_source_nodes() {
    // Index only the reported source arrays, not the whole metagenome.
    const size_t k = graph_->k();
    for (const auto &source : sources_) {
        auto sequence = source.repeat + source.spacer + source.repeat;
        for (const auto &oriented : {sequence, reverse_complement(sequence)})
            for (size_t i = 0; i + k <= oriented.size(); ++i) {
                auto node = find_node(oriented.substr(i, k));
                if (node != missing) {
                    array_nodes_[source.array].insert(node);
                    // A pure spacer k-mer can be a genuine target elsewhere.
                    // Only repeat/junction overlap indicates array context.
                    bool pure_spacer = i >= source.repeat.size() &&
                                       i + k <= source.repeat.size() + source.spacer.size();
                    if (!pure_spacer) array_repeat_nodes_[source.array].insert(node);
                }
            }
    }
}

void ProtospacerWorkflow::load_supplied_targets() {
    if (options_.protospacers.empty())
        return;
    std::ifstream input(options_.protospacers);
    std::string line;
    std::set<std::string> ids;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (line.front() == '>') {
            std::istringstream header(line.substr(1));
            std::string id;
            header >> id;
            if (id.empty() || !ids.insert(id).second)
                throw std::runtime_error("empty or duplicate supplied target ID");
            supplied_.emplace_back(id, "");
        } else {
            if (supplied_.empty() || !dna(line))
                throw std::runtime_error("invalid supplied target FASTA");
            supplied_.back().second += line;
        }
    }
    if (input.bad() || supplied_.empty())
        throw std::runtime_error("empty or unreadable supplied target FASTA");
    auto array_sources = std::move(sources_);
    sources_.clear();
    for (const auto &[id, sequence] : supplied_) {
        if (!dna(sequence))
            throw std::runtime_error("empty supplied target: " + id);
        size_t matched = 0;
        std::set<std::tuple<std::string, std::string, std::string, std::string>> contexts;
        for (const auto &source : array_sources)
            for (bool reverse : {false, true}) {
                auto spacer = reverse ? reverse_complement(source.spacer) : source.spacer;
                auto repeat = reverse ? reverse_complement(source.repeat) : source.repeat;
                for (size_t at = spacer.find(sequence); at != std::string::npos;
                     at = spacer.find(sequence, at + 1)) {
                    auto left = spacer.substr(0, at), right = spacer.substr(at + sequence.size());
                    if (!contexts.emplace(source.array, repeat, left, right).second)
                        continue;
                    SpacerSource s;
                    s.id = id + ":" + std::to_string(++matched);
                    s.array = source.array;
                    s.spacer = s.original_spacer = sequence;
                    s.repeat = s.original_repeat = repeat;
                    s.source_left = left;
                    s.source_right = right;
                    sources_.push_back(std::move(s));
                }
            }
        if (!matched) {
            // Keep an explicit input outcome; do not invent missing source context.
            SpacerSource s;
            s.id = id + ":0";
            s.array = "unassociated";
            s.spacer = s.original_spacer = sequence;
            sources_.push_back(std::move(s));
        }
    }
}

void ProtospacerWorkflow::group_queries() {
    std::map<std::string, std::vector<size_t>> grouped;
    for (size_t i = 0; i < sources_.size(); ++i)
        grouped[canonical(sources_[i].spacer)].push_back(i);
    queries_.assign(grouped.begin(), grouped.end());
    query_source_nodes_.resize(queries_.size());
    query_repeat_nodes_.resize(queries_.size());
    for (size_t q = 0; q < queries_.size(); ++q)
        for (auto index : queries_[q].second) {
            const auto &nodes = array_nodes_[sources_[index].array];
            query_source_nodes_[q].insert(nodes.begin(), nodes.end());
            const auto &repeats = array_repeat_nodes_[sources_[index].array];
            query_repeat_nodes_[q].insert(repeats.begin(), repeats.end());
        }
    results_.resize(queries_.size());
    query_status_.resize(queries_.size(), "not_processed");
    query_work_.resize(queries_.size());
    mapping_work_.resize(queries_.size());
    query_limited_.resize(queries_.size());
    query_timings_.resize(queries_.size());
}

std::vector<TargetCandidate> ProtospacerWorkflow::map_spacer(size_t query, size_t &work,
                                                             bool &limited) const {
    const auto &sequence = queries_[query].first;
    const size_t k = graph_->k(), length = sequence.size();
    std::vector<TargetCandidate> found;
    if (length < k) {
        if (options_.protospacers.empty() || k - length > options_.max_flank)
            return found;
        // A supplied short target is an interval inside a full graph k-mer.
        // Keep the complete target unchanged; padding is explicit, not alignment.
        auto starts = find_starts(sequence, 0, work, limited);
        for (auto [node, errors] : starts) {
            TargetCandidate candidate;
            candidate.query = query;
            candidate.end = length;
            candidate.anchor = {node};
            candidate.anchor_prefix = k - length;
            candidate.target = sequence;
            candidate.identity = candidate.coverage = 1.;
            found.push_back(std::move(candidate));
            if (found.size() >= options_.max_placements) {
                limited |= starts.size() > found.size();
                break;
            }
        }
        return found;
    }
    const size_t minimum = std::max(k, size_t(std::ceil(options_.coverage * length - 1e-9)));
    std::map<std::pair<std::string, size_t>, std::vector<std::pair<uint64_t, size_t>>> cache;
    auto step = [&]() {
        if (interrupted_ || work++ >= options_.placement_work) {
            limited = true;
            return false;
        }
        return true;
    };
    std::set<std::pair<size_t, std::vector<uint64_t>>> aligned_states;
    auto finish_right = [&](TargetCandidate start, size_t allowed) {
        std::vector<TargetCandidate> pending{std::move(start)};
        while (!pending.empty() && step()) {
            auto current = std::move(pending.back());
            pending.pop_back();
            if (!aligned_states.emplace(current.begin, current.anchor).second)
                continue;
            size_t aligned = current.end - current.begin;
            if (aligned >= minimum && current.mismatches <= aligned * (1 - options_.identity) + 1e-9) {
                auto candidate = current;
                candidate.target = spell(candidate.anchor);
                candidate.identity = 1. - double(candidate.mismatches) / aligned;
                candidate.coverage = double(aligned) / length;
                found.push_back(std::move(candidate));
                if (found.size() >= options_.max_placements * 256) {
                    limited = true;
                    return;
                }
            }
            if (current.end == length)
                continue;
            for (const auto &arc : neighbors(current.anchor.back(), false)) {
                size_t changes = current.mismatches + (added_base(arc, false) != sequence[current.end]);
                if (changes > allowed)
                    continue;
                auto next = current;
                next.anchor.push_back(arc.node);
                ++next.end;
                next.mismatches = changes;
                pending.push_back(std::move(next));
            }
        }
    };
    for (size_t begin = 0; begin + minimum <= length && !interrupted_; ++begin) {
        const size_t allowed = size_t(std::floor((1 - options_.identity) * (length - begin) + 1e-9));
        // An exact k-mer lookup remains available even if inexact mapping reaches
        // its cap. It is checked with the same connected alignment afterwards.
        auto exact = find_node(sequence.substr(begin, k));
        if (exact != missing && work < options_.placement_work) {
            TargetCandidate start;
            start.query = query;
            start.begin = begin;
            start.end = begin + k;
            start.anchor = {exact};
            finish_right(std::move(start), allowed);
        }
        if (work >= options_.placement_work)
            break;
        // Partition the minimum covered interval. If at most E substitutions
        // occur, one of these blocks has at most floor(E/parts) substitutions.
        // Unlike a single whole-kmer search, this includes seeds at both ends.
        size_t parts = std::max(size_t(2), (minimum + k - 1) / k);
        for (size_t part = 0; part < parts && work < options_.placement_work; ++part) {
            size_t first = begin + minimum * part / parts, last = begin + minimum * (part + 1) / parts;
            auto key = std::make_pair(sequence.substr(first, last - first), allowed / parts);
            if (!cache.count(key)) {
                bool seed_limited = false;
                cache[key] = find_starts(key.first, key.second, work, seed_limited);
                limited |= seed_limited;
            }
            for (auto [node, seed_errors] : cache[key]) {
                if (!step())
                    break;
                (void)seed_errors;
                const auto word = label(node);
                int64_t at = int64_t(last) - int64_t(k);
                size_t errors = 0;
                for (size_t j = 0; j < word.size(); ++j)
                    if (at + int64_t(j) >= int64_t(begin))
                        errors += word[j] != sequence[at + j];
                if (errors > allowed)
                    continue;
                TargetCandidate start;
                start.query = query;
                start.begin = begin;
                start.end = last;
                start.mismatches = errors;
                start.anchor = {node};
                std::vector<TargetCandidate> pending{std::move(start)};
                while (!pending.empty() && step()) {
                    auto current = std::move(pending.back());
                    pending.pop_back();
                    if (at < int64_t(begin)) {
                        // The seed k-mer begins before the aligned query interval.
                        // Move forward until its complete label belongs to the query.
                        if (current.end == begin + k) {
                            current.anchor = {current.anchor.back()};
                            finish_right(std::move(current), allowed);
                            continue;
                        }
                        for (const auto &arc : neighbors(current.anchor.back(), false)) {
                            size_t changes =
                                current.mismatches + (added_base(arc, false) != sequence[current.end]);
                            if (changes > allowed)
                                continue;
                            auto next = current;
                            next.anchor = {arc.node};
                            ++next.end;
                            next.mismatches = changes;
                            pending.push_back(std::move(next));
                        }
                    } else {
                        // Walk backwards from this seed to the interval start.
                        size_t current_begin = last - k - (current.anchor.size() - 1);
                        if (current_begin == begin) {
                            std::reverse(current.anchor.begin(), current.anchor.end());
                            finish_right(std::move(current), allowed);
                            continue;
                        }
                        for (const auto &arc : neighbors(current.anchor.back(), true)) {
                            size_t changes =
                                current.mismatches + (added_base(arc, true) != sequence[current_begin - 1]);
                            if (changes > allowed)
                                continue;
                            auto next = current;
                            next.anchor.push_back(arc.node);
                            next.mismatches = changes;
                            pending.push_back(std::move(next));
                        }
                    }
                }
                if (work >= options_.placement_work)
                    break;
            }
        }
    }
    std::sort(found.begin(), found.end(), [](const auto &a, const auto &b) {
        if (a.coverage != b.coverage)
            return a.coverage > b.coverage;
        if (a.identity != b.identity)
            return a.identity > b.identity;
        return std::tie(a.begin, a.anchor) < std::tie(b.begin, b.anchor);
    });
    std::vector<TargetCandidate> unique;
    for (auto &candidate : found) {
        bool contained = false;
        for (const auto &other : unique)
            if (candidate.begin >= other.begin && candidate.end <= other.end &&
                std::search(other.anchor.begin(), other.anchor.end(), candidate.anchor.begin(),
                            candidate.anchor.end()) != other.anchor.end()) {
                contained = true;
                break;
            }
        if (!contained)
            unique.push_back(std::move(candidate));
        if (unique.size() > options_.max_placements) {
            unique.pop_back();
            limited = true;
            break;
        }
    }
    return unique;
}
int ProtospacerWorkflow::source_context(size_t query, const TargetCandidate &candidate,
                                        const std::string &outward, bool backwards) const {
    bool pending = false;
    const auto &sequence = queries_[query].first;
    for (size_t index : queries_[query].second) {
        const auto &source = sources_[index];
        std::string repeat = source.repeat;
        auto left = source.source_left, right = source.source_right;
        if (source.spacer != sequence) {
            repeat = reverse_complement(repeat);
            std::swap(left, right);
            left = reverse_complement(left);
            right = reverse_complement(right);
        }
        auto remainder =
            backwards ? left + sequence.substr(0, candidate.begin) : sequence.substr(candidate.end) + right;
        if (backwards) {
            std::reverse(repeat.begin(), repeat.end());
            std::reverse(remainder.begin(), remainder.end());
        }
        if (repeat.empty())
            return 0;
        int state = repeat_context(outward, remainder, repeat, options_.repeat_identity);
        if (state == 1)
            return 1;
        pending |= state == 0;
    }
    return pending ? 0 : 2;
}

size_t ProtospacerWorkflow::source_clearance(size_t query, const TargetCandidate &c, bool backwards) const {
    size_t end = 0;
    for (auto index : queries_[query].second) {
        const auto &s = sources_[index];
        bool flipped = s.spacer != queries_[query].first;
        size_t remainder = backwards ? c.begin : queries_[query].first.size() - c.end;
        remainder += (backwards != flipped) ? s.source_left.size() : s.source_right.size();
        size_t edits = size_t(std::floor((1 - options_.repeat_identity) * s.repeat.size() + 1e-9));
        end = std::max(end, remainder + s.repeat.size() + edits + size_t(graph_->k()));
    }
    size_t padding = backwards ? c.anchor_prefix : 0;
    return end > padding ? end - padding : 0;
}

void ProtospacerWorkflow::trim_walk(ExtensionWalk &walk, size_t length, const std::string &reason) const {
    length = std::min(length, walk.bases.size());
    walk.bases.resize(length);
    if (walk.all_counts.size() > length) walk.all_counts.resize(length);
    if (!walk.nodes.empty())
        walk.nodes.resize(length + 1);
    walk.stop = reason;
    if (walk.ambiguous && length <= walk.first_ambiguous_base) {
        walk.ambiguous = false;
        walk.competing_routes = walk.first_ambiguous_base = 0;
        walk.runner_up_score = -1.;
    }
    while (!walk.sample_positions.empty() && walk.sample_positions.back() > length) {
        walk.sample_positions.pop_back();
        walk.counts.pop_back();
        walk.capacities.pop_back();
    }
    if (length < walk.clear_at) {
        walk.external = false;
        walk.clear_at = 0;
    }
    if (!length) {
        walk.crosses_source = false;
        walk.source_until = walk.source_last = 0;
        walk.seed_flow = walk.remaining_flow = 0.;
        walk.flow_decisions.clear();
    } else if (!walk.flow_decisions.empty()) {
        walk.flow_decisions.erase(std::remove_if(walk.flow_decisions.begin(), walk.flow_decisions.end(),
            [&](const FlowDecision &d) { return d.position > length || (d.position == length && d.chosen >= 0); }),
            walk.flow_decisions.end());
        walk.remaining_flow = walk.seed_flow;
        for (const auto &d : walk.flow_decisions)
            if (d.chosen >= 0) walk.remaining_flow = d.remaining;
    }
    walk.flow = walk.capacities.empty() ? 0 : walk.capacities.back();
    walk.profile_samples = 0;
    walk.deviation = 0;
    walk.transition_energy = 0;
    walk.transitions = 0;
    if (walk.branch_costs.size() > length)
        walk.branch_costs.resize(length);
    walk.branch_loss = 0.;
    walk.branch_observations = 0;
    update_profile(walk, true);
}

void ProtospacerWorkflow::update_profile(ExtensionWalk &walk, bool) const {
    // Consecutive observed edge counts, not unitig medians. Unitig boundaries
    // must not change the evidence for an otherwise identical graph path.
    while (walk.profile_samples < walk.counts.size()) {
        size_t i = walk.profile_samples++;
        if (i && (walk.sample_positions.empty() ||
                  walk.sample_positions[i] == walk.sample_positions[i - 1] + 1)) {
            double change = std::log2((walk.counts[i] + 1.) / (walk.counts[i - 1] + 1.));
            walk.transition_energy += change * change;
            ++walk.transitions;
            // A branch is scored only when both ends have downstream evidence.
            // This excludes the masked source-exit fork, but includes downstream
            // forks observed while waiting to confirm source clearance.
            size_t position = walk.sample_positions.empty() ? i + 1 : walk.sample_positions[i];
            if (position && position <= walk.branch_costs.size()) {
                double cost = walk.branch_costs[position - 1];
                walk.branch_loss += cost;
                walk.branch_observations += cost > 0.;
            }
        }
    }
    walk.deviation = walk.transitions ? walk.transition_energy / walk.transitions : 0.;
}

double ProtospacerWorkflow::walk_score(const ExtensionWalk &walk) const {
    // Conditional graph-walk surprisal per scored downstream edge. Unbranched
    // edges contribute zero. Normalize by sequence length, not fork count, so
    // extra easy forks cannot dilute a weak choice. Raw count continuity is a
    // separate tie-break, never a reward traded against branch support.
    // This score is not a probability of genomic or viral correctness.
    return walk.counts.empty() ? 0. : walk.branch_loss / walk.counts.size();
}

void ProtospacerWorkflow::rank_walks(std::vector<ExtensionWalk> &walks, size_t limit, bool diverse) const {
    if (walks.empty())
        return;
    auto complete = [&](const auto &w) {
        return w.bases.size() + w.anchor_padding >= options_.max_flank;
    };
    std::sort(walks.begin(), walks.end(), [&](const auto &a, const auto &b) {
        if (a.external != b.external)
            return a.external > b.external;
        if (complete(a) != complete(b))
            return complete(a) > complete(b);
        if (walk_score(a) != walk_score(b))
            return walk_score(a) < walk_score(b);
        if (a.deviation != b.deviation)
            return a.deviation < b.deviation;
        if (a.bases.size() != b.bases.size())
            return a.bases.size() > b.bases.size();
        return a.nodes < b.nodes;
    });
    if (diverse && walks.size() > limit) {
        std::set<std::tuple<uint64_t, bool, bool>> endpoints;
        std::vector<size_t> chosen;
        std::vector<bool> used(walks.size());
        // Keep distinct graph continuations before spending remaining slots on
        // rejoined variants. Each retained state keeps its own full history.
        for (size_t i = 0; i < walks.size() && chosen.size() < limit; ++i)
            if (endpoints.emplace(walks[i].nodes.back(), walks[i].external, walks[i].crosses_source).second) {
                chosen.push_back(i);
                used[i] = true;
            }
        for (size_t i = 0; i < walks.size() && chosen.size() < limit; ++i)
            if (!used[i]) {
                chosen.push_back(i);
                used[i] = true;
            }
        std::vector<ExtensionWalk> reordered;
        reordered.reserve(walks.size());
        for (size_t i : chosen)
            reordered.push_back(std::move(walks[i]));
        for (size_t i = 0; i < walks.size(); ++i)
            if (!used[i])
                reordered.push_back(std::move(walks[i]));
        walks = std::move(reordered);
    }
}

void ProtospacerWorkflow::save_walks(const TargetCandidate &candidate, bool backwards,
                                     const std::vector<ExtensionWalk> &walks, size_t &serial) const {
    if (!options_.save_search || walks.empty())
        return;
    auto query = numbered("Q", candidate.query + 1);
    auto path = options_.output / "search_paths" / (query + ".tsv.gz");
    bool header = !fs::exists(path);
    gzFile file = gzopen(path.c_str(), "ab1");
    if (!file)
        throw std::runtime_error("cannot open search audit: " + path.string());
    auto write = [&](const std::string &s) {
        if (gzwrite(file, s.data(), unsigned(s.size())) != int(s.size()))
            throw std::runtime_error("cannot write search audit: " + path.string());
    };
    try {
        if (header)
            write("path_id\tquery_id\tplacement\tside\tquery_begin\tquery_"
                  "end\tidentity\tquery_coverage\t"
                  "flank_bp\texternal\tclear_"
                  "at\tstop\tlimited\tflow\tsequence\tnodes_outward\t"
                  "count_positions\tcounts\tcoverage_deviation\tanchor_prefix_bp\tcrosses_source_"
                  "array\tidentified_before_extension\ttransition_energy\tadjacent_observations\tbranch_loss\tbranch_observations\tbranch_costs\n");
        for (size_t i = 0; i < walks.size(); ++i) {
            const auto &walk = walks[i];
            auto flank = walk.bases;
            if (backwards)
                std::reverse(flank.begin(), flank.end());
            std::ostringstream row;
            auto side = backwards ? "L" : "R";
            row << query << '_' << (candidate.allow_crossing ? "C" : "S") << '_'
                << numbered("A", candidate.placement + 1) << '_' << side << '_' << numbered("P", ++serial)
                << '\t' << query << '\t' << candidate.placement << '\t' << side << '\t' << candidate.begin
                << '\t' << candidate.end << '\t' << candidate.identity << '\t' << candidate.coverage << '\t'
                << walk.bases.size() << '\t' << walk.external << '\t' << walk.clear_at << '\t' << walk.stop
                << '\t' << walk.limited << '\t' << walk.flow << '\t'
                << (backwards ? flank + spell(candidate.anchor) : spell(candidate.anchor) + flank);
            auto list = [&](const auto &values) {
                row << '\t';
                for (size_t j = 0; j < values.size(); ++j) {
                    if (j)
                        row << ',';
                    row << values[j];
                }
            };
            list(walk.nodes);
            list(walk.sample_positions);
            list(walk.counts);
            row << '\t' << walk.deviation << '\t' << candidate.anchor_prefix << '\t' << walk.crosses_source
                << '\t' << candidate.identified << '\t' << walk.transition_energy << '\t' << walk.transitions
                << '\t' << walk.branch_loss << '\t' << walk.branch_observations;
            list(walk.branch_costs);
            row << '\n';
            write(row.str());
        }
    } catch (...) {
        gzclose(file);
        throw;
    }
    if (gzclose(file) != Z_OK)
        throw std::runtime_error("cannot finish search audit: " + path.string());
}

bool ProtospacerWorkflow::extension_barrier(ExtensionWalk &walk, const TargetCandidate &candidate,
                                            bool backwards) const {
    if (!options_.filter_artifacts)
        return false;
    // Include the anchor boundary: a technical motif can straddle seed and flank.
    std::string tail;
    if (walk.bases.size() < 64) {
        tail = candidate.anchor_prefix ? label(candidate.anchor.front()) : candidate.target;
        if (backwards)
            std::reverse(tail.begin(), tail.end());
    }
    tail += walk.bases.substr(walk.bases.size() > 64 ? walk.bases.size() - 64 : 0);
    if (tail.size() > 64)
        tail.erase(0, tail.size() - 64);
    size_t position = adapter_position(tail, backwards);
    if (position != std::string::npos) {
        size_t retained =
            walk.bases.size() + position > tail.size() ? walk.bases.size() + position - tail.size() : 0;
        trim_walk(walk, retained, "adapter_sequence");
        return true;
    }
    if (tail.size() == 64 && walk.bases.size() % 8 == 0 && !informative(tail)) {
        trim_walk(walk, walk.bases.size() > 64 ? walk.bases.size() - 64 : 0, "low_complexity_extension");
        return true;
    }
    return false;
}

void ProtospacerWorkflow::filter_selected_extension(TargetCandidate &candidate) const {
    if (!candidate.detected || !candidate.extended) return;
    const bool backwards = !candidate.left.bases.empty();
    auto &walk = backwards ? candidate.left : candidate.right;
    // One read contributes read_length-k consecutive (k+1)-mer observations.
    // Use complete disjoint blocks of unmasked downstream counts only.
    const size_t k = options_.context_bases - 1;
    const size_t span = std::max(options_.context_bases,
                                options_.read_length > k ? options_.read_length - k : size_t(1));
    candidate.coverage_blocks = walk.counts.size() / span;
    for (size_t block = 0; block < candidate.coverage_blocks; ++block) {
        const size_t start = block * span;
        std::vector<double> counts(walk.counts.begin() + start, walk.counts.begin() + start + span);
        auto middle = counts.begin() + counts.size() / 2;
        std::nth_element(counts.begin(), middle, counts.end());
        double median = *middle;
        if (counts.size() % 2 == 0)
            median = (median + *std::max_element(counts.begin(), middle)) / 2.;
        candidate.coverage_min = block ? std::min(candidate.coverage_min, median) : median;
        candidate.coverage_max = block ? std::max(candidate.coverage_max, median) : median;
    }
    if (options_.max_coverage_fold == 0. || candidate.coverage_blocks < 2 ||
        (candidate.coverage_min > 0. && candidate.coverage_max <= options_.max_coverage_fold * candidate.coverage_min))
        return;
    // Validate the independently selected route; do not search for a different
    // smooth-looking route after the strongest one fails this check.
    // Retain the exact target and its original discovery evidence. Keep the
    // withheld route/counts in a separate audit, never in the accepted FASTA.
    candidate.withheld_sequence = candidate.sequence;
    candidate.withheld_left = backwards;
    candidate.withheld = std::move(walk);
    walk = ExtensionWalk{};
    walk.stop = "coverage_discontinuity";
    candidate.sequence = (candidate.anchor_prefix ? label(candidate.anchor.front()).substr(0, candidate.anchor_prefix) : "") + candidate.target;
    candidate.extended = false;
    candidate.status = "detected";
    candidate.reason = "coverage_discontinuity";
}

bool ProtospacerWorkflow::advance_walk(ExtensionWalk &walk, const GraphArc &arc, size_t query,
                                       const TargetCandidate &candidate, bool backwards, double branch_cost) const {
    if (std::find(candidate.anchor.begin(), candidate.anchor.end(), arc.node) != candidate.anchor.end()) {
        walk.stop = "anchor_revisit";
        return false;
    }
    bool source_node = candidate.allow_crossing && query_source_nodes_[query].count(arc.node);
    if (source_node && std::find(walk.nodes.begin(), walk.nodes.end(), arc.node) != walk.nodes.end()) {
        walk.stop = "source_node_revisit";
        return false;
    }
    bool repeat_node = candidate.allow_crossing && query_repeat_nodes_[query].count(arc.node);
    if (repeat_node && walk.external && candidate.allow_crossing) {
        walk.stop = "source_array_reentry";
        return false;
    }
    if (size_t(std::count(walk.nodes.begin(), walk.nodes.end(), arc.node)) >= options_.max_visits) {
        walk.stop = "repeat_visit_limit";
        return false;
    }
    walk.nodes.push_back(arc.node);
    walk.bases += added_base(arc, backwards);
    walk.all_counts.push_back(arc.count);
    walk.branch_costs.push_back(branch_cost);
    if (extension_barrier(walk, candidate, backwards))
        return false;
    if (!walk.external) {
        auto outward = walk.bases;
        if (backwards && candidate.anchor_prefix) {
            auto padding = label(candidate.anchor.front()).substr(0, candidate.anchor_prefix);
            std::reverse(padding.begin(), padding.end());
            outward = padding + outward;
        }
        int state = walk.crosses_source ? 1 : source_context(query, candidate, outward, backwards);
        if (state == 1) {
            if (!candidate.allow_crossing) {
                trim_walk(walk, 0, "source_array_context");
                return false;
            }
            walk.crosses_source = true;
            walk.source_until = source_clearance(query, candidate, backwards);
            if (repeat_node)
                walk.source_last = walk.bases.size();
            if (walk.bases.size() >= walk.source_until &&
                walk.bases.size() >= walk.source_last + options_.context_bases) {
                walk.external = true;
                walk.clear_at = walk.bases.size();
            }
            if (!walk.external && walk.bases.size() >= walk.source_until + options_.context_bases) {
                walk.stop = "source_passage_not_exited";
                return false;
            }
        }
        if (state == 2) {
            walk.external = true;
            walk.clear_at = walk.bases.size();
        }
    }
    if (walk.external) {
        size_t first = walk.crosses_source ? std::max(walk.source_last,
            walk.source_until > size_t(graph_->k()) ? walk.source_until - graph_->k() : 0)
                                          : size_t(graph_->k());
        size_t recorded = walk.sample_positions.empty() ? first : walk.sample_positions.back();
        for (size_t i = recorded; i < walk.all_counts.size(); ++i) {
            if (walk.all_counts[i] >= 65535) {
                trim_walk(walk, i, "saturated_multiplicity");
                return false;
            }
            walk.counts.push_back(walk.all_counts[i]);
            walk.sample_positions.push_back(i + 1);
            walk.flow = walk.capacities.empty() ? walk.all_counts[i] : std::min(walk.flow, walk.all_counts[i]);
            walk.capacities.push_back(walk.flow);
        }
    }
    if (walk.bases.size() + (backwards ? candidate.anchor_prefix : 0) >= options_.max_flank) {
        walk.stop = "maximum_flank";
        return false;
    }
    return true;
}

bool ProtospacerWorkflow::retain_completed(std::vector<ExtensionWalk> &walks) const {
    walks.erase(std::remove_if(walks.begin(), walks.end(),
                               [&](const auto &w) {
                                   return !w.external ||
                                          w.bases.size() < std::max(options_.context_bases - 1, w.clear_at) +
                                                               options_.context_bases;
                               }),
                walks.end());
    bool has_extension = std::any_of(walks.begin(), walks.end(),
                                     [&](const auto &w) { return w.bases.size() + w.anchor_padding >= options_.min_flank && w.counts.size() >= options_.min_external; });
    if (has_extension)
        walks.erase(std::remove_if(walks.begin(), walks.end(),
                                   [&](const auto &w) { return w.bases.size() + w.anchor_padding < options_.min_flank || w.counts.size() < options_.min_external; }),
                    walks.end());
    if (walks.size() <= options_.beam_width * 4)
        return false;
    rank_walks(walks, options_.beam_width * 4, false);
    walks.resize(options_.beam_width * 4);
    return true;
}

bool ProtospacerWorkflow::same_route(const ExtensionWalk &a, const ExtensionWalk &b) const {
    // A small bubble that rejoins the same oriented node is local sequence
    // uncertainty. It must not erase all the supported sequence after the join.
    const size_t bubble = options_.context_bases + 1;
    size_t i = 0, j = 0;
    while (i < a.nodes.size() && j < b.nodes.size()) {
        if (a.nodes[i] == b.nodes[j]) {
            ++i;
            ++j;
            continue;
        }
        size_t ai = a.nodes.size(), bj = b.nodes.size(), distance = 2 * bubble;
        for (size_t x = i; x < std::min(a.nodes.size(), i + bubble); ++x)
            for (size_t y = j; y < std::min(b.nodes.size(), j + bubble); ++y)
                if (a.nodes[x] == b.nodes[y] && x - i + y - j < distance) {
                    ai = x;
                    bj = y;
                    distance = x - i + y - j;
                }
        if (ai == a.nodes.size() || bj == b.nodes.size())
            return false;
        // Require real shared downstream context, not only a repeat junction.
        size_t shared = 0;
        while (ai + shared < a.nodes.size() && bj + shared < b.nodes.size() &&
               a.nodes[ai + shared] == b.nodes[bj + shared])
            ++shared;
        if (shared < options_.context_bases)
            return false;
        i = ai + shared;
        j = bj + shared;
    }
    return true;
}

ExtensionWalk ProtospacerWorkflow::choose_walk(std::vector<ExtensionWalk> walks) const {
    for (auto &walk : walks)
        update_profile(walk, true);
    retain_completed(walks);
    if (walks.empty()) {
        ExtensionWalk empty;
        empty.stop = "no_supported_external_context";
        return empty;
    }
    bool full = std::any_of(walks.begin(), walks.end(), [&](const auto &w) {
        return w.bases.size() + w.anchor_padding >= options_.max_flank;
    });
    if (full)
        walks.erase(std::remove_if(walks.begin(), walks.end(),
                                   [&](const auto &w) {
                                       return w.bases.size() + w.anchor_padding < options_.max_flank;
                                   }),
                    walks.end());
    rank_walks(walks, walks.size(), false);
    ExtensionWalk best = walks.front();
    size_t common = best.nodes.size();
    for (size_t i = 1; i < walks.size(); ++i) {
        const auto &other = walks[i];
        if (walk_score(other) > walk_score(best) * options_.branch_ratio + 1e-12)
            continue;
        if (same_route(best, other))
            continue;
        size_t n = 0;
        while (n < std::min(best.nodes.size(), other.nodes.size()) && best.nodes[n] == other.nodes[n])
            ++n;
        if (n == std::min(best.nodes.size(), other.nodes.size()))
            continue;
        common = std::min(common, n);
        ++best.competing_routes;
        if (best.runner_up_score < 0 || walk_score(other) < best.runner_up_score)
            best.runner_up_score = walk_score(other);
    }
    best.ambiguous = common < best.nodes.size();
    best.first_ambiguous_base = best.ambiguous ? (common ? common - 1 : 0) : 0;
    // A repeated source junction is unphased even when the other arm is too
    // short to enter the completed-candidate pool. Do not report certainty
    // merely because those alternatives were filtered later.
    if (best.crosses_source && best.pre_context_fork) {
        best.first_ambiguous_base = best.ambiguous
                                        ? std::min(best.first_ambiguous_base, best.first_context_fork)
                                        : best.first_context_fork;
        best.ambiguous = true;
    }
    return best;
}
ExtensionWalk ProtospacerWorkflow::extend_flank(size_t query, const TargetCandidate &candidate,
                                                bool backwards, size_t &work, size_t work_limit) const {
    ExtensionWalk seed;
    seed.anchor_padding = backwards ? candidate.anchor_prefix : 0;
    seed.nodes = {backwards ? candidate.anchor.front() : candidate.anchor.back()};
    std::vector<ExtensionWalk> frontier{seed}, finished;
    std::map<uint64_t, std::vector<GraphArc>> adjacency;
    bool limited = false;
    size_t serial = 0;
    // Every frontier contains paths at exactly the same physical distance.
    // Compression into unitigs cannot let one path jump ahead of its competitors.
    while (!frontier.empty() && !interrupted_) {
        std::vector<ExtensionWalk> next, ended, pruned;
        for (auto &current : frontier) {
            if (work >= work_limit) {
                current.limited = true;
                current.stop = "query_work_limit";
                ended.push_back(std::move(current));
                limited = true;
                continue;
            }
            auto node = current.nodes.back();
            if (!adjacency.count(node)) {
                ++work;
                adjacency.emplace(node, neighbors(node, backwards));
            }
            const auto &arcs = adjacency.at(node);
            if (arcs.empty()) {
                current.stop = "dead_end";
                ended.push_back(std::move(current));
                continue;
            }
            double total = 0.;
            for (const auto &arc : arcs) total += arc.count;
            for (const auto &arc : arcs) {
                auto child = current;
                if (arcs.size() > 1 && !current.external) {
                    if (!child.pre_context_fork) child.first_context_fork = current.bases.size();
                    child.pre_context_fork = true;
                }
                if (work++ >= work_limit) {
                    child.stop = "query_work_limit";
                    child.limited = limited = true;
                    ended.push_back(std::move(child));
                    continue;
                }
                ++child.work;
                // Store the raw choice cost now; update_profile applies the
                // same downstream mask to branch and count-change evidence.
                double cost = arcs.size() > 1 && arc.count > 0
                                ? -std::log2(arc.count / total) : 0.;
                bool continuing = advance_walk(child, arc, query, candidate, backwards, cost);
                update_profile(child);
                if (continuing)
                    next.push_back(std::move(child));
                else
                    ended.push_back(std::move(child));
            }
        }
        // Do not compare an unscored source passage with mature downstream
        // paths. Give it a separate bounded frontier until one read's effective
        // k-mer span has been observed outside the source region.
        size_t evidence_span = std::max(options_.context_bases,
            options_.read_length > size_t(graph_->k()) ? options_.read_length - graph_->k() : size_t(1));
        std::vector<ExtensionWalk> pending, mature;
        for (auto &walk : next)
            (walk.external && walk.counts.size() >= evidence_span ? mature : pending).push_back(std::move(walk));
        next.clear();
        auto retain = [&](std::vector<ExtensionWalk> &pool, size_t width, const char *reason) {
            rank_walks(pool, width, true);
            if (pool.size() > width) {
                for (size_t i = width; i < pool.size(); ++i) {
                    pool[i].limited = true;
                    pool[i].stop = reason;
                    pruned.push_back(std::move(pool[i]));
                }
                pool.resize(width);
                limited = true;
            }
            for (auto &walk : pool) next.push_back(std::move(walk));
        };
        retain(pending, options_.context_width, "context_path_limit");
        retain(mature, options_.beam_width, "beam_limit");
        save_walks(candidate, backwards, pruned, serial);
        save_walks(candidate, backwards, ended, serial);
        for (auto &walk : ended)
            finished.push_back(std::move(walk));
        if (retain_completed(finished)) limited = true;
        frontier = std::move(next);
    }
    for (auto &w : frontier) {
        w.stop = "interrupted";
        w.limited = true;
        finished.push_back(std::move(w));
    }
    auto chosen = choose_walk(std::move(finished));
    chosen.limited |= limited;
    return chosen;
}

void ProtospacerWorkflow::assemble_candidate(TargetCandidate &c) const {
    auto left = c.left.bases;
    std::reverse(left.begin(), left.end());
    auto padding = c.anchor_prefix ? label(c.anchor.front()).substr(0, c.anchor_prefix) : "";
    c.sequence = left + padding + c.target + c.right.bases;
    auto valid = [&](const ExtensionWalk &w) {
        return w.external &&
               w.bases.size() >= std::max(size_t(graph_->k()), w.clear_at) + options_.context_bases;
    };
    c.detected = c.identified || valid(c.left) || valid(c.right);
    if (!c.detected) {
        c.status = "rejected";
        c.reason = "no_independent_target_context";
        return;
    }
    // Crossing recovery must retain the independent context that established
    // discovery. The actual recovery route is saved in graph_paths.tsv.
    if (!c.allow_crossing) {
        c.context_left = std::string(padding.rbegin(), padding.rend()) + c.left.bases;
        c.context_right = c.right.bases;
        c.context_left_stop = c.left.stop;
        c.context_right_stop = c.right.stop;
        c.clear_left = c.left.clear_at + c.anchor_prefix;
        c.clear_right = c.right.clear_at;
    }
    // Validate the complete selected walk against the actual graph before output.
    std::vector<uint64_t> full;
    if (!c.left.nodes.empty())
        full.assign(c.left.nodes.rbegin(), c.left.nodes.rend());
    else
        full.push_back(c.anchor.front());
    full.insert(full.end(), c.anchor.begin() + 1, c.anchor.end());
    if (!c.right.nodes.empty())
        full.insert(full.end(), c.right.nodes.begin() + 1, c.right.nodes.end());
    if (spell(full) != c.sequence)
        throw std::runtime_error("selected path spelling mismatch");
    c.status = "detected";
    c.reason = "target_context_passed";
    if (!c.identified && !informative(c.target)) {
        c.detected = false;
        c.status = "unresolved";
        c.reason = "uninformative_target_sequence";
        return;
    }
    for (bool back : {true, false}) {
        auto &w = back ? c.left : c.right;
        if (w.bases.size() + w.anchor_padding < options_.min_flank || w.counts.size() < options_.min_external || w.flow <= 0 || !informative(w.bases))
            trim_walk(w, 0, "insufficient_supported_extension");
    }
    // The graph's pooled edge counts do not establish linkage between two
    // independently selected flanks. Report one actual supported flank.
    if (!c.left.bases.empty() && !c.right.bases.empty()) {
        std::vector<ExtensionWalk> sides{c.left, c.right};
        rank_walks(sides, 2, false);
        if (sides.front().nodes == c.left.nodes)
            trim_walk(c.right, 0, "alternative_flank");
        else
            trim_walk(c.left, 0, "alternative_flank");
    }
    left = c.left.bases;
    std::reverse(left.begin(), left.end());
    c.sequence = left + padding + c.target + c.right.bases;
    c.extended = !c.left.bases.empty() || !c.right.bases.empty();
    if (c.extended) {
        c.status = "extended";
        c.reason = "supported_extension";
    } else if (c.reason == "target_context_passed")
        c.reason = "target_detected_extension_unresolved";
}

void ProtospacerWorkflow::select_candidate(std::vector<TargetCandidate> &candidates) const {
    std::sort(candidates.begin(), candidates.end(), [&](const auto &a, const auto &b) {
        if (a.detected != b.detected)
            return a.detected > b.detected;
        if (a.extended != b.extended)
            return a.extended > b.extended;
        if (a.coverage != b.coverage)
            return a.coverage > b.coverage;
        if (a.identity != b.identity)
            return a.identity > b.identity;
        auto flank = [](const auto &c) -> const ExtensionWalk & {
            return c.left.bases.empty() ? c.right : c.left;
        };
        const auto &aw = flank(a), &bw = flank(b);
        bool af = aw.bases.size() + aw.anchor_padding >= options_.max_flank;
        bool bf = bw.bases.size() + bw.anchor_padding >= options_.max_flank;
        if (af != bf) return af > bf;
        if (walk_score(aw) != walk_score(bw)) return walk_score(aw) < walk_score(bw);
        if (aw.deviation != bw.deviation) return aw.deviation < bw.deviation;
        if (aw.bases.size() != bw.bases.size()) return aw.bases.size() > bw.bases.size();
        return a.anchor < b.anchor;
    });
    if (!candidates.empty() && candidates.front().detected) {
        auto &first = candidates.front();
        bool ambiguous = false;
        for (size_t i = 1; i < candidates.size(); ++i) {
            auto &other = candidates[i];
            if (other.detected && other.extended == first.extended && other.coverage == first.coverage &&
                other.identity == first.identity && other.target != first.target)
                ambiguous = true;
            other.detected = false;
            other.extended = false;
            other.status = "alternative";
            other.reason = "alternative_alignment";
        }
        if (ambiguous) {
            first.detected = false;
            first.extended = false;
            first.status = "unresolved";
            first.reason = "multiple_target_placements";
        }
    }
}

void ProtospacerWorkflow::process_query(size_t query) {
    if (options_.filter_artifacts && adapter_position(queries_[query].first) != std::string::npos) {
        query_status_[query] = "adapter_like_query";
        return;
    }
    if (options_.protospacers.empty() && !informative(queries_[query].first)) {
        query_status_[query] = "uninformative_query";
        return;
    }
    if (!options_.protospacers.empty() &&
        std::all_of(queries_[query].second.begin(), queries_[query].second.end(),
                    [&](size_t i) { return sources_[i].repeat.empty(); })) {
        query_status_[query] = "supplied_target_without_source_association";
        return;
    }
    auto start = std::chrono::steady_clock::now();
    size_t work = 0;
    bool limited = false;
    auto anchors = map_spacer(query, work, limited);
    query_timings_[query][0] =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    mapping_work_[query] = work;
    query_limited_[query] = limited;
    if (anchors.empty()) {
        query_status_[query] = queries_[query].first.size() < size_t(graph_->k()) ? "shorter_than_graph_k"
                               : limited                                          ? "placement_work_limit"
                                                                                  : "no_alignment";
        query_work_[query] = work;
        return;
    }
    std::vector<TargetCandidate> checked;
    for (size_t first = 0; first < anchors.size() && !interrupted_;) {
        size_t end = first + 1;
        while (end < anchors.size() && anchors[end].coverage == anchors[first].coverage &&
               anchors[end].identity == anchors[first].identity)
            ++end;
        size_t remaining = options_.max_work - std::min(work, options_.max_work);
        // Finish the best alignment tier before spending work on inferior mappings.
        // Reserve half the remaining budget for lower tiers if this one cannot
        // extend.
        size_t tier_limit = work + (end < anchors.size() ? remaining / 2 : remaining);
        bool has_complete_extension = false;
        for (size_t i = first; i < end && !interrupted_; ++i) {
            auto &candidate = anchors[i];
            candidate.placement = i;
            candidate.limited = limited;
            candidate.identified = !options_.protospacers.empty();
            // First preserve any extension supported without a source passage,
            // in both discovery and supplied-target modes.
            candidate.allow_crossing = false;
            if (options_.filter_artifacts &&
                adapter_position(candidate.anchor_prefix ? label(candidate.anchor.front())
                                                         : candidate.target) != std::string::npos) {
                candidate.reason = "adapter_like_target";
                checked.push_back(std::move(candidate));
                continue;
            }
            size_t budget = (tier_limit - std::min(work, tier_limit)) / (end - i);
            size_t left_limit = work + budget / 2, right_limit = work + budget;
            start = std::chrono::steady_clock::now();
            candidate.left = extend_flank(query, candidate, true, work, left_limit);
            candidate.right = extend_flank(query, candidate, false, work, right_limit);
            query_timings_[query][1] +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            candidate.limited |=
                candidate.left.limited || candidate.right.limited || work >= options_.max_work;
            query_limited_[query] |= candidate.limited;
            start = std::chrono::steady_clock::now();
            assemble_candidate(candidate);
            query_timings_[query][2] +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            has_complete_extension |= candidate.extended &&
                        std::max(candidate.left.bases.size() + candidate.anchor_prefix,
                                 candidate.right.bases.size()) >= options_.max_flank;
            checked.push_back(std::move(candidate));
        }
        if (has_complete_extension || work >= options_.max_work)
            break;
        first = end;
    }
    // Independent target detection precedes source crossing. An existing
    // partial/wrong strict extension must not prevent this comparison.
    if (options_.source_array_crossing && !interrupted_) {
        size_t eligible = std::count_if(checked.begin(), checked.end(), [](const auto &c) { return c.detected; });
        size_t extra = options_.max_work;
        size_t originals = checked.size();
        for (size_t i = 0; i < originals && eligible && !interrupted_; ++i) {
            if (!checked[i].detected) continue;
            auto candidate = checked[i];
            candidate.identified = candidate.allow_crossing = true;
            size_t budget = extra / eligible--;
            size_t before = work;
            start = std::chrono::steady_clock::now();
            candidate.left = extend_flank(query, candidate, true, work, before + budget / 2);
            candidate.right = extend_flank(query, candidate, false, work, before + budget);
            extra -= std::min(extra, work - before);
            query_timings_[query][1] +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            candidate.limited |= candidate.left.limited || candidate.right.limited;
            query_limited_[query] |= candidate.limited;
            start = std::chrono::steady_clock::now();
            assemble_candidate(candidate);
            query_timings_[query][2] +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            checked.push_back(std::move(candidate));
        }
    }
    select_candidate(checked);
    start = std::chrono::steady_clock::now();
    if (!checked.empty()) filter_selected_extension(checked.front());
    query_timings_[query][2] += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    query_status_[query] = checked.empty() ? "interrupted" : checked.front().reason;
    query_work_[query] = work;
    results_[query] = std::move(checked);
}

void ProtospacerWorkflow::consolidate() {
    std::map<std::string, std::string> targets, regions;
    std::set<std::string> strict_regions;
    for (auto &group : results_)
        for (auto &c : group)
            if (c.detected) {
                auto key = canonical(c.target);
                targets.emplace(key, "");
                if (c.extended) {
                    regions.emplace(canonical(c.sequence), "");
                    if (!c.allow_crossing)
                        strict_regions.insert(canonical(c.sequence));
                }
            }
    size_t number = 0;
    for (auto &item : targets)
        item.second = numbered("P", ++number);
    number = 0;
    for (auto &item : regions)
        item.second = numbered("E", ++number);
    // Consolidate strict-pass predictions first. A new recovery must not replace
    // an already accepted output with a longer, less certain sequence.
    // Never merge two partially overlapping paths into an invented sequence.
    std::vector<std::string> order;
    for (const auto &item : regions)
        order.push_back(item.first);
    std::sort(order.begin(), order.end(), [&](const auto &a, const auto &b) {
        if (strict_regions.count(a) != strict_regions.count(b))
            return strict_regions.count(a) > strict_regions.count(b);
        return a.size() != b.size() ? a.size() > b.size() : a < b;
    });
    std::vector<std::string> representatives;
    for (const auto &sequence : order) {
        bool contained = false;
        for (const auto &other : representatives)
            if (other.find(sequence) != std::string::npos ||
                other.find(reverse_complement(sequence)) != std::string::npos) {
                regions[sequence] = regions[other];
                contained = true;
                break;
            }
        if (!contained)
            representatives.push_back(sequence);
    }
    for (auto &group : results_)
        for (auto &c : group)
            if (c.detected) {
                c.target_id = targets.at(canonical(c.target));
                if (c.extended)
                    c.extension_id = regions.at(canonical(c.sequence));
            }
}

void ProtospacerWorkflow::write_outputs() {
    fs::path staging = options_.output / "staging";
    fs::create_directory(staging);
    auto open = [&](const std::string &name) {
        std::ofstream f(staging / name);
        f.exceptions(std::ios::badbit | std::ios::failbit);
        return f;
    };
    std::map<std::string, std::string> targets, candidate_targets, regions;
    for (const auto &group : results_)
        for (const auto &c : group)
            if (c.detected) {
                candidate_targets[c.target_id] = canonical(c.target);
                if (c.extended) {
                    targets[c.target_id] = canonical(c.target);
                    if (c.sequence.size() > regions[c.extension_id].size())
                        regions[c.extension_id] = canonical(c.sequence);
                }
            }
    {
        auto target_fasta = open("Protospacers.fasta"), extension_fasta = open("ContiguousPaths.fasta");
        auto candidate_fasta = open("CandidateProtospacers.fasta");
        for (const auto &[id, s] : candidate_targets)
            candidate_fasta << '>' << id << '\n' << s << '\n';
        for (const auto &[id, s] : targets)
            target_fasta << '>' << id << '\n' << s << '\n';
        for (const auto &[id, s] : regions)
            extension_fasta << '>' << id << '\n' << s << '\n';
    }
    if (!supplied_.empty()) {
        auto fasta = open("SuppliedProtospacers.fasta"), table = open("supplied_targets.tsv");
        table << "input_id\tsequence\tquery_id\n";
        std::map<std::string, size_t> index;
        for (size_t q = 0; q < queries_.size(); ++q)
            index[queries_[q].first] = q;
        for (const auto &[id, sequence] : supplied_) {
            fasta << '>' << id << '\n' << sequence << '\n';
            table << id << '\t' << sequence << '\t' << numbered("Q", index.at(canonical(sequence)) + 1)
                  << '\n';
        }
    }
    {
        auto input = open("sources.tsv"), queries = open("queries.tsv"), candidates = open("candidates.tsv"),
             rejected = open("unresolved_and_rejected.tsv");
        input << "source_id\tarray_id\toriginal_spacer\toriginal_"
                 "repeat\tspacer\trepeat\ttrim_left\ttrim_"
                 "right\n";
        for (const auto &s : sources_)
            input << s.id << '\t' << s.array << '\t' << s.original_spacer << '\t' << s.original_repeat << '\t'
                  << s.spacer << '\t' << s.repeat << '\t' << s.trim_left << '\t' << s.trim_right << '\n';
        queries << "query_id\tspacer\tsource_count\tstatus\twork\tsearch_limited\n";
        auto metrics = open("query_metrics.tsv");
        metrics << "query_id\tmapping_work\textension_work\tmapping_"
                   "seconds\textension_seconds\tfilter_seconds\n";
        candidates << "query_id\ttarget_id\tstatus\treason\tquery_begin\tquery_"
                      "end\tidentity\tquery_"
                      "coverage\tsequence\tleft_bp\tright_bp\tleft_stop\tright_"
                      "stop\tsearch_limited\t"
                      "left_deviation\tright_deviation\tleft_crosses_source_"
                      "array\tright_crosses_source_array\tidentified_before_"
                      "extension\tanchor_prefix_bp\tcoverage_blocks\tcoverage_block_min\tcoverage_block_max\n";
        rejected << "query_id\tstatus\treason\n";
        for (size_t i = 0; i < queries_.size(); ++i) {
            auto id = numbered("Q", i + 1);
            queries << id << '\t' << queries_[i].first << '\t' << queries_[i].second.size() << '\t'
                    << query_status_[i] << '\t' << query_work_[i] << '\t' << int(query_limited_[i]) << '\n';
            metrics << id << '\t' << mapping_work_[i] << '\t' << query_work_[i] - mapping_work_[i];
            for (double seconds : query_timings_[i])
                metrics << '\t' << seconds;
            metrics << '\n';
            if (results_[i].empty())
                rejected << id << "\tunresolved\t" << query_status_[i] << '\n';
            for (const auto &c : results_[i]) {
                candidates << id << '\t' << c.target_id << '\t' << c.status << '\t' << c.reason << '\t'
                           << c.begin << '\t' << c.end << '\t' << c.identity << '\t' << c.coverage << '\t'
                           << c.target << '\t' << c.left.bases.size() + c.anchor_prefix << '\t'
                           << c.right.bases.size() << '\t' << c.left.stop << '\t' << c.right.stop << '\t'
                           << c.limited << '\t' << c.left.deviation << '\t' << c.right.deviation << '\t'
                           << c.left.crosses_source << '\t' << c.right.crosses_source << '\t' << c.identified
                           << '\t' << c.anchor_prefix << '\t' << c.coverage_blocks << '\t' << c.coverage_min << '\t' << c.coverage_max << '\n';
                if (!c.extended)
                    rejected << id << '\t' << c.status << '\t' << c.reason << '\n';
            }
        }
    }
    {
        auto audit = open("withheld_extensions.tsv");
        audit << "query_id\ttarget_id\treason\tside\tsequence\tcoverage_blocks\tblock_min\tblock_max\tnodes_outward\tcount_positions\tcounts\n";
        for (const auto &group : results_)
            for (const auto &c : group) {
                if (c.withheld_sequence.empty()) continue;
                audit << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t' << c.reason << '\t'
                      << (c.withheld_left ? "left" : "right") << '\t' << c.withheld_sequence << '\t'
                      << c.coverage_blocks << '\t' << c.coverage_min << '\t' << c.coverage_max;
                auto values = [&](const auto &items) {
                    audit << '\t';
                    for (size_t i = 0; i < items.size(); ++i) {
                        if (i) audit << ',';
                        audit << items[i];
                    }
                };
                values(c.withheld.nodes);
                values(c.withheld.sample_positions);
                values(c.withheld.counts);
                audit << '\n';
            }
    }
    {
        auto links = open("associations.tsv"), candidate_links = open("candidate_associations.tsv"),
             evidence = open("coverage_evidence.tsv");
        const char *header = "source_id\tquery_id\ttarget_id\textension_id\tquery_begin\tquery_"
                 "end\tidentity\tquery_"
                 "coverage\ttarget_strand\tseed_begin\tseed_end\tsupported_"
                 "fragment_begin\tsupported_"
                 "fragment_end\tleft_bp\tright_bp\tleft_stop\tright_stop\tsearch_"
                 "limited\n";
        links << header;
        candidate_links << header;
        evidence << "query_id\ttarget_id\textension_id\tside\tadded_base_"
                    "position\tcount\tpath_capacity\n";
        for (const auto &group : results_)
            for (const auto &c : group)
                if (c.detected) {
                    size_t offset = 0, seed_begin = 0, seed_end = 0;
                    bool reverse = canonical(c.target) != c.target;
                    if (c.extended) {
                        const auto &region = regions.at(c.extension_id);
                        reverse = false;
                        offset = region.find(c.sequence);
                        if (offset == std::string::npos) {
                            reverse = true;
                            offset = region.find(reverse_complement(c.sequence));
                        }
                        if (offset == std::string::npos)
                            throw std::runtime_error("consolidation lost a supported fragment");
                        seed_begin =
                            offset + (reverse ? c.right.bases.size() : c.left.bases.size() + c.anchor_prefix);
                        seed_end = seed_begin + c.target.size();
                        if (region.substr(seed_begin, c.target.size()) !=
                            (reverse ? reverse_complement(c.target) : c.target))
                            throw std::runtime_error("extension target coordinates are inconsistent");
                    }
                    for (auto index : queries_[c.query].second) {
                        const auto &source = sources_[index];
                        bool source_flip = source.spacer != queries_[c.query].first;
                        size_t begin = source_flip ? source.spacer.size() - c.end : c.begin;
                        size_t end = source_flip ? source.spacer.size() - c.begin : c.end;
                        std::ostringstream row;
                        row << source.id << '\t' << numbered("Q", c.query + 1) << '\t' << c.target_id
                              << '\t' << c.extension_id << '\t' << begin + source.trim_left << '\t'
                              << end + source.trim_left << '\t' << c.identity << '\t' << c.coverage << '\t'
                              << ((reverse ^ source_flip) ? '-' : '+') << '\t';
                        if (c.extended)
                            row << seed_begin << '\t' << seed_end << '\t' << offset << '\t'
                                  << offset + c.sequence.size();
                        else
                            row << "\t\t\t";
                        row << '\t' << c.left.bases.size() + c.anchor_prefix << '\t' << c.right.bases.size()
                              << '\t' << c.left.stop << '\t' << c.right.stop << '\t' << c.limited << '\n';
                        candidate_links << row.str();
                        if (c.extended) links << row.str();
                    }
                    for (bool back : {true, false}) {
                        const auto &w = back ? c.left : c.right;
                        for (size_t j = 0; j < w.counts.size(); ++j)
                            evidence << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t' << c.extension_id << '\t'
                                     << (back ? "left" : "right") << '\t'
                                     << w.sample_positions[j] + (back ? c.anchor_prefix : 0) << '\t'
                                     << w.counts[j] << '\t' << w.capacities[j] << '\n';
                    }
                }
    }
    {
        auto contexts = open("target_context.tsv");
        contexts << "query_id\ttarget_id\tleft_outward\tright_outward\tleft_clear_"
                    "at\tright_clear_at\tleft_"
                    "stop_before_filter\tright_stop_before_filter\n";
        auto paths = open("graph_paths.tsv");
        paths << "query_id\ttarget_id\textension_id\tanchor_nodes\tleft_nodes_"
                 "outward\tright_nodes_outward\tanchor_prefix_bp\n";
        for (const auto &group : results_)
            for (const auto &c : group)
                if (c.detected) {
                    contexts << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t' << c.context_left
                             << '\t' << c.context_right << '\t' << c.clear_left << '\t' << c.clear_right
                             << '\t' << c.context_left_stop << '\t' << c.context_right_stop << '\n';
                    paths << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t' << c.extension_id;
                    for (const auto *nodes : {&c.anchor, &c.left.nodes, &c.right.nodes}) {
                        paths << '\t';
                        for (size_t i = 0; i < nodes->size(); ++i) {
                            if (i)
                                paths << ',';
                            paths << (*nodes)[i];
                        }
                    }
                    paths << '\t' << c.anchor_prefix << '\n';
                }
    }
    {
        auto scores = open("path_selection.tsv");
        scores << "query_id\tplacement\tphase\tselected\tside\tlength\tstop\traw_transition_energy\tadjacent_observations\tmean_energy\tbranch_loss\tbranch_observations\tambiguous\tcompeting_routes\tfirst_ambiguous_base\trunner_up_path_cost\tindependent_bases\tunphased_source_branch\tpath_cost\n";
        for (const auto &group : results_)
            for (const auto &c : group)
                for (bool back : {true, false}) {
                    const auto &w = back ? c.left : c.right;
                    scores << numbered("Q", c.query + 1) << '\t' << c.placement << '\t'
                           << (c.allow_crossing ? "crossing" : "strict") << '\t' << c.detected << '\t'
                           << (back ? "left" : "right") << '\t' << w.bases.size() + w.anchor_padding
                           << '\t' << w.stop << '\t' << w.transition_energy << '\t' << w.transitions
                           << '\t' << w.deviation << '\t' << w.branch_loss << '\t' << w.branch_observations
                           << '\t' << w.ambiguous << '\t' << w.competing_routes << '\t' << w.first_ambiguous_base
                           << '\t' << w.runner_up_score << '\t' << w.counts.size() << '\t'
                           << (w.crosses_source && w.pre_context_fork) << '\t' << walk_score(w) << '\n';
                }
    }
    for (const auto &f : fs::directory_iterator(staging))
        fs::rename(f.path(), options_.output / f.path().filename());
    fs::remove(staging);
}

void ProtospacerWorkflow::write_summary(const std::string &status, const std::string &error) const {
    if (!options_.extension_seeds.empty()) {
        std::ofstream out(options_.output / "summary.json");
        out << "{\"implementation\":\"protospacer-seeded-local-v1\",\"status\":" << json(status)
            << ",\"error\":" << json(error) << "}\n";
        return;
    }
    std::set<std::string> arrays, targets, candidate_targets, extensions;
    size_t detected = 0, extended = 0, limited = 0;
    for (const auto &s : sources_)
        arrays.insert(s.array);
    for (const auto &group : results_) {
        bool d = false, e = false, l = false;
        for (const auto &c : group) {
            d |= c.detected;
            e |= c.extended;
            l |= c.limited;
            if (c.detected)
                candidate_targets.insert(c.target_id);
            if (c.extended) {
                targets.insert(c.target_id);
                extensions.insert(c.extension_id);
            }
        }
        detected += d;
        extended += e;
        limited += l;
    }
    limited = std::count(query_limited_.begin(), query_limited_.end(), uint8_t(1));
    fs::path temporary = options_.output / "run_summary.json.tmp";
    std::ofstream out(temporary);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << "{\n  \"implementation\": \"protospacer-workflow-v9\",\n  \"status\": " << json(status)
        << ",\n  \"error\": " << json(error) << ",\n  \"arrays\": " << json(options_.arrays.string())
        << ",\n  \"graph\": " << json(options_.graph.string())
        << ",\n  \"graph_k\": " << (graph_ ? graph_->k() : 0) << ",\n  \"array_count\": " << arrays.size()
        << ",\n  \"spacer_observations\": " << sources_.size()
        << ",\n  \"distinct_spacers\": " << queries_.size() << ",\n  \"detected_queries\": " << detected
        << ",\n  \"protospacer_regions\": " << targets.size() << ",\n  \"extended_queries\": " << extended
        << ",\n  \"candidate_protospacer_regions\": " << candidate_targets.size()
        << ",\n  \"candidate_only_regions\": " << candidate_targets.size() - targets.size()
        << ",\n  \"protospacer_output_policy\": \"accepted_extensions_only\""
        << ",\n  \"unique_extensions\": " << extensions.size()
        << ",\n  \"queries_without_extension\": " << queries_.size() - extended
        << ",\n  \"limited_queries\": " << limited << ",\n  \"threads\": " << options_.threads
        << ",\n  \"beam_width\": " << options_.beam_width
        << ",\n  \"max_context_paths\": " << options_.context_width
        << ",\n  \"read_length\": " << options_.read_length
        << ",\n  \"min_external_bases\": " << options_.min_external
        << ",\n  \"max_flank_length\": " << options_.max_flank
        << ",\n  \"min_extension_flank\": " << options_.min_flank
        << ",\n  \"min_identity\": " << options_.identity
        << ",\n  \"min_query_coverage\": " << options_.coverage
        << ",\n  \"array_repeat_identity\": " << options_.repeat_identity
        << ",\n  \"boundary_support\": " << options_.boundary_support
        << ",\n  \"correct_spacer_boundaries\": " << (options_.correct_boundaries ? "true" : "false")
        << ",\n  \"save_search\": " << (options_.save_search ? "true" : "false")
        << ",\n  \"filter_extension_artifacts\": " << (options_.filter_artifacts ? "true" : "false")
        << ",\n  \"extension_complexity_window\": 64,\n  "
           "\"extension_complexity_stride\": 8"
        << ",\n  \"max_query_states\": " << options_.max_work
        << ",\n  \"max_placement_states\": " << options_.placement_work
        << ",\n  \"max_placements\": " << options_.max_placements
        << ",\n  \"branch_ratio\": " << options_.branch_ratio
        << ",\n  \"max_coverage_fold\": " << options_.max_coverage_fold
        << ",\n  \"coverage_filter_stage\": \"after_selection_no_fallback\""
        << ",\n  \"max_node_visits\": " << options_.max_visits
        << ",\n  \"extension_method\": "
           "\"conditional_branch_support_per_scored_edge\""
        << ",\n  \"context_bases\": " << options_.context_bases
        << ",\n  \"probabilities_of_correctness\": false,\n  "
           "\"reference_used_for_selection\": false,\n  "
           "\"source_array_crossing\": "
        << (options_.source_array_crossing ? "true" : "false")
        << ",\n  \"supplied_targets\": " << supplied_.size()
        << ",\n  \"target_mode\": " << json(options_.protospacers.empty() ? "discovery" : "supplied")
        << ",\n  \"source_crossing_establishes_detection\": false,\n  "
           "\"timings_seconds\": {";
    bool first = true;
    for (const auto &[name, seconds] : timings_) {
        if (!first)
            out << ',';
        out << '\n' << "    " << json(name) << ": " << seconds;
        first = false;
    }
    out << "\n  }\n}\n";
    out.close();
    fs::rename(temporary, options_.output / "run_summary.json");
    if (status == "complete")
        std::cout << "Candidate protospacer regions (diagnostic): " << candidate_targets.size()
                  << "\nProtospacer regions with an accepted extension: " << targets.size()
                  << "\nSaved extension regions: " << extensions.size()
                  << "\nQueries without an extension: " << queries_.size() - extended << '\n';
}
int ProtospacerWorkflow::run() {
    auto stage_start = std::chrono::steady_clock::now();
    auto finish = [&](const std::string &name) {
        auto now = std::chrono::steady_clock::now();
        timings_[name] = std::chrono::duration<double>(now - stage_start).count();
        stage_start = now;
    };
    write_summary("running");
    std::cout << "Loading existing graph (read-only)...\n" << std::flush;
    graph_ = std::make_unique<SDBG>();
    graph_->LoadFromFile(options_.graph.string().c_str());
    finish("load_graph");
    options_.context_bases = graph_->k() + 1;
    verify_index();
    finish("verify_index");
    load_sources();
    correct_boundaries();
    prepare_source_nodes();
    load_supplied_targets();
    group_queries();
    if (!options_.extension_seeds.empty()) {
        finish("prepare_spacers");
        return run_seeded();
    }
    if (options_.save_search)
        fs::create_directory(options_.output / "search_paths");
    finish("prepare_spacers");
    std::cout << sources_.size() << " source observations; " << queries_.size()
              << " distinct corrected spacers\n"
              << std::flush;
    write_summary("running");
    std::atomic<size_t> completed{0};
    std::vector<std::string> errors(queries_.size());
    omp_set_num_threads(options_.threads);
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < queries_.size(); ++i) {
        if (interrupted_) {
            query_status_[i] = "interrupted";
            continue;
        }
        try {
            process_query(i);
        } catch (const std::exception &e) {
            errors[i] = e.what();
        }
        size_t done = ++completed;
        if (done % 10 == 0 || done == queries_.size()) {
#pragma omp critical(workflow_progress)
            std::cout << "Processed " << done << '/' << queries_.size() << " spacers\n" << std::flush;
        }
    }
    for (size_t i = 0; i < errors.size(); ++i)
        if (!errors[i].empty())
            throw std::runtime_error(numbered("Q", i + 1) + ": " + errors[i]);
    finish("detect_extend_filter");
    consolidate();
    finish("consolidate");
    write_outputs();
    finish("write_outputs");
    write_summary(interrupted_ ? "interrupted" : "complete");
    std::cout << (interrupted_ ? "Workflow interrupted; partial outputs saved. Counts: "
                               : "Workflow complete. Final counts: ")
              << (options_.output / "run_summary.json") << '\n';
    return interrupted_ ? 128 + interrupted_ : 0;
}
double ProtospacerWorkflow::count_loss(double count, double mean, const std::string &method,
                                       double dispersion) {
    mean = std::max(mean, 1e-9);
    if (method == "huber" || method == "flow" || method == "lookahead") {
        double z = std::abs(std::log((count + .5) / (mean + .5)));
        return z <= 1. ? .5 * z * z : z - .5;
    }
    double term = count > 0. ? count * std::log(count / mean) : 0.;
    if (method == "negbin")
        return std::max(
            0., 2. * (term - (count + dispersion) * std::log((count + dispersion) / (mean + dispersion))));
    return std::max(0., 2. * (term - count + mean));
}
double ProtospacerWorkflow::sequence_similarity(const std::string &a, const std::string &b, double minimum) {
    size_t maximum = std::max(a.size(), b.size());
    if (!maximum)
        return 1.;
    size_t band = size_t(maximum * (1. - minimum) + 1e-9);
    if (maximum - std::min(a.size(), b.size()) > band)
        return 0.;
    std::vector<size_t> previous(b.size() + 1, maximum + 1), next(previous);
    for (size_t j = 0; j <= std::min(band, b.size()); ++j)
        previous[j] = j;
    for (size_t i = 1; i <= a.size(); ++i) {
        std::fill(next.begin(), next.end(), maximum + 1);
        if (i <= band)
            next[0] = i;
        for (size_t j = std::max(size_t(1), i > band ? i - band : size_t(1));
             j <= std::min(b.size(), i + band); ++j)
            next[j] = std::min({previous[j] + 1, next[j - 1] + 1, previous[j - 1] + (a[i - 1] != b[j - 1])});
        previous.swap(next);
    }
    return 1. - double(previous.back()) / maximum;
}
bool ProtospacerWorkflow::same_graph_region(const std::string &a, const std::string &b, size_t k) {
    // Direct-to-representative grouping, never transitive chaining. A shared
    // oriented k-mer is the same de Bruijn graph node; require a long identical
    // node run as well as high identity across the overlapping sequence.
    if (a.size() < k || b.size() < k)
        return a == b;
    std::map<std::string, std::vector<int64_t>> index;
    for (size_t i = 0; i + k <= a.size(); ++i)
        index[a.substr(i, k)].push_back(i);
    for (const auto &s : {b, reverse_complement(b)}) {
        std::set<int64_t> offsets;
        for (size_t j = 0; j + k <= s.size(); ++j) {
            auto it = index.find(s.substr(j, k));
            if (it != index.end())
                for (auto i : it->second)
                    offsets.insert(i - int64_t(j));
        }
        for (auto offset : offsets) {
            size_t ab = std::max(int64_t(0), offset), bb = std::max(int64_t(0), -offset);
            size_t overlap = std::min(a.size() - ab, s.size() - bb);
            if (10 * overlap < 9 * std::min(a.size(), s.size()))
                continue;
            size_t identical = 0, longest = 0, run = 0;
            for (size_t i = 0; i < overlap; ++i) {
                bool equal = a[ab + i] == s[bb + i];
                identical += equal;
                run = equal ? run + 1 : 0;
                longest = std::max(longest, run);
            }
            if (100 * identical >= 95 * overlap && longest >= k &&
                10 * (longest - k + 1) >= 9 * (std::min(a.size(), s.size()) - k + 1))
                return true;
        }
    }
    return false;
}
LocalSegment ProtospacerWorkflow::probe_segment(uint64_t from, GraphArc arc, size_t query, bool backwards,
                                                size_t maximum_edges, size_t &work, size_t limit) const {
    LocalSegment segment;
    std::set<uint64_t> visited{from};
    while (work < limit && segment.arcs.size() < maximum_edges && !interrupted_) {
        ++work;
        if (!visited.insert(arc.node).second)
            break;
        segment.arcs.push_back(arc);
        segment.bases += added_base(arc, backwards);
        if (!query_repeat_nodes_[query].count(from) && !query_repeat_nodes_[query].count(arc.node)) {
            segment.depth += arc.count;
            ++segment.depth_samples;
        }
        ++work;
        auto outgoing = neighbors(arc.node, backwards);
        ++work;
        segment.terminal_indegree = neighbors(arc.node, !backwards).size();
        segment.dead_end = outgoing.empty();
        if (outgoing.size() != 1 || segment.terminal_indegree != 1) {
            segment.complete = true;
            break;
        }
        from = arc.node;
        arc = outgoing.front();
    }
    if (segment.depth_samples)
        segment.depth /= segment.depth_samples;
    return segment;
}
void ProtospacerWorkflow::score_seeded_walk(ExtensionWalk &walk) const {
    const size_t k = graph_->k();
    const size_t span = std::max(k + 1, options_.read_length > k ? options_.read_length - k : size_t(1));
    if (walk.counts.size() < walk.coverage_seen) {
        walk.coverage_baseline = walk.coverage_loss = 0.;
        walk.coverage_seen = 0;
    }
    if (walk.counts.size() < span)
        return;
    if (!walk.coverage_seen) {
        auto values = std::vector<double>(walk.counts.begin(), walk.counts.begin() + span);
        std::sort(values.begin(), values.end());
        walk.coverage_baseline = .5 * (values[(span - 1) / 2] + values[span / 2]);
        // Calibrate on the initial read-span; assess count consistency afterward.
        walk.coverage_seen = span;
    }
    while (walk.coverage_seen < walk.counts.size())
        walk.coverage_loss += count_loss(walk.counts[walk.coverage_seen++], walk.coverage_baseline,
                                         options_.coverage_method, options_.nb_dispersion);
}
double ProtospacerWorkflow::seeded_walk_score(const ExtensionWalk &walk) const {
    const size_t k = graph_->k();
    const size_t span = std::max(k + 1, options_.read_length > k ? options_.read_length - k : size_t(1));
    // A flat count-1 path must still pay for a weak branch decision. Overlapping
    // k-mer counts share reads; give a read-span one unit of coverage weight.
    // This composite ranking score is not a calibrated path probability.
    double support = walk.counts.empty() ? 0. : std::log(2.) * walk.branch_loss / walk.counts.size();
    size_t tested = walk.coverage_seen > span ? walk.coverage_seen - span : 0;
    double scale = (options_.coverage_method == "huber" || options_.coverage_method == "flow" || options_.coverage_method == "lookahead") ? 1. : .5;
    return support + (tested ? scale * walk.coverage_loss / (tested * double(span)) : 0.);
}

std::pair<int, double> ProtospacerWorkflow::flow_choice(double incoming,
                                                       const std::vector<double> &outgoing) {
    // Excess flow: incoming support minus flow that may leave by other arms.
    // Real counts are not conserved; additionally cap support by the chosen arm.
    double total = std::accumulate(outgoing.begin(), outgoing.end(), 0.);
    int chosen = -1;
    double support = 0.;
    for (size_t i = 0; i < outgoing.size(); ++i) {
        double next = std::min(outgoing[i], incoming - (total - outgoing[i]));
        if (next <= 1e-9) continue;
        if (chosen >= 0) return {-2, 0.}; // Multiple supported continuations: do not guess.
        chosen = int(i);
        support = next;
    }
    return {chosen, support};
}

double ProtospacerWorkflow::support_cost(const TargetCandidate &candidate) {
    // Worst-side evidence: adding a poor flank must not improve the score by
    // averaging it with a well-supported flank or simply making a longer record.
    double cost = 0.;
    for (const auto *w : {&candidate.left, &candidate.right}) {
        if (w->bases.empty()) continue;
        double side = w->uses_lookahead ? w->selection_metric : w->seed_flow > 0.
            ? 1. - std::clamp(w->remaining_flow / w->seed_flow, 0., 1.)
            : (w->counts.empty() ? 1. : w->branch_loss / w->counts.size());
        cost = std::max(cost, side);
    }
    return cost;
}

bool ProtospacerWorkflow::supported_representative(const TargetCandidate &a, const TargetCandidate &b) {
    if (support_cost(a) != support_cost(b)) return support_cost(a) < support_cost(b);
    if (a.limited != b.limited) return !a.limited;
    bool aa = a.left.ambiguous || a.right.ambiguous, ba = b.left.ambiguous || b.right.ambiguous;
    if (aa != ba) return !aa;
    double av = std::max(a.left.deviation, a.right.deviation);
    double bv = std::max(b.left.deviation, b.right.deviation);
    if (av != bv) return av < bv;
    if (a.sequence.size() != b.sequence.size()) return a.sequence.size() > b.sequence.size();
    return std::tie(a.sequence, a.query) < std::tie(b.sequence, b.query);
}

ExtensionWalk ProtospacerWorkflow::extend_flow(size_t query, const TargetCandidate &candidate,
                                              bool backwards, size_t &work, size_t limit) const {
    ExtensionWalk walk;
    walk.anchor_padding = backwards ? candidate.anchor_prefix : 0;
    walk.nodes = {backwards ? candidate.anchor.front() : candidate.anchor.back()};
    // The complete anchor constrains the initial flow. A one-node anchor has
    // only node-level evidence, so use the smaller incident directional total.
    double initial = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < candidate.anchor.size(); ++i) {
        ++work;
        for (const auto &arc : neighbors(candidate.anchor[i - 1], false))
            if (arc.node == candidate.anchor[i]) initial = std::min(initial, arc.count);
    }
    if (!std::isfinite(initial)) {
        initial = std::numeric_limits<double>::infinity();
        for (bool direction : {false, true}) {
            ++work;
            double total = 0.;
            for (const auto &arc : neighbors(walk.nodes.back(), direction)) total += arc.count;
            initial = std::min(initial, total);
        }
    }
    walk.seed_flow = walk.remaining_flow = initial;
    const size_t k = graph_->k();
    while (!interrupted_ && work < limit && walk.bases.size() + walk.anchor_padding < options_.max_flank) {
        auto node = walk.nodes.back();
        ++work;
        auto arcs = neighbors(node, backwards);
        if (arcs.empty()) { walk.stop = "dead_end"; break; }
        size_t remaining = options_.max_flank - walk.bases.size() - walk.anchor_padding;
        std::vector<LocalSegment> arms;
        std::vector<size_t> retained;
        for (size_t i = 0; i < arcs.size(); ++i) {
            arms.push_back(probe_segment(node, arcs[i], query, backwards, remaining, work, limit));
            bool tip = arcs.size() > 1 && arms.back().complete && arms.back().dead_end &&
                       arms.back().arcs.size() < 2 * k;
            if (tip) ++walk.local_tips;
            else retained.push_back(i);
        }
        if (interrupted_ || work >= limit) {
            walk.stop = interrupted_ ? "interrupted" : "query_work_limit";
            walk.limited = true;
            break;
        }
        if (retained.empty()) { walk.stop = "local_short_tips"; break; }
        auto unmasked = [&](size_t i) {
            return !arms[i].arcs.empty() && arms[i].depth_samples == arms[i].arcs.size();
        };
        // Resolve only an actually reconverging, clearly deeper sequence bubble.
        // Tied variants remain unresolved. Never infer a tip from a probe bound.
        if (retained.size() > 1 && std::all_of(retained.begin(), retained.end(), [&](size_t i) {
                return arms[i].complete && !arms[i].dead_end && unmasked(i);
            })) {
            auto strongest = *std::max_element(retained.begin(), retained.end(), [&](size_t a, size_t b) {
                return arms[a].depth < arms[b].depth;
            });
            double second = 0.;
            bool bubble = true, similar = true;
            auto start_label = label(node);
            if (backwards) std::reverse(start_label.begin(), start_label.end());
            for (auto i : retained) {
                if (i != strongest) second = std::max(second, arms[i].depth);
                bubble &= arms[i].arcs.back().node == arms[strongest].arcs.back().node &&
                          arms[i].terminal_indegree == retained.size();
                similar &= arms[i].arcs.size() <= k + 2 ||
                    (arms[i].arcs.size() <= size_t(std::lround(20. * k / .98)) &&
                     sequence_similarity(start_label + arms[strongest].bases,
                                         start_label + arms[i].bases, .98) >= .98);
            }
            if (bubble && similar && arms[strongest].depth > options_.branch_ratio * second) {
                walk.local_bubbles += retained.size() - 1;
                retained = {strongest};
            }
        }
        std::vector<double> weights;
        for (auto i : retained) {
            // Unitig mean depth avoids treating each overlapping k-mer as an
            // independent observation. Masked source-repeat means are not used.
            weights.push_back(unmasked(i) ? arms[i].depth : arcs[i].count);
        }
        auto [choice, next] = flow_choice(walk.remaining_flow, weights);
        walk.flow_decisions.push_back({walk.bases.size(), node, walk.remaining_flow, next, weights, choice});
        if (choice < 0) {
            walk.stop = choice == -2 ? "unresolved_flow_branch" : "insufficient_seed_flow";
            // The ambiguity is at the excluded next edge, not inside this prefix.
            break;
        }
        auto selected = retained[size_t(choice)];
        const auto &arm = arms[selected];
        if (arm.arcs.empty()) { walk.stop = "probe_cycle"; break; }
        walk.remaining_flow = next; // Merges and source-repeat depth never replenish it.
        double total = 0.;
        for (const auto &arc : arcs) total += arc.count;
        bool continuing = true;
        for (size_t j = 0; j < arm.arcs.size(); ++j) {
            if (work++ >= limit) { walk.limited = true; walk.stop = "query_work_limit"; continuing = false; break; }
            double branch = j == 0 && arcs.size() > 1 && arcs[selected].count > 0.
                ? -std::log2(arcs[selected].count / total) : 0.;
            continuing = advance_walk(walk, arm.arcs[j], query, candidate, backwards, branch);
            update_profile(walk);
            score_seeded_walk(walk);
            if (!continuing) break;
        }
        if (!continuing) break;
    }
    if (interrupted_) { walk.stop = "interrupted"; walk.limited = true; }
    else if (work >= limit) { walk.stop = "query_work_limit"; walk.limited = true; }
    // An artifact or source barrier may shorten the selected arm. Remove flow
    // decisions outside the retained prefix and restore its actual support.
    trim_walk(walk, walk.bases.size(), walk.stop);
    score_seeded_walk(walk);
    return walk;
}

std::vector<double> ProtospacerWorkflow::future_depths(size_t query, uint64_t from, const GraphArc &first,
                                                     bool backwards, size_t distance, size_t &work,
                                                     size_t limit) const {
    std::vector<double> depths;
    std::deque<std::tuple<uint64_t, GraphArc, size_t>> pending{{from, first, 0}};
    std::set<uint64_t> seen;
    while (!pending.empty() && seen.size() < 256 && limit - std::min(work, limit) >= 3 && !interrupted_) {
        auto [node, arc, used] = pending.front();
        pending.pop_front();
        if (used >= distance || !seen.insert(arc.edge).second) continue;
        // probe_segment accounts three graph operations per edge. Reserve the
        // final two so even a probe ending mid-unitig obeys this smaller budget.
        auto segment = probe_segment(node, arc, query, backwards, distance - used, work, limit - 2);
        if (segment.arcs.empty()) continue;
        if (segment.depth_samples == segment.arcs.size() && segment.depth_samples >= size_t(graph_->k() + 1))
            depths.push_back(segment.depth);
        used += segment.arcs.size();
        if (used < distance && work < limit) {
            ++work;
            auto end = segment.arcs.back().node;
            for (const auto &next : neighbors(end, backwards)) pending.emplace_back(end, next, used);
        }
    }
    return depths;
}

void ProtospacerWorkflow::reserve_branch_families(std::vector<ExtensionWalk> &walks, size_t width) {
    if (walks.size() <= width) return;
    // Input order is the caller's reference-free preference. Round robin gives
    // each of the first two unresolved branch choices a protected quota.
    std::map<std::vector<uint64_t>, size_t> index;
    std::vector<std::vector<size_t>> families;
    for (size_t i = 0; i < walks.size(); ++i) {
        auto [it, inserted] = index.emplace(walks[i].branch_family, families.size());
        if (inserted) families.emplace_back();
        families[it->second].push_back(i);
    }
    std::vector<ExtensionWalk> kept;
    for (size_t round = 0; kept.size() < width; ++round) {
        bool any = false;
        for (const auto &family : families) {
            if (round < family.size() && kept.size() < width) {
                kept.push_back(std::move(walks[family[round]]));
                any = true;
            }
        }
        if (!any) break;
    }
    walks = std::move(kept);
}

bool ProtospacerWorkflow::compatible_flanks(const TargetCandidate &seed, const ExtensionWalk &left,
                                           const ExtensionWalk &right) const {
    if (left.coverage_baseline > 0. && right.coverage_baseline > 0. && options_.max_coverage_fold > 0. &&
        std::max(left.coverage_baseline, right.coverage_baseline) > options_.max_coverage_fold *
            std::min(left.coverage_baseline, right.coverage_baseline)) return false;
    std::map<uint64_t, size_t> visits;
    for (auto node : seed.anchor) ++visits[node];
    for (const auto *walk : {&left, &right})
        for (size_t i = 1; i < walk->nodes.size(); ++i) {
            auto node = walk->nodes[i];
            if (++visits[node] > (query_source_nodes_[seed.query].count(node) ? 1 : options_.max_visits)) return false;
        }
    return true;
}

ExtensionWalk ProtospacerWorkflow::extend_seeded(size_t query, const TargetCandidate &candidate,
                                                 bool backwards, size_t &work, size_t limit,
                                                 std::vector<ExtensionWalk> *alternatives) const {
    const bool lookahead = options_.coverage_method == "lookahead";
    ExtensionWalk seed;
    seed.uses_lookahead = lookahead;
    seed.anchor_padding = backwards ? candidate.anchor_prefix : 0;
    seed.nodes = {backwards ? candidate.anchor.front() : candidate.anchor.back()};
    std::vector<ExtensionWalk> frontier{seed}, finished;
    std::map<uint64_t, std::vector<GraphArc>> adjacency;
    std::map<uint64_t, std::vector<LocalSegment>> segments;
    std::map<std::pair<uint64_t, size_t>, std::vector<double>> future;
    // Optional branch probes must not consume the budget needed to extend the
    // candidates themselves. Each flank uses at most one eighth of the query
    // budget for these probes; each new probe uses at most 512 graph operations.
    size_t future_work = 0;
    const size_t future_budget = options_.max_work / 8;
    bool limited = false;
    const size_t k = graph_->k();
    const size_t span = std::max(k + 1, options_.read_length > k ? options_.read_length - k : size_t(1));
    auto continuity = [&](const ExtensionWalk &w) {
        double low = 0., high = 0.;
        size_t blocks = 0;
        for (size_t b = 0; b + span <= w.counts.size(); b += span) {
            auto v = std::vector<double>(w.counts.begin() + b, w.counts.begin() + b + span);
            std::sort(v.begin(), v.end());
            double m = .5 * (v[(span - 1) / 2] + v[span / 2]);
            low = blocks ? std::min(low, m) : m;
            high = std::max(high, m);
            ++blocks;
        }
        return options_.max_coverage_fold == 0. || blocks < 2 ||
               (low > 0. && high <= options_.max_coverage_fold * low);
    };
    auto score = [&](const ExtensionWalk &w) { return seeded_walk_score(w); };
    auto better = [&](const auto &a, const auto &b) {
        if (a.external != b.external)
            return a.external > b.external;
        if (score(a) != score(b))
            return score(a) < score(b);
        if (!a.coverage_seen && !b.coverage_seen && a.deviation != b.deviation)
            return a.deviation < b.deviation;
        if (a.bases.size() != b.bases.size())
            return a.bases.size() > b.bases.size();
        return a.nodes < b.nodes;
    };
    auto keep = [&](std::vector<ExtensionWalk> &pool, size_t width, bool diverse) {
        std::sort(pool.begin(), pool.end(), better);
        if (pool.size() <= width)
            return;
        limited = true;
        if (lookahead) {
            reserve_branch_families(pool, width);
            return;
        }
        if (diverse) {
            std::set<std::tuple<uint64_t, bool, bool>> endpoints;
            std::vector<ExtensionWalk> unique, extra;
            for (auto &w : pool)
                (endpoints.emplace(w.nodes.back(), w.external, w.crosses_source).second ? unique : extra)
                    .push_back(std::move(w));
            pool = std::move(unique);
            for (auto &w : extra)
                pool.push_back(std::move(w));
        }
        pool.resize(width);
    };
    while (!frontier.empty() && !interrupted_) {
        std::vector<ExtensionWalk> pending, mature;
        for (auto &current : frontier) {
            if (work >= limit) {
                current.stop = "query_work_limit";
                current.limited = limited = true;
                finished.push_back(std::move(current));
                continue;
            }
            auto node = current.nodes.back();
            if (!adjacency.count(node)) {
                ++work;
                adjacency[node] = neighbors(node, backwards);
            }
            const auto &arcs = adjacency.at(node);
            if (arcs.empty()) {
                current.stop = "dead_end";
                finished.push_back(std::move(current));
                continue;
            }
            std::vector<size_t> retained;
            for (size_t i = 0; i < arcs.size(); ++i)
                retained.push_back(i);
            size_t tips = 0, bubbles = 0, weak = 0;
            if (arcs.size() > 1) {
                if (lookahead && current.counts.size() >= std::max(options_.min_external, 2 * span)) {
                    auto checkpoint = current;
                    checkpoint.stop = "supported_prefix";
                    finished.push_back(std::move(checkpoint));
                }
                if (!segments.count(node)) {
                    std::vector<LocalSegment> probed;
                    const size_t remaining =
                        options_.max_flank - current.bases.size() - current.anchor_padding;
                    for (const auto &arc : arcs)
                        probed.push_back(probe_segment(node, arc, query, backwards, remaining, work, limit));
                    segments.emplace(node, std::move(probed));
                }
                const auto &arms = segments.at(node);
                // Same 2*k short-tip rule. A probe cut by our length/work bound
                // is not a dead end and must never be classified as a tip.
                retained.erase(std::remove_if(retained.begin(), retained.end(),
                                              [&](size_t i) {
                                                  bool tip = arms[i].complete && arms[i].dead_end &&
                                                             arms[i].arcs.size() < 2 * size_t(graph_->k());
                                                  tips += tip;
                                                  return tip;
                                              }),
                               retained.end());
                bool depths_known = std::all_of(retained.begin(), retained.end(), [&](size_t i) {
                    return arms[i].complete && arms[i].depth_samples &&
                           arms[i].depth_samples == arms[i].arcs.size();
                });
                if (retained.size() > 1 && depths_known) {
                    uint64_t end = arms[retained.front()].arcs.back().node;
                    bool bubble = std::all_of(retained.begin(), retained.end(), [&](size_t i) {
                        return !arms[i].dead_end && arms[i].arcs.back().node == end &&
                               arms[i].terminal_indegree == retained.size();
                    });
                    if (bubble) {
                        auto strongest =
                            *std::min_element(retained.begin(), retained.end(), [&](size_t a, size_t b) {
                                return arms[a].depth != arms[b].depth ? arms[a].depth > arms[b].depth
                                                                      : arcs[a].node < arcs[b].node;
                            });
                        bool simple = std::all_of(retained.begin(), retained.end(), [&](size_t i) {
                            return arms[i].arcs.size() <= size_t(graph_->k() + 2);
                        });
                        auto start_label = label(node);
                        if (backwards)
                            std::reverse(start_label.begin(), start_label.end());
                        bool similar = std::all_of(retained.begin(), retained.end(), [&](size_t i) {
                            return arms[i].arcs.size() <= size_t(std::lround(20. * graph_->k() / .98)) &&
                                   sequence_similarity(start_label + arms[strongest].bases,
                                                       start_label + arms[i].bases, .98) >= .98;
                        });
                        double second = 0.;
                        for (auto j : retained) if (j != strongest) second = std::max(second, arms[j].depth);
                        if ((simple || similar) && (!lookahead || arms[strongest].depth > options_.branch_ratio * second)) {
                            bubbles = retained.size() - 1;
                            retained = {strongest};
                        }
                    }
                    if (retained.size() > 1) {
                        double total = 0.;
                        for (auto i : retained)
                            total += arms[i].depth;
                        retained.erase(std::remove_if(retained.begin(), retained.end(),
                                                      [&](size_t i) {
                            // A source-exit bridge can have count 1 even when
                            // its eventual external continuation has count 40.
                            // Do not prune it before external calibration, or
                            // discard a rare arm consistent with this target's
                            // depth merely because a pooled neighbor is abundant.
                            bool below_local_depth = current.coverage_seen &&
                                arms[i].depth * options_.max_coverage_fold < current.coverage_baseline;
                            bool low = !lookahead && options_.max_coverage_fold > 0. && below_local_depth &&
                                arms[i].depth <= .1 * total;
                                                          weak += low;
                                                          return low;
                                                      }),
                                       retained.end());
                    }
                }
            }
            if (retained.empty()) {
                current.stop = "local_short_tips";
                current.local_tips += tips;
                finished.push_back(std::move(current));
            }
            std::vector<double> weights(arcs.size(), 0.);
            double weight_total = 0.;
            for (size_t i : retained) {
                weights[i] = arcs[i].count;
                if (lookahead && retained.size() > 1) {
                    size_t distance = std::min(2 * span, options_.max_flank - current.bases.size() - current.anchor_padding);
                    auto key = std::make_pair(arcs[i].edge, distance);
                    if (!future.count(key)) {
                        size_t before = work;
                        size_t budget = std::min(size_t(512), future_budget - std::min(future_work, future_budget));
                        future[key] = future_depths(query, node, arcs[i], backwards, distance,
                                                    work, std::min(limit, work + budget));
                        future_work += work - before;
                    }
                    const auto &depths = future.at(key);
                    if (current.coverage_seen) {
                        double penalty = count_loss(arcs[i].count, current.coverage_baseline, "huber", 10.);
                        if (!depths.empty()) {
                            penalty = std::numeric_limits<double>::infinity();
                            for (auto depth : depths) penalty = std::min(penalty, count_loss(depth, current.coverage_baseline, "huber", 10.));
                        }
                        weights[i] = std::exp(-penalty);
                    } else if (!depths.empty()) weights[i] = *std::max_element(depths.begin(), depths.end());
                }
                weight_total += weights[i];
            }
            for (size_t i : retained) {
                auto child = current;
                if (lookahead && retained.size() > 1 && child.branch_family.size() < 2)
                    child.branch_family.push_back(arcs[i].node);
                child.local_tips += tips;
                child.local_bubbles += bubbles;
                child.local_weak_links += weak;
                if (work++ >= limit) {
                    child.stop = "query_work_limit";
                    child.limited = limited = true;
                    finished.push_back(std::move(child));
                    continue;
                }
                if (retained.size() > 1 && !child.external) {
                    if (!child.pre_context_fork)
                        child.first_context_fork = child.bases.size();
                    child.pre_context_fork = true;
                }
                double total_count = 0.;
                for (const auto &arc : arcs) total_count += arc.count;
                double branch_cost = arcs.size() > 1 && arcs[i].count > 0. && total_count > 0.
                    ? -std::log2(arcs[i].count / total_count) : 0.;
                if (lookahead && retained.size() > 1 && weight_total > 0.)
                    branch_cost = -std::log2(std::max(1e-12, weights[i] / weight_total));
                bool continuing = advance_walk(child, arcs[i], query, candidate, backwards, branch_cost);
                update_profile(child);
                score_seeded_walk(child);
                if (lookahead && continuing && child.counts.size() >= 2 * span && child.counts.size() % span == 0) {
                    auto checkpoint = child;
                    checkpoint.stop = "supported_prefix";
                    finished.push_back(std::move(checkpoint));
                }
                if (!continuing)
                    finished.push_back(std::move(child));
                else
                    (child.coverage_seen ? mature : pending).push_back(std::move(child));
            }
        }
        keep(pending, options_.context_width, true);
        keep(mature, options_.beam_width, true);
        frontier = std::move(pending);
        for (auto &w : mature)
            frontier.push_back(std::move(w));
        // Keep complete walks separately: a short, perfectly flat dead end
        // must not evict a completed 750-base candidate during the search.
        std::vector<ExtensionWalk> full, partial;
        for (auto &w : finished)
            if (w.external && w.counts.size() >= options_.context_bases && w.flow > 0. && continuity(w))
                (w.bases.size() + w.anchor_padding >= options_.max_flank ? full : partial)
                    .push_back(std::move(w));
        keep(full, options_.beam_width * 4, false);
        keep(partial, options_.beam_width * 4, false);
        finished = std::move(full);
        for (auto &w : partial)
            finished.push_back(std::move(w));
    }
    if (interrupted_) {
        for (auto &w : frontier) {
            w.stop = "interrupted";
            w.limited = true;
            finished.push_back(std::move(w));
        }
    }
    // Apply the same read-span continuity gate to each candidate before final
    // selection. Source repeat counts are already excluded by advance_walk.
    finished.erase(std::remove_if(finished.begin(), finished.end(),
                                  [&](const auto &w) {
                                      return !w.external || w.counts.size() < options_.context_bases ||
                                             w.flow <= 0. || !continuity(w);
                                  }),
                   finished.end());
    if (finished.empty()) {
        seed.stop = "no_coverage_supported_extension";
        seed.limited = limited;
        return seed;
    }
    bool full = std::any_of(finished.begin(), finished.end(), [&](const auto &w) {
        return w.bases.size() + w.anchor_padding >= options_.max_flank;
    });
    if (full && !lookahead)
        finished.erase(std::remove_if(finished.begin(), finished.end(),
                                      [&](const auto &w) {
                                          return w.bases.size() + w.anchor_padding < options_.max_flank;
                                      }),
                       finished.end());
    std::sort(finished.begin(), finished.end(), better);
    if (lookahead) {
        // Preserve actual terminal routes as well as independently scored prefixes.
        std::vector<ExtensionWalk> terminal, prefixes;
        for (auto &w : finished) {
            w.selection_metric = score(w);
            (w.stop == "supported_prefix" ? prefixes : terminal).push_back(w);
        }
        std::stable_sort(terminal.begin(), terminal.end(), [&](const auto &a, const auto &b) {
            bool af = a.bases.size() + a.anchor_padding >= options_.max_flank;
            bool bf = b.bases.size() + b.anchor_padding >= options_.max_flank;
            return af != bf ? af > bf : better(a, b);
        });
        reserve_branch_families(terminal, 16);
        reserve_branch_families(prefixes, 16);
        if (alternatives) {
            std::set<std::string> seen;
            for (auto *pool : {&terminal, &prefixes})
                for (auto &w : *pool) if (seen.insert(w.bases).second) alternatives->push_back(w);
            std::sort(alternatives->begin(), alternatives->end(), better);
            for (auto &w : *alternatives) {
                w.limited |= limited;
                for (const auto &other : *alternatives) {
                    if (same_route(w, other)) continue;
                    ++w.competing_routes;
                    if (w.competing_routes == 1 || score(other) < w.runner_up_score)
                        w.runner_up_score = score(other);
                    if (score(other) <= score(w) * options_.branch_ratio + 1e-12)
                        w.ambiguous = true;
                }
            }
        }
    }
    auto best = std::move(finished.front());
    best.limited |= limited;
    if (finished.size() > 1) {
        best.runner_up_score = score(finished[1]);
        best.ambiguous = best.runner_up_score <= score(best) * options_.branch_ratio + 1e-12;
        best.competing_routes = finished.size() - 1;
    }
    return best;
}
std::vector<TargetCandidate> ProtospacerWorkflow::load_extension_seeds() const {
    auto table = [&](const std::string &name) {
        std::ifstream input(options_.extension_seeds / name);
        if (!input)
            throw std::runtime_error("cannot open seed table " + name);
        auto split = [](const std::string &line) {
            std::vector<std::string> cells;
            size_t first = 0;
            do {
                auto end = line.find('\t', first);
                cells.push_back(line.substr(first, end - first));
                if (end == std::string::npos)
                    break;
                first = end + 1;
            } while (true);
            return cells;
        };
        std::string line;
        std::getline(input, line);
        auto header = split(line);
        std::vector<std::map<std::string, std::string>> rows;
        while (std::getline(input, line)) {
            auto cells = split(line);
            if (cells.size() != header.size())
                throw std::runtime_error("malformed seed table " + name);
            std::map<std::string, std::string> row;
            for (size_t i = 0; i < cells.size(); ++i)
                row[header[i]] = cells[i];
            rows.push_back(std::move(row));
        }
        return rows;
    };
    std::map<std::string, size_t> query_index;
    for (size_t q = 0; q < queries_.size(); ++q)
        query_index[numbered("Q", q + 1)] = q;
    auto queries = table("queries.tsv");
    if (queries.size() != queries_.size())
        throw std::runtime_error("seed input query set differs from arrays");
    for (const auto &r : queries)
        if (!query_index.count(r.at("query_id")) ||
            queries_[query_index.at(r.at("query_id"))].first != r.at("spacer"))
            throw std::runtime_error("seed query sequence differs from corrected array input");
    std::map<std::string, std::pair<std::string, std::string>> source_index;
    for (const auto &s : sources_)
        source_index[s.id] = {s.spacer, s.repeat};
    auto source_rows = table("sources.tsv");
    if (source_rows.size() != sources_.size())
        throw std::runtime_error("seed source set differs from arrays");
    for (const auto &r : source_rows)
        if (!source_index.count(r.at("source_id")) ||
            source_index.at(r.at("source_id")) != std::make_pair(r.at("spacer"), r.at("repeat")))
            throw std::runtime_error("seed source boundaries differ from array input");
    std::map<std::string, std::vector<uint64_t>> anchors;
    for (const auto &r : table("graph_paths.tsv")) {
        auto &nodes = anchors[r.at("query_id")];
        if (!nodes.empty())
            throw std::runtime_error("duplicate saved query placement");
        std::istringstream input(r.at("anchor_nodes"));
        for (std::string word; std::getline(input, word, ',');) {
            auto node = std::stoull(word);
            if (node >= graph_->size())
                throw std::runtime_error("seed node outside graph");
            nodes.push_back(node);
        }
    }
    std::vector<TargetCandidate> seeds;
    std::set<size_t> seen;
    for (const auto &r : table("candidates.tsv")) {
        if (r.at("status") != "detected" && r.at("status") != "extended")
            continue;
        TargetCandidate c;
        c.query = query_index.at(r.at("query_id"));
        if (!seen.insert(c.query).second)
            throw std::runtime_error("multiple saved starts for one query are not "
                                     "supported by this seed format");
        c.anchor = anchors.at(r.at("query_id"));
        c.target = r.at("sequence");
        c.target_id = r.at("target_id");
        c.begin = std::stoull(r.at("query_begin"));
        c.end = std::stoull(r.at("query_end"));
        c.anchor_prefix = std::stoull(r.at("anchor_prefix_bp"));
        c.identity = std::stod(r.at("identity"));
        c.coverage = std::stod(r.at("query_coverage"));
        c.identified = c.detected = true;
        c.allow_crossing = options_.source_array_crossing;
        if (c.anchor.empty() || c.begin >= c.end || c.end > queries_[c.query].first.size() ||
            c.end - c.begin != c.target.size() || !dna(c.target))
            throw std::runtime_error("invalid saved target coordinates");
        auto spelled = spell(c.anchor);
        if (spelled.size() != c.anchor_prefix + c.target.size() ||
            spelled.substr(c.anchor_prefix) != c.target)
            throw std::runtime_error("saved seed does not spell its target on this graph");
        seeds.push_back(std::move(c));
    }
    if (seeds.size() != anchors.size())
        throw std::runtime_error("saved candidate and path counts differ");
    return seeds;
}
int ProtospacerWorkflow::run_seeded() {
    auto seeds = load_extension_seeds();
    std::map<std::string, std::string> targets;
    for (const auto &c : seeds) {
        auto sequence = canonical(c.target);
        if (targets.count(c.target_id) && targets.at(c.target_id) != sequence)
            throw std::runtime_error("target ID has conflicting sequences");
        targets[c.target_id] = sequence;
    }
    {
        std::ofstream out(options_.output / "protospacers.fasta");
        out.exceptions(std::ios::badbit | std::ios::failbit);
        for (const auto &[id, sequence] : targets)
            out << '>' << id << '\n' << sequence << '\n';
    }
    std::vector<std::string> methods = options_.coverage_method == "all"
                                           ? std::vector<std::string>{"poisson", "negbin", "huber"}
                                           : std::vector<std::string>{options_.coverage_method};
    std::ofstream evidence(options_.output / "extensions.tsv");
    evidence.exceptions(std::ios::badbit | std::ios::failbit);
    evidence << "method\tquery_id\ttarget_id\taccepted\tviral_contig_"
                "id\trelation_to_representative\treason\tsources\t"
                "anchor_first\tanchor_last\tseed_begin\tseed_end\tleft_bp\tright_"
                "bp\tleft_stop\tright_stop\t"
                "left_baseline\tright_baseline\tmean_loss\tsearch_"
                "limited\tambiguous\twork\t"
                "tips_removed\tbubbles_removed\tweak_links_"
                "removed\tsequence\tleft_counts\tright_counts\t"
                "left_scored_positions\tright_scored_positions\tleft_branch_bits\tright_branch_bits\tselection_cost\tleft_seed_flow\tright_seed_flow\tleft_remaining_flow\tright_remaining_flow\tsupport_cost\tleft_flow_decisions\tright_flow_decisions\treported\trepresentative_rank\n";
    std::vector<std::string> summaries;
    std::unique_ptr<gzFile_s, int (*)(gzFile)> candidate_audit(nullptr, &gzclose);
    if (options_.coverage_method == "lookahead") {
        candidate_audit.reset(gzopen((options_.output / "candidate_paths.tsv.gz").c_str(), "wb"));
        if (!candidate_audit) throw std::runtime_error("cannot open candidate audit");
        gzputs(candidate_audit.get(), "query_id\ttarget_id\tside\tindex\tscore\tsequence\tseed_begin\tseed_end\tstop\tlimited\n");
    }
    auto write_metadata = [&](const std::string &status) {
        std::ofstream out(options_.output / "summary.json");
        out.exceptions(std::ios::badbit | std::ios::failbit);
        out << "{\n\"implementation\":\"protospacer-seeded-local-v5\",\"status\":" << json(status)
            << ",\n\"graph\":" << json(options_.graph.string()) << ",\"graph_k\":" << graph_->k()
            << ",\n\"seed_directory\":" << json(options_.extension_seeds.string())
            << ",\"input_starts\":" << seeds.size() << ",\"input_target_sequences\":" << targets.size()
            << ",\"original_corrected_spacers\":" << queries_.size()
            << ",\n\"max_flank_per_side\":" << options_.max_flank
            << ",\"min_contig_length\":" << options_.min_contig_length
            << ",\"beam_width\":" << options_.beam_width << ",\"context_width\":" << options_.context_width
            << ",\"max_work_per_start\":" << options_.max_work << ",\"read_length\":" << options_.read_length
            << ",\"nb_dispersion\":" << options_.nb_dispersion
            << ",\"huber_delta\":1,\"huber_pseudocount\":0.5"
            << ",\"max_coverage_fold\":" << options_.max_coverage_fold
            << ",\"source_array_crossing\":" << (options_.source_array_crossing ? "true" : "false")
            << ",\"threads\":" << options_.threads
            << ",\n\"score_policy\":\"mean branch surprisal plus read-span-weighted coverage loss; initial window used for calibration only\""
            << ",\n\"lookahead_policy\":\"two read-spans, 256 unitig-entry edges, 512 operations per probe; at most one eighth of query work per flank for probing; first-two-fork family quotas; no hard weak-link pruning; independent prefix selection; top 50 regions\""
            << ",\n\"flow_policy\":\"unique positive excess-flow continuation; cap by arm depth; no budget replenishment; stop before unresolved split\""
            << ",\n\"representative_policy\":\"lowest worst-flank support cost, then uncertainty, then length; no reference\""
            << ",\n\"weak_link_policy\":\"only after external calibration, and only below the local baseline by max_coverage_fold as well as <=10% neighbor depth\""
            << ",\n\"reference_used_for_selection\":false,\"whole_graph_assembly\":"
               "false,\"global_graph_modified\":false,"
               "\n\"baseline_policy\":\"per-path median of first read-span of "
               "external edge counts; then frozen\","
               "\n\"grouping_policy\":\"direct representative; >=90% identical "
               "contiguous graph-node run and >=95% overlapping sequence "
               "identity\","
               "\n\"methods\":[";
        for (size_t i = 0; i < summaries.size(); ++i) {
            if (i)
                out << ',';
            out << summaries[i];
        }
        out << "],\n\"initialization_seconds\":{";
        bool first = true;
        for (const auto &[name, time] : timings_) {
            if (!first)
                out << ',';
            first = false;
            out << json(name) << ':' << time;
        }
        out << "}\n}\n";
    };
    write_metadata("running");
    std::cout << "Seed-local extension: " << seeds.size() << " saved starts; " << targets.size()
              << " target sequences. No graph assembly.\n"
              << std::flush;
    for (const auto &method : methods) {
        options_.coverage_method = method;
        auto start = std::chrono::steady_clock::now();
        auto predictions = seeds;
        std::vector<size_t> work(seeds.size());
        std::vector<std::string> errors(seeds.size());
        std::atomic<size_t> done{0};
        omp_set_num_threads(options_.threads);
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < predictions.size(); ++i) {
            auto &c = predictions[i];
            if (interrupted_) {
                c.reason = "interrupted";
                c.limited = true;
                continue;
            }
            try {
                std::vector<ExtensionWalk> left_pool, right_pool;
                if (method == "lookahead") {
                    c.left = extend_seeded(c.query, c, true, work[i], options_.max_work / 2, &left_pool);
                    c.right = extend_seeded(c.query, c, false, work[i], options_.max_work, &right_pool);
                    ExtensionWalk empty_left, empty_right;
                    empty_left.nodes = {c.anchor.front()}; empty_right.nodes = {c.anchor.back()};
                    empty_left.stop = empty_right.stop = "side_not_selected";
                    empty_left.uses_lookahead = empty_right.uses_lookahead = true;
                    left_pool.push_back(empty_left); right_pool.push_back(empty_right);
                    std::vector<std::tuple<double, size_t, size_t, size_t>> pairs;
                    size_t span = std::max(size_t(graph_->k() + 1), options_.read_length > size_t(graph_->k()) ? options_.read_length - graph_->k() : size_t(1));
                    for (size_t l = 0; l < left_pool.size(); ++l) for (size_t r = 0; r < right_pool.size(); ++r) {
                        size_t length = left_pool[l].bases.size() + right_pool[r].bases.size() + c.target.size() + c.anchor_prefix;
                        if (length < options_.min_contig_length || left_pool[l].counts.size() + right_pool[r].counts.size() < std::max(options_.min_external, 2 * span)) continue;
                        pairs.emplace_back(std::max(left_pool[l].selection_metric, right_pool[r].selection_metric),
                            size_t(2 * options_.max_flank + c.target.size() - length), l, r);
                    }
                    std::sort(pairs.begin(), pairs.end());
                    for (const auto &[score, length_order, l, r] : pairs)
                        if (compatible_flanks(c, left_pool[l], right_pool[r])) {
                            c.left = left_pool[l]; c.right = right_pool[r]; break;
                        }
                } else if (method == "flow") {
                    c.left = extend_flow(c.query, c, true, work[i], options_.max_work / 2);
                    c.right = extend_flow(c.query, c, false, work[i], options_.max_work);
                } else {
                    c.left = extend_seeded(c.query, c, true, work[i], options_.max_work / 2);
                    c.right = extend_seeded(c.query, c, false, work[i], options_.max_work);
                }
                c.limited = c.left.limited || c.right.limited;
                // Do not join incompatible depths or introduce a forbidden
                // source-array revisit by concatenating independently explored sides.
                auto prefer_left = [&] {
                    if (method == "flow") {
                        double l = c.left.seed_flow > 0. ? c.left.remaining_flow / c.left.seed_flow : 0.;
                        double r = c.right.seed_flow > 0. ? c.right.remaining_flow / c.right.seed_flow : 0.;
                        if (l != r) return l > r;
                    }
                    if (c.left.bases.size() != c.right.bases.size())
                        return c.left.bases.size() > c.right.bases.size();
                    return seeded_walk_score(c.left) <= seeded_walk_score(c.right);
                };
                std::map<uint64_t, size_t> visits;
                for (auto n : c.anchor)
                    ++visits[n];
                bool overlap = false;
                for (const auto *w : {&c.left, &c.right})
                    for (size_t j = 1; j < w->nodes.size(); ++j) {
                        auto n = w->nodes[j];
                        size_t allowed = query_source_nodes_[c.query].count(n) ? 1 : options_.max_visits;
                        if (++visits[n] > allowed)
                            overlap = true;
                    }
                bool incompatible = c.left.coverage_baseline > 0. && c.right.coverage_baseline > 0. &&
                                    options_.max_coverage_fold > 0. &&
                                    std::max(c.left.coverage_baseline, c.right.coverage_baseline) >
                                        options_.max_coverage_fold *
                                            std::min(c.left.coverage_baseline, c.right.coverage_baseline);
                if (overlap || incompatible) {
                    auto &remove = prefer_left() ? c.right : c.left;
                    trim_walk(remove, 0, overlap ? "opposite_flank_revisit" : "incompatible_flank_depths");
                    score_seeded_walk(remove);
                    if (method == "lookahead") remove.selection_metric = seeded_walk_score(remove);
                }
                auto left = c.left.bases;
                std::reverse(left.begin(), left.end());
                auto padding = c.anchor_prefix ? label(c.anchor.front()).substr(0, c.anchor_prefix) : "";
                c.sequence = left + padding + c.target + c.right.bases;
                std::vector<uint64_t> path(c.left.nodes.rbegin(), c.left.nodes.rend());
                if (path.empty())
                    path.push_back(c.anchor.front());
                path.insert(path.end(), c.anchor.begin() + 1, c.anchor.end());
                if (!c.right.nodes.empty())
                    path.insert(path.end(), c.right.nodes.begin() + 1, c.right.nodes.end());
                if (spell(path) != c.sequence)
                    throw std::runtime_error("seeded extension is not the saved graph walk");
                if (c.left.bases.size() + c.anchor_prefix > options_.max_flank ||
                    c.right.bases.size() > options_.max_flank)
                    throw std::runtime_error("flank length exceeds requested maximum");
                const size_t k = graph_->k();
                const size_t span =
                    std::max(k + 1, options_.read_length > k ? options_.read_length - k : size_t(1));
                c.extended =
                    !interrupted_ && c.sequence.size() >= options_.min_contig_length &&
                    c.left.counts.size() + c.right.counts.size() >= std::max(options_.min_external, 2 * span);
                c.reason = interrupted_ ? "interrupted"
                           : c.extended ? "accepted"
                           : c.sequence.size() < options_.min_contig_length
                               ? "shorter_than_minimum_contig"
                               : "insufficient_external_coverage";
                if (method == "lookahead") {
                    std::ostringstream audit;
                    for (bool backwards : {true, false}) {
                        const auto &pool = backwards ? left_pool : right_pool;
                        std::set<std::string> seen;
                        for (size_t a = 0; a < pool.size(); ++a) {
                            const auto &w = pool[a];
                            auto bases = w.bases;
                            if (backwards) std::reverse(bases.begin(), bases.end());
                            auto sequence = backwards ? bases + spell(c.anchor) : spell(c.anchor) + bases;
                            if (sequence.size() < options_.min_contig_length || !seen.insert(sequence).second) continue;
                            std::vector<uint64_t> path;
                            if (backwards) {
                                path.assign(w.nodes.rbegin(), w.nodes.rend());
                                path.insert(path.end(), c.anchor.begin() + 1, c.anchor.end());
                            } else {
                                path = c.anchor;
                                if (!w.nodes.empty()) path.insert(path.end(), w.nodes.begin() + 1, w.nodes.end());
                            }
                            if (spell(path) != sequence) throw std::runtime_error("candidate audit walk is not connected");
                            size_t begin = (backwards ? w.bases.size() : 0) + c.anchor_prefix;
                            audit << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t'
                                  << (backwards ? "left" : "right") << '\t' << a + 1 << '\t'
                                  << w.selection_metric << '\t' << sequence << '\t' << begin << '\t'
                                  << begin + c.target.size() << '\t' << w.stop << '\t' << w.limited << '\n';
                        }
                    }
                    auto text = audit.str();
                    int written = 0;
#pragma omp critical(candidate_audit_write)
                    written = gzwrite(candidate_audit.get(), text.data(), unsigned(text.size()));
                    if (written != int(text.size())) throw std::runtime_error("cannot write candidate audit");
                }
            } catch (const std::exception &e) {
                errors[i] = e.what();
            }
            auto completed = ++done;
            if (completed % 10 == 0 || completed == predictions.size()) {
#pragma omp critical(seeded_progress)
                std::cout << method << ": " << completed << '/' << predictions.size() << " starts\n"
                          << std::flush;
            }
        }
        for (size_t i = 0; i < errors.size(); ++i)
            if (!errors[i].empty())
                throw std::runtime_error(numbered("Q", predictions[i].query + 1) + ": " + errors[i]);
        double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
        const size_t span = std::max(size_t(graph_->k() + 1), options_.read_length > size_t(graph_->k())
            ? options_.read_length - graph_->k() : size_t(1));
        auto loss = [&](const TargetCandidate &c) {
            size_t n = (c.left.coverage_seen > span ? c.left.coverage_seen - span : 0)
                + (c.right.coverage_seen > span ? c.right.coverage_seen - span : 0);
            return n ? (c.left.coverage_loss + c.right.coverage_loss) / n : 0.;
        };
        auto selection_cost = [&](const TargetCandidate &c) {
            size_t n = c.left.counts.size() + c.right.counts.size();
            return n ? (seeded_walk_score(c.left) * c.left.counts.size()
                + seeded_walk_score(c.right) * c.right.counts.size()) / n : 0.;
        };
        std::vector<size_t> ordered;
        for (size_t i = 0; i < predictions.size(); ++i)
            if (predictions[i].extended)
                ordered.push_back(i);
        std::sort(ordered.begin(), ordered.end(), [&](size_t a, size_t b) {
            const auto &x = predictions[a];
            const auto &y = predictions[b];
            return supported_representative(x, y);
        });
        std::map<std::string, std::set<size_t>> representative_nodes;
        std::vector<size_t> representatives;
        std::vector<std::string> relation(predictions.size(), "unextended");
        size_t k = graph_->k();
        for (size_t i : ordered) {
            auto &c = predictions[i];
            std::set<size_t> possible;
            // A >=90% shared run in a >=500 bp contig necessarily includes a
            // 64-spaced query seed. Representatives index every k-mer.
            for (size_t p = 0; p + k <= c.sequence.size(); p += 64) {
                auto it = representative_nodes.find(canonical(c.sequence.substr(p, k)));
                if (it != representative_nodes.end())
                    possible.insert(it->second.begin(), it->second.end());
            }
            for (size_t r : possible) {
                const auto &rep = predictions[representatives[r]];
                if (same_graph_region(rep.sequence, c.sequence, k)) {
                    c.extension_id = rep.extension_id;
                    relation[i] =
                        rep.sequence.find(c.sequence) != std::string::npos ||
                                rep.sequence.find(reverse_complement(c.sequence)) != std::string::npos
                            ? "contained"
                            : "same_graph_region_variant";
                    break;
                }
            }
            if (c.extension_id.empty()) {
                size_t r = representatives.size();
                representatives.push_back(i);
                c.extension_id = numbered("V", r + 1);
                relation[i] = "representative";
                for (size_t p = 0; p + k <= c.sequence.size(); ++p)
                    representative_nodes[canonical(c.sequence.substr(p, k))].insert(r);
            }
        }
        std::ofstream fasta(options_.output / ("viral_contigs_" + method + ".fasta"));
        fasta.exceptions(std::ios::badbit | std::ios::failbit);
        std::map<std::string, size_t> ranks;
        for (size_t rank = 0; rank < representatives.size(); ++rank) {
            size_t i = representatives[rank];
            ranks[predictions[i].extension_id] = rank + 1;
            if (method != "lookahead" || rank < 50)
                fasta << '>' << predictions[i].extension_id << '\n' << predictions[i].sequence << '\n';
        }
        size_t accepted = 0, paired = 0, limited = 0;
        std::set<std::string> accepted_targets;
        for (size_t i = 0; i < predictions.size(); ++i) {
            const auto &c = predictions[i];
            accepted += c.extended;
            limited += c.limited;
            paired += c.extended && !c.left.bases.empty() && !c.right.bases.empty();
            if (c.extended)
                accepted_targets.insert(c.target_id);
            evidence << method << '\t' << numbered("Q", c.query + 1) << '\t' << c.target_id << '\t'
                     << c.extended << '\t' << c.extension_id << '\t' << relation[i] << '\t' << c.reason
                     << '\t';
            bool first = true;
            for (auto s : queries_[c.query].second) {
                if (!first)
                    evidence << ';';
                first = false;
                evidence << sources_[s].id;
            }
            size_t begin = c.left.bases.size() + c.anchor_prefix;
            evidence << '\t' << c.anchor.front() << '\t' << c.anchor.back() << '\t' << begin << '\t'
                     << begin + c.target.size() << '\t' << begin << '\t' << c.right.bases.size() << '\t'
                     << c.left.stop << '\t' << c.right.stop << '\t' << c.left.coverage_baseline << '\t'
                     << c.right.coverage_baseline << '\t' << loss(c) << '\t' << c.limited << '\t'
                     << (c.left.ambiguous || c.right.ambiguous) << '\t' << work[i] << '\t'
                     << c.left.local_tips + c.right.local_tips << '\t'
                     << c.left.local_bubbles + c.right.local_bubbles << '\t'
                     << c.left.local_weak_links + c.right.local_weak_links << '\t' << c.sequence;
            auto values = [&](const auto &v) {
                evidence << '\t';
                for (size_t j = 0; j < v.size(); ++j) {
                    if (j)
                        evidence << ',';
                    evidence << v[j];
                }
            };
            values(c.left.all_counts);
            values(c.right.all_counts);
            values(c.left.sample_positions);
            values(c.right.sample_positions);
            evidence << '\t' << c.left.branch_loss << '\t' << c.right.branch_loss << '\t' << selection_cost(c);
            evidence << '\t' << c.left.seed_flow << '\t' << c.right.seed_flow
                     << '\t' << c.left.remaining_flow << '\t' << c.right.remaining_flow
                     << '\t' << support_cost(c);
            for (const auto *w : {&c.left, &c.right}) {
                evidence << '\t';
                bool first_decision = true;
                for (const auto &d : w->flow_decisions) {
                    if (!first_decision) evidence << ';';
                    first_decision = false;
                    evidence << d.position << ':' << d.node << ':' << d.incoming << ':';
                    for (size_t j = 0; j < d.outgoing.size(); ++j) {
                        if (j) evidence << ',';
                        evidence << d.outgoing[j];
                    }
                    evidence << ':' << d.chosen << ':' << d.remaining;
                }
            }
            size_t rank = c.extended ? ranks.at(c.extension_id) : 0;
            evidence << '\t' << (rank > 0 && (method != "lookahead" || rank <= 50)) << '\t' << rank << '\n';
        }
        evidence.flush();
        double grouping = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::ostringstream result;
        result << "{\"method\":" << json(method) << ",\"processed_starts\":" << done.load()
               << ",\"accepted_starts\":" << accepted
               << ",\"accepted_target_sequences\":" << accepted_targets.size()
               << ",\"viral_contig_representatives\":" << representatives.size()
               << ",\"reported_representatives\":" << (method == "lookahead" ? std::min(size_t(50), representatives.size()) : representatives.size())
               << ",\"paired_flank_starts\":" << paired << ",\"limited_starts\":" << limited
               << ",\"extension_seconds\":" << seconds << ",\"grouping_output_seconds\":" << grouping << '}';
        summaries.push_back(result.str());
        write_metadata(interrupted_ ? "interrupted" : "running");
        std::cout << method << ": " << accepted << " accepted starts; " << representatives.size()
                  << " representative regions\n"
                  << std::flush;
        if (interrupted_)
            break;
    }
    if (candidate_audit) {
        if (gzclose(candidate_audit.release()) != Z_OK) throw std::runtime_error("cannot finish candidate audit");
    }
    write_metadata(interrupted_ ? "interrupted" : "complete");
    return interrupted_ ? 128 + interrupted_ : 0;
}
void ProtospacerWorkflow::self_test() {
    auto require = [](bool condition, const std::string &reason) {
        if (!condition)
            throw std::runtime_error("self-test failed: " + reason);
    };
    require(flow_choice(10., {9., 1.}) == std::make_pair(-2, 0.), "two positive flow continuations stop");
    require(flow_choice(9., {9., 1.}) == std::make_pair(0, 8.), "unique positive excess-flow branch");
    require(flow_choice(1., {9., 1.}) == std::make_pair(-1, 0.), "zero excess flow cannot extend");
    require(flow_choice(4., {100.}) == std::make_pair(0, 4.), "shared repeat does not inflate seed flow");
    require(flow_choice(4., {2.}) == std::make_pair(0, 2.), "real-data flow capped at arm depth");
    std::vector<ExtensionWalk> family_test(20);
    for (auto &walk : family_test) walk.branch_family = {1, 2};
    family_test.back().branch_family = {1, 3};
    reserve_branch_families(family_test, 4);
    require(family_test.size() == 4 && family_test[1].branch_family == std::vector<uint64_t>({1, 3}),
            "minority early branch retains a beam slot despite worse rank");
    TargetCandidate short_supported, long_weak;
    short_supported.sequence.assign(600, 'A'); long_weak.sequence.assign(1500, 'A');
    short_supported.left.bases.assign(570, 'A'); long_weak.left.bases.assign(1470, 'A');
    short_supported.left.seed_flow = long_weak.left.seed_flow = 20.;
    short_supported.left.remaining_flow = 18.; long_weak.left.remaining_flow = 2.;
    require(supported_representative(short_supported, long_weak), "support outranks representative length");
    for (const auto &method : {"poisson", "negbin", "huber"}) {
        require(count_loss(10, 10, method, 10) == 0., "coverage loss minimum");
        require(count_loss(3, 10, method, 10) > 0. && count_loss(30, 10, method, 10) > 0.,
                "both multiplicity deviations have positive loss");
    }
    require(std::abs(count_loss(0, 10, "poisson", 10) - 20.) < 1e-10, "zero-count Poisson limit");
    require(std::abs(count_loss(0, 10, "negbin", 10) - 20. * std::log(2.)) < 1e-10,
            "negative-binomial deviance formula");
    require(count_loss(30, 10, "negbin", 10) < count_loss(30, 10, "poisson", 10),
            "overdispersed model tolerates extra count variation");
    require(std::abs(count_loss(10.5 * std::exp(2.) - .5, 10, "huber", 10) - 1.5) < 1e-10,
            "Huber linear tail on log multiplicity ratio");
    require(sequence_similarity(std::string(100, 'A'), std::string(98, 'A') + "CC", .98) >= .98,
            "complex bubble similarity inclusive boundary");
    require(sequence_similarity(std::string(100, 'A'), std::string(97, 'A') + "CCC", .98) < .98,
            "complex bubble rejects excessive differences");
    std::string core;
    uint64_t random_state = 37;
    for (size_t i = 0; i < 733; ++i) {
        random_state = random_state * 6364136223846793005ULL + 1;
        core += "ACGT"[(random_state >> 32) & 3];
    }
    require(same_graph_region("AAAA" + core, "CCCC" + core, 23), "shared long graph region grouped");
    require(same_graph_region(core, reverse_complement(core), 23), "reverse-complement region grouped");
    require(!same_graph_region(core, core.substr(0, 100) + std::string(633, 'A'), 23),
            "short shared repeat does not merge different regions");
    require(reverse_complement("AACG") == "CGTT", "reverse complement");
    require(canonical("AACG") == canonical("CGTT"), "strand canonicalization");
    require(repeat_context("ACGTACGT", "", "ACGTACGT", .85) == 1, "source context rejection");
    require(repeat_context("ACG", "", "ACGTACGT", .85) == 0, "incomplete context");
    require(repeat_context("TTTTTTTTTT", "", "ACGTACGT", .85) == 2, "different context");
    require(repeat_context("GGACGTACGT", "GG", "ACGTACGT", .85) == 1, "unmatched spacer remainder");
    require(triplet_entropy("AAAAAAAAAA") == 0, "uninformative sequence");
    require(json("a\n\"b") == "\"a\\u000a\\\"b\"", "JSON escaping");
    require(!informative("TGTGGGGGGGGGGGGGGGTGTGAAA"), "homopolymer target rejection");
    ProtospacerWorkflow workflow(WorkflowOptions{});
    ExtensionWalk trimmed_flow;
    trimmed_flow.bases = "ACGTAC";
    trimmed_flow.nodes.resize(7);
    trimmed_flow.seed_flow = 20.;
    trimmed_flow.remaining_flow = 10.;
    trimmed_flow.flow_decisions = {{0, 1, 20., 18., {18., 2.}, 0}, {3, 2, 18., 10., {10., 8.}, 0}};
    workflow.trim_walk(trimmed_flow, 3, "test_prefix");
    require(trimmed_flow.flow_decisions.size() == 1 && trimmed_flow.remaining_flow == 18.,
            "trim restores support of retained prefix, excluding the next edge");
    workflow.options_.max_coverage_fold = 4.;
    auto coverage_case = [](const std::vector<double> &counts) {
        TargetCandidate c;
        c.target = "ACGT";
        c.detected = c.extended = true;
        c.context_right = "original discovery evidence";
        c.right.counts = c.right.all_counts = counts;
        c.right.capacities.assign(counts.size(), 1.);
        c.right.sample_positions.resize(counts.size());
        std::iota(c.right.sample_positions.begin(), c.right.sample_positions.end(), size_t(1));
        c.right.bases.assign(counts.size(), 'A');
        c.right.nodes.resize(counts.size() + 1);
        c.sequence = c.target + c.right.bases;
        return c;
    };
    std::vector<double> stable_counts(254, 10.);
    stable_counts[160] = 1.;
    stable_counts[161] = 200.;
    auto stable_coverage = coverage_case(stable_counts);
    workflow.filter_selected_extension(stable_coverage);
    require(stable_coverage.extended, "isolated count dips and spikes do not reject a stable path");
    auto jump_counts = std::vector<double>(254, 10.);
    std::fill(jump_counts.begin() + 127, jump_counts.end(), 41.);
    auto jump_coverage = coverage_case(jump_counts);
    workflow.filter_selected_extension(jump_coverage);
    require(jump_coverage.detected && !jump_coverage.extended && jump_coverage.target == "ACGT" &&
                jump_coverage.context_right == "original discovery evidence" &&
                jump_coverage.right.bases.empty() && jump_coverage.withheld.counts == jump_counts &&
                jump_coverage.withheld_sequence.size() == 258 &&
                jump_coverage.reason == "coverage_discontinuity",
            "coverage rejection preserves the target, discovery context and "
            "complete withheld-route evidence");
    std::fill(jump_counts.begin() + 127, jump_counts.end(), 40.);
    auto boundary_coverage = coverage_case(jump_counts);
    workflow.filter_selected_extension(boundary_coverage);
    require(boundary_coverage.extended, "coverage fold bound is inclusive");
    std::fill(jump_counts.begin(), jump_counts.begin() + 127, 200.);
    auto drop_coverage = coverage_case(jump_counts);
    workflow.filter_selected_extension(drop_coverage);
    require(!drop_coverage.extended, "persistent abundance drops are detected as well as increases");
    auto masked_coverage = coverage_case(jump_counts);
    masked_coverage.right.counts.assign(254, 40.);
    workflow.filter_selected_extension(masked_coverage);
    require(masked_coverage.extended, "masked source counts do not calibrate external coverage");
    auto incomplete_counts = std::vector<double>(254, 10.);
    incomplete_counts.insert(incomplete_counts.end(), 126, 1000.);
    auto incomplete = coverage_case(incomplete_counts);
    workflow.filter_selected_extension(incomplete);
    require(incomplete.extended && incomplete.coverage_blocks == 2,
            "partial trailing block does not establish a persistent change");
    workflow.options_.max_coverage_fold = 0.;
    auto disabled_coverage = coverage_case(jump_counts);
    workflow.filter_selected_extension(disabled_coverage);
    require(disabled_coverage.extended, "coverage filter can be disabled explicitly");
    workflow.options_.max_coverage_fold = 2.;
    const std::string adapter = "AGATCGGAAGAGCGTCGT";
    require(adapter_position("TT" + adapter) == 2 && adapter_position(reverse_complement(adapter)) == 0,
            "adapter core recognition on both strands");
    auto mutated_adapter = adapter;
    mutated_adapter.back() = 'A';
    require(adapter_position("TT" + mutated_adapter) == 2,
            "one substitution in a complete technical adapter remains detectable");
    require(adapter_position(adapter.substr(1)) == std::string::npos,
            "short adapter fragment alone is not rejected");
    const std::string shared_adapter = "AGATCGGAAGAGC";
    require(adapter_position(shared_adapter + "GTAGAGTAGGG") == 0,
            "shared adapter prefix survives different downstream adapter variants");
    require(adapter_position(reverse_complement(shared_adapter) + "ACGT") == 0 &&
                adapter_position(std::string(shared_adapter.rbegin(), shared_adapter.rend()), true) == 0,
            "shared adapter recognized in both strand and traversal orientations");
    auto altered_prefix = shared_adapter;
    altered_prefix[6] = 'A';
    require(adapter_position(altered_prefix) == std::string::npos &&
                adapter_position(shared_adapter.substr(1)) == std::string::npos,
            "short prefix requires all thirteen bases without mismatches");
    TargetCandidate test_target;
    test_target.target = "ACGTTAGCCATGATCGTACCTAGGTCATGC";
    ExtensionWalk technical;
    technical.bases = std::string(200, 'A') + adapter;
    technical.nodes.resize(technical.bases.size() + 1);
    require(workflow.extension_barrier(technical, test_target, false) && technical.bases.size() == 200 &&
                technical.nodes.size() == 201 && technical.stop == "adapter_sequence",
            "adapter branch retains only its preceding graph prefix");
    auto reverse_adapter = std::string(adapter.rbegin(), adapter.rend());
    technical.bases = std::string(200, 'A') + reverse_adapter;
    technical.nodes.resize(technical.bases.size() + 1);
    require(workflow.extension_barrier(technical, test_target, true) && technical.bases.size() == 200,
            "backward traversal recognizes the reversed motif");
    test_target.target += adapter.substr(0, 10);
    technical.bases = adapter.substr(10);
    technical.nodes.resize(technical.bases.size() + 1);
    require(workflow.extension_barrier(technical, test_target, false) && technical.bases.empty(),
            "adapter crossing the seed boundary is stopped");
    technical.bases = std::string(264, 'C');
    technical.nodes.resize(265);
    require(workflow.extension_barrier(technical, test_target, false) && technical.bases.size() == 200 &&
                technical.stop == "low_complexity_extension",
            "low-complexity continuation stops before its uninformative window");
    auto walk = [](size_t length, double flow, uint64_t branch) {
        ExtensionWalk w;
        w.external = true;
        w.clear_at = 1;
        w.flow = flow;
        w.counts.assign(length, flow);
        w.capacities.assign(length, flow);
        w.stop = "dead_end";
        w.bases.assign(length, 'A');
        w.nodes.resize(length + 1);
        std::iota(w.nodes.begin(), w.nodes.end(), branch);
        w.nodes[0] = 0;
        return w;
    };
    auto long_path = walk(250, 10, 1000), short_path = walk(80, 10, 2000);
    require(workflow.choose_walk({long_path, short_path}).bases.size() == 250,
            "short dead end cannot erase supported extension");
    auto weak_path = walk(250, 2, 3000);
    require(workflow.choose_walk({long_path, weak_path}).ambiguous,
            "absolute abundance alone cannot resolve otherwise equal paths");
    auto supported_fork = walk(750, 10, 10000), smooth_wrong_fork = walk(750, 200, 20000);
    supported_fork.branch_observations = smooth_wrong_fork.branch_observations = 10;
    supported_fork.branch_loss = 1.4;
    smooth_wrong_fork.branch_loss = 2.3;
    // Use real count changes: choose_walk recomputes deviation from the profile.
    supported_fork.counts[300] = 11.;
    workflow.update_profile(supported_fork, true);
    workflow.update_profile(smooth_wrong_fork, true);
    require(supported_fork.transition_energy > 0. && smooth_wrong_fork.transition_energy == 0. &&
                workflow.choose_walk({supported_fork, smooth_wrong_fork}).nodes == supported_fork.nodes,
            "supported branch decisions outweigh a small measured count change");
    auto tied_path = walk(250, 10, 4000);
    require(workflow.choose_walk({long_path, tied_path}).ambiguous &&
                workflow.choose_walk({long_path, tied_path}).bases.size() == 250,
            "a competing route is annotated without destroying the selected "
            "actual path");
    ExtensionWalk continuous, split;
    std::vector<double> trace{8, 8, 9, 8, 1, 8, 9, 8};
    continuous.counts = trace;
    workflow.update_profile(continuous, true);
    for (double c : trace) {
        split.counts.push_back(c);
        workflow.update_profile(split);
    }
    require(std::abs(continuous.deviation - split.deviation) < 1e-12 && continuous.deviation > 0.,
            "raw path evidence is invariant to arbitrary unitig segmentation");
    ExtensionWalk masked;
    masked.bases.assign(7, 'A');
    masked.nodes.resize(8);
    masked.all_counts.assign(7, 10.);
    masked.counts = {10., 11., 12.};
    masked.capacities = {10., 10., 10.};
    masked.sample_positions = {4, 5, 7};
    masked.branch_costs = {9., 9., 9., 6., 2., 9., 8.};
    workflow.update_profile(masked);
    require(masked.branch_loss == 2. && masked.branch_observations == 1 && masked.transitions == 1,
            "branch and continuity evidence use the same adjacent downstream mask");
    workflow.trim_walk(masked, 4, "test_trim");
    require(masked.branch_loss == 0. && masked.branch_observations == 0 && masked.transitions == 0,
            "trimming removes both fork and continuity evidence beyond the saved "
            "prefix");
    auto bubble_a = walk(250, 10, 1000), bubble_b = bubble_a;
    for (size_t i = 50; i < 74; ++i)
        bubble_b.nodes[i] += 10000;
    require(workflow.choose_walk({bubble_a, bubble_b}).bases.size() == 250,
            "rejoining single-variant bubble preserves downstream context");
    bubble_b.nodes[180] += 10000;
    require(workflow.same_route(bubble_a, bubble_b), "separate rejoining bubbles remain one route class");
    for (size_t i = 180; i < bubble_b.nodes.size(); ++i)
        bubble_b.nodes[i] += 20000;
    require(!workflow.same_route(bubble_a, bubble_b), "persistent divergence cannot be merged as a bubble");
    auto longer = walk(750, 9, 5000);
    std::vector<ExtensionWalk> comparable{long_path, longer};
    workflow.rank_walks(comparable, 2, false);
    require(comparable.front().bases.size() == 750, "prefer length within comparable coverage support");
    auto sibling = long_path;
    std::vector<ExtensionWalk> families{long_path, sibling, longer};
    workflow.rank_walks(families, 2, true);
    require(families[0].nodes.back() != families[1].nodes.back(), "beam keeps distinct current endpoints");
    auto short_dead_end = walk(30, 100, 6000), valid_extension = walk(750, 9, 7000);
    valid_extension.deviation = 1;
    std::vector<ExtensionWalk> terminal_buffer(300, short_dead_end);
    terminal_buffer.push_back(valid_extension);
    workflow.retain_completed(terminal_buffer);
    require(terminal_buffer.size() == 1 && terminal_buffer.front().bases.size() == 750,
            "ineligible dead ends cannot displace a completed extension");
    TargetCandidate exact, supported;
    exact.detected = true;
    exact.identity = exact.coverage = 1.;
    exact.target = "ACGT";
    supported.detected = supported.extended = true;
    supported.identity = .94;
    supported.coverage = 1.;
    supported.target = "ACGA";
    std::vector<TargetCandidate> placements{exact, supported};
    workflow.select_candidate(placements);
    require(placements.front().extended && placements.front().identity == .94,
            "an unextended exact placement cannot suppress a qualifying extended "
            "placement");
    exact.extended = true;
    placements = {supported, exact};
    workflow.select_candidate(placements);
    require(placements.front().identity == 1., "an already extended better alignment remains preferred");
    supported.identity = 1.;
    placements = {supported, exact};
    workflow.select_candidate(placements);
    require(!placements.front().detected && placements.front().reason == "multiple_target_placements",
            "equally aligned different extended targets remain unresolved");
    auto padded = walk(196, 10, 8000), unpadded = walk(199, 10, 9000);
    padded.anchor_padding = 4;
    workflow.options_.min_external = 150;
    std::vector<ExtensionWalk> boundary{padded, unpadded};
    workflow.retain_completed(boundary);
    require(boundary.size() == 1 && boundary.front().anchor_padding == 4,
            "total minimum flank includes padding while independent support has "
            "its own threshold");
    workflow.options_.min_external = 200;
    TargetCandidate strict, recovered, contained;
    strict.detected = strict.extended = true;
    strict.target = "ACGT";
    strict.sequence = "ACGTACCGGTTAACCGGTT";
    recovered = strict;
    recovered.allow_crossing = true;
    recovered.sequence = "CCC" + strict.sequence + "GGG";
    contained = strict;
    contained.sequence = strict.sequence.substr(2, 10);
    workflow.results_ = {{strict}, {recovered}, {contained}};
    workflow.consolidate();
    require(workflow.results_[0][0].extension_id != workflow.results_[1][0].extension_id,
            "recovery cannot replace an accepted strict output during "
            "containment consolidation");
    require(workflow.results_[0][0].extension_id == workflow.results_[2][0].extension_id,
            "strict-pass containment still consolidates as before");
    std::cout << "Workflow arithmetic and context checks passed. No dataset was run.\n";
}
} // namespace mcaat::proto
