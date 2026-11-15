#include "protospacer_detection/workflow.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace mcaat::proto {
namespace fs = std::filesystem;
static constexpr uint64_t missing = std::numeric_limits<uint64_t>::max();

volatile std::sig_atomic_t ProtospacerWorkflow::interrupted_ = 0;
void ProtospacerWorkflow::interrupt(int signal) { interrupted_ = signal; }

ProtospacerWorkflow::ProtospacerWorkflow(WorkflowOptions options) : options_(std::move(options)) {}

void ProtospacerWorkflow::print_help() {
    std::cout
        << "MCAAT protospacer workflow v2\n"
           "mcaat phage-curate --arrays DIR --graph PREFIX --output NEW_DIR [options]\n"
           "One workflow: spacer alignment -> source-context check -> extension -> filtering -> output.\n"
           "--threads 1 --beam-width 64 --max-flank-length 750 --min-extension-flank 200\n"
           "--min-identity 0.90 --min-query-coverage 0.85 --array-repeat-identity 0.85\n"
           "--max-query-states 2000000 --max-placement-states 400000 --max-placements 128\n"
           "--branch-ratio 1.5 --max-node-visits 2\n"
           "--correct-spacer-boundaries true --boundary-support 0.85\n"
           "The flank limit is a maximum. Ambiguous branches stop at their supported common prefix.\n"
           "Source-array traversal never establishes a target. Graph files are never modified.\n"
           "No supplied-target preservation mode and no reference-guided path selection.\n"
           "Output must be new or empty. All counts and search limits are recorded in run_summary.json.\n"
           "--self-test runs arithmetic/context checks without a graph or biological dataset.\n";
}

WorkflowOptions ProtospacerWorkflow::parse_options(int argc, char **argv) {
    WorkflowOptions o;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (++i >= argc)
            throw std::invalid_argument("missing value for " + flag);
        std::string value = argv[i];
        if (flag == "--arrays")
            o.arrays = value;
        else if (flag == "--graph")
            o.graph = value;
        else if (flag == "--output")
            o.output = value;
        else if (flag == "--correct-spacer-boundaries") {
            if (value != "true" && value != "false")
                throw std::invalid_argument(flag + " requires true or false");
            o.correct_boundaries = value == "true";
        } else if (flag == "--threads" || flag == "--beam-width" || flag == "--max-flank-length" ||
                   flag == "--min-extension-flank" || flag == "--max-query-states" ||
                   flag == "--max-placement-states" || flag == "--max-placements" ||
                   flag == "--max-node-visits") {
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
        } else if (flag == "--min-identity" || flag == "--min-query-coverage" ||
                   flag == "--array-repeat-identity" || flag == "--boundary-support" ||
                   flag == "--branch-ratio") {
            size_t used = 0;
            double n = std::stod(value, &used);
            if (used != value.size() || !std::isfinite(n))
                throw std::invalid_argument("invalid " + flag);
            if (flag == "--branch-ratio") {
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
    if (!fs::is_directory(o.arrays))
        throw std::invalid_argument("array directory does not exist");
    if (fs::is_directory(o.graph))
        o.graph /= "graph";
    if (!fs::is_regular_file(o.graph.string() + ".sdbg_info"))
        throw std::invalid_argument("graph metadata does not exist");
    if (o.threads > 1024 || o.beam_width > 256 || o.max_placements > 256 || o.max_visits > 4 ||
        o.max_flank > 10000 || o.min_flank > o.max_flank || o.placement_work > o.max_work)
        throw std::invalid_argument("inconsistent or excessive workflow limits");
    if (fs::exists(o.output) && (!fs::is_directory(o.output) || !fs::is_empty(o.output)))
        throw std::invalid_argument("output must be new or empty; existing results are preserved");
    o.arrays = fs::absolute(o.arrays);
    o.graph = fs::absolute(o.graph);
    o.output = fs::absolute(o.output);
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

double ProtospacerWorkflow::median(std::vector<double> v) {
    if (v.empty())
        return 0;
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
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

void ProtospacerWorkflow::group_queries() {
    std::map<std::string, std::vector<size_t>> grouped;
    for (size_t i = 0; i < sources_.size(); ++i)
        grouped[canonical(sources_[i].spacer)].push_back(i);
    queries_.assign(grouped.begin(), grouped.end());
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
    if (length < k)
        return found;
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
        if (source.spacer != sequence)
            repeat = reverse_complement(repeat);
        auto remainder = backwards ? sequence.substr(0, candidate.begin) : sequence.substr(candidate.end);
        if (backwards) {
            std::reverse(repeat.begin(), repeat.end());
            std::reverse(remainder.begin(), remainder.end());
        }
        int state = repeat_context(outward, remainder, repeat, options_.repeat_identity);
        if (state == 1)
            return 1;
        pending |= state == 0;
    }
    return pending ? 0 : 2;
}

void ProtospacerWorkflow::trim_walk(ExtensionWalk &walk, size_t length, const std::string &reason) const {
    length = std::min(length, walk.bases.size());
    walk.bases.resize(length);
    if (!walk.nodes.empty())
        walk.nodes.resize(length + 1);
    walk.stop = reason;
    while (!walk.sample_positions.empty() && walk.sample_positions.back() > length) {
        walk.sample_positions.pop_back();
        walk.counts.pop_back();
        walk.capacities.pop_back();
    }
    if (length < walk.clear_at) {
        walk.external = false;
        walk.clear_at = 0;
    }
    walk.flow = walk.capacities.empty() ? 0 : walk.capacities.back();
}

Unitig ProtospacerWorkflow::build_unitig(uint64_t start, const GraphArc &first, bool backwards, size_t &work,
                                         size_t work_limit) const {
    Unitig unitig;
    std::set<uint64_t> seen{start};
    GraphArc arc = first;
    while (unitig.arcs.size() < options_.max_flank && !interrupted_) {
        if (work++ >= work_limit) {
            unitig.limited = true;
            break;
        }
        unitig.arcs.push_back(arc);
        if (!seen.insert(arc.node).second)
            break;
        if (neighbors(arc.node, !backwards).size() != 1)
            break;
        auto next = neighbors(arc.node, backwards);
        if (next.size() != 1)
            break;
        arc = next.front();
    }
    return unitig;
}

void ProtospacerWorkflow::update_flow(ExtensionWalk &walk, size_t before) const {
    if (walk.counts.size() == before)
        return;
    // The median is a robust unitig capacity, not an independent-read count.
    // Max-min path flow tolerates a shared repeat's increased coverage.
    std::vector<double> counts(walk.counts.begin() + before, walk.counts.end());
    double capacity = median(std::move(counts));
    walk.flow = before ? std::min(walk.flow, capacity) : capacity;
    for (size_t i = before; i < walk.capacities.size(); ++i)
        walk.capacities[i] = walk.flow;
    ++walk.unitigs;
}

bool ProtospacerWorkflow::advance_walk(ExtensionWalk &walk, const GraphArc &arc, size_t query,
                                       const TargetCandidate &candidate, bool backwards) const {
    if (std::find(candidate.anchor.begin(), candidate.anchor.end(), arc.node) != candidate.anchor.end()) {
        walk.stop = "anchor_revisit";
        return false;
    }
    if (size_t(std::count(walk.nodes.begin(), walk.nodes.end(), arc.node)) >= options_.max_visits) {
        walk.stop = "repeat_visit_limit";
        return false;
    }
    walk.nodes.push_back(arc.node);
    walk.bases += added_base(arc, backwards);
    if (!walk.external) {
        int state = source_context(query, candidate, walk.bases, backwards);
        if (state == 1) {
            trim_walk(walk, 0, "source_array_context");
            return false;
        }
        if (state == 2) {
            walk.external = true;
            walk.clear_at = walk.bases.size();
        }
    }
    if (walk.external && walk.bases.size() > std::max(size_t(graph_->k()), walk.clear_at)) {
        if (arc.count >= 65535) {
            walk.stop = "saturated_multiplicity";
            return false;
        }
        walk.counts.push_back(arc.count);
        walk.sample_positions.push_back(walk.bases.size());
        walk.capacities.push_back(0);
    }
    if (walk.bases.size() >= options_.max_flank) {
        walk.stop = "maximum_flank";
        return false;
    }
    return true;
}

ExtensionWalk ProtospacerWorkflow::choose_walk(std::vector<ExtensionWalk> walks) const {
    walks.erase(std::remove_if(walks.begin(), walks.end(),
                               [&](const auto &w) {
                                   return !w.external ||
                                          w.bases.size() < std::max(options_.context_bases - 1, w.clear_at) +
                                                               options_.context_bases;
                               }),
                walks.end());
    if (walks.empty()) {
        ExtensionWalk empty;
        empty.stop = "no_supported_external_context";
        return empty;
    }
    bool has_extension = std::any_of(walks.begin(), walks.end(),
                                     [&](const auto &w) { return w.bases.size() >= options_.min_flank; });
    if (has_extension)
        walks.erase(std::remove_if(walks.begin(), walks.end(),
                                   [&](const auto &w) { return w.bases.size() < options_.min_flank; }),
                    walks.end());
    std::sort(walks.begin(), walks.end(), [](const auto &a, const auto &b) {
        if (a.flow != b.flow)
            return a.flow > b.flow;
        if (a.bases.size() != b.bases.size())
            return a.bases.size() > b.bases.size();
        return a.nodes < b.nodes;
    });
    ExtensionWalk best = walks.front();
    size_t common = best.nodes.size();
    for (size_t i = 1; i < walks.size(); ++i) {
        const auto &other = walks[i];
        if (other.flow * options_.branch_ratio < best.flow)
            continue;
        size_t n = 0;
        while (n < std::min(best.nodes.size(), other.nodes.size()) && best.nodes[n] == other.nodes[n])
            ++n;
        // Prefixes of the same path are not competing biological routes.
        if (n == std::min(best.nodes.size(), other.nodes.size()))
            continue;
        common = std::min(common, n);
    }
    if (common < best.nodes.size())
        trim_walk(best, common ? common - 1 : 0, "ambiguous_flow");
    return best;
}

ExtensionWalk ProtospacerWorkflow::extend_flank(size_t query, const TargetCandidate &candidate,
                                                bool backwards, size_t &work, size_t work_limit) const {
    ExtensionWalk seed;
    seed.nodes = {backwards ? candidate.anchor.front() : candidate.anchor.back()};
    std::vector<ExtensionWalk> frontier{seed}, finished;
    std::map<uint64_t, Unitig> unitigs;
    bool limited = false;
    while (!frontier.empty() && !interrupted_) {
        std::vector<ExtensionWalk> next;
        for (auto &current : frontier) {
            if (work >= work_limit) {
                current.limited = true;
                current.stop = "query_work_limit";
                finished.push_back(std::move(current));
                limited = true;
                continue;
            }
            auto arcs = neighbors(current.nodes.back(), backwards);
            if (arcs.empty()) {
                current.stop = "dead_end";
                finished.push_back(std::move(current));
                continue;
            }
            for (const auto &arc : arcs) {
                if (!unitigs.count(arc.edge))
                    unitigs.emplace(arc.edge,
                                    build_unitig(current.nodes.back(), arc, backwards, work, work_limit));
                const auto &segment = unitigs.at(arc.edge);
                auto child = current;
                size_t before = child.counts.size();
                bool continuing = true;
                for (const auto &edge : segment.arcs) {
                    if (work++ >= work_limit) {
                        child.stop = "query_work_limit";
                        child.limited = true;
                        limited = true;
                        continuing = false;
                        break;
                    }
                    ++child.work;
                    if (!advance_walk(child, edge, query, candidate, backwards)) {
                        continuing = false;
                        break;
                    }
                }
                if (child.counts.size() >= before)
                    update_flow(child, before);
                if (segment.limited) {
                    child.limited = true;
                    child.stop = "query_work_limit";
                    continuing = false;
                    limited = true;
                }
                if (continuing && !segment.arcs.empty())
                    next.push_back(std::move(child));
                else
                    finished.push_back(std::move(child));
            }
        }
        // Unitigs have different lengths. Rank progress only after capacity;
        // every retained state still traverses actual graph arcs.
        std::sort(next.begin(), next.end(), [](const auto &a, const auto &b) {
            if (a.external != b.external)
                return a.external > b.external;
            if (a.flow != b.flow)
                return a.flow > b.flow;
            if (a.bases.size() != b.bases.size())
                return a.bases.size() > b.bases.size();
            return a.nodes < b.nodes;
        });
        if (next.size() > options_.beam_width) {
            // Save pruned prefixes as unresolved alternatives. They can prevent
            // an unsupported unique call but are never treated as completed paths.
            for (size_t i = options_.beam_width; i < next.size(); ++i) {
                next[i].limited = true;
                next[i].stop = "beam_limit";
                finished.push_back(std::move(next[i]));
            }
            next.resize(options_.beam_width);
            limited = true;
        }
        frontier = std::move(next);
        if (finished.size() > options_.beam_width * 32) {
            for (auto &w : frontier) {
                w.limited = true;
                w.stop = "candidate_limit";
                finished.push_back(std::move(w));
            }
            frontier.clear();
            limited = true;
        }
    }
    for (auto &w : frontier)
        finished.push_back(std::move(w));
    auto chosen = choose_walk(std::move(finished));
    chosen.limited |= limited;
    return chosen;
}

void ProtospacerWorkflow::assemble_candidate(TargetCandidate &c) const {
    auto left = c.left.bases;
    std::reverse(left.begin(), left.end());
    c.sequence = left + c.target + c.right.bases;
    auto valid = [&](const ExtensionWalk &w) {
        return w.external &&
               w.bases.size() >= std::max(size_t(graph_->k()), w.clear_at) + options_.context_bases;
    };
    c.detected = valid(c.left) || valid(c.right);
    if (!c.detected) {
        c.status = "rejected";
        c.reason = "no_independent_target_context";
        return;
    }
    c.context_left = c.left.bases;
    c.context_right = c.right.bases;
    c.context_left_stop = c.left.stop;
    c.context_right_stop = c.right.stop;
    c.clear_left = c.left.clear_at;
    c.clear_right = c.right.clear_at;
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
    if (!informative(c.target)) {
        c.detected = false;
        c.status = "unresolved";
        c.reason = "uninformative_target_sequence";
        return;
    }
    for (bool back : {true, false}) {
        auto &w = back ? c.left : c.right;
        if (w.bases.size() < options_.min_flank || w.flow <= 0 || !informative(w.bases))
            trim_walk(w, 0, "insufficient_supported_extension");
    }
    // Each flank is supported separately. Do not silently join incompatible
    // abundance contexts merely because their seed sequence is shared.
    if (!c.left.bases.empty() && !c.right.bases.empty()) {
        double l = c.left.flow, r = c.right.flow;
        bool conflict = std::max(l, r) > 4 * std::min(l, r);
        std::set<uint64_t> nodes(c.left.nodes.begin() + 1, c.left.nodes.end());
        for (auto n : c.right.nodes)
            conflict |= nodes.count(n) > 0;
        if (conflict) {
            if (l >= r)
                trim_walk(c.right, 0, "incompatible_flanks");
            else
                trim_walk(c.left, 0, "incompatible_flanks");
            c.reason = "incompatible_flanks";
        }
    }
    left = c.left.bases;
    std::reverse(left.begin(), left.end());
    c.sequence = left + c.target + c.right.bases;
    c.extended = !c.left.bases.empty() || !c.right.bases.empty();
    if (c.extended) {
        c.status = "extended";
        c.reason = "supported_extension";
    } else if (c.reason == "target_context_passed")
        c.reason = "target_detected_extension_unresolved";
}

void ProtospacerWorkflow::process_query(size_t query) {
    if (!informative(queries_[query].first)) {
        query_status_[query] = "uninformative_query";
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
    for (size_t i = 0; i < anchors.size() && !interrupted_; ++i) {
        auto &candidate = anchors[i];
        candidate.limited = limited;
        size_t budget = (options_.max_work - std::min(work, options_.max_work)) / (anchors.size() - i);
        size_t left_limit = work + budget / 2, right_limit = work + budget;
        start = std::chrono::steady_clock::now();
        candidate.left = extend_flank(query, candidate, true, work, left_limit);
        candidate.right = extend_flank(query, candidate, false, work, right_limit);
        query_timings_[query][1] +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        candidate.limited |= candidate.left.limited || candidate.right.limited || work >= options_.max_work;
        query_limited_[query] |= candidate.limited;
        start = std::chrono::steady_clock::now();
        assemble_candidate(candidate);
        query_timings_[query][2] +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        checked.push_back(std::move(candidate));
        if (work >= options_.max_work)
            break;
    }
    std::sort(checked.begin(), checked.end(), [](const auto &a, const auto &b) {
        if (a.detected != b.detected)
            return a.detected > b.detected;
        if (a.coverage != b.coverage)
            return a.coverage > b.coverage;
        if (a.identity != b.identity)
            return a.identity > b.identity;
        return a.anchor < b.anchor;
    });
    if (!checked.empty() && checked.front().detected) {
        auto &first = checked.front();
        bool ambiguous = false;
        for (size_t i = 1; i < checked.size(); ++i) {
            auto &other = checked[i];
            if (other.detected && other.coverage == first.coverage && other.identity == first.identity &&
                other.target != first.target)
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
    query_status_[query] = checked.empty() ? "interrupted" : checked.front().reason;
    query_work_[query] = work;
    results_[query] = std::move(checked);
}

void ProtospacerWorkflow::consolidate() {
    std::map<std::string, std::string> targets, regions;
    for (auto &group : results_)
        for (auto &c : group)
            if (c.detected) {
                auto key = canonical(c.target);
                targets.emplace(key, "");
                if (c.extended)
                    regions.emplace(canonical(c.sequence), "");
            }
    size_t number = 0;
    for (auto &item : targets)
        item.second = numbered("P", ++number);
    number = 0;
    for (auto &item : regions)
        item.second = numbered("E", ++number);
    // Assign each contained fragment to a fixed longest representative.
    // Never merge two partially overlapping paths into an invented sequence.
    std::vector<std::string> order;
    for (const auto &item : regions)
        order.push_back(item.first);
    std::sort(order.begin(), order.end(), [](const auto &a, const auto &b) {
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
    std::map<std::string, std::string> targets, regions;
    for (const auto &group : results_)
        for (const auto &c : group)
            if (c.detected) {
                targets[c.target_id] = canonical(c.target);
                if (c.extended && c.sequence.size() > regions[c.extension_id].size())
                    regions[c.extension_id] = canonical(c.sequence);
            }
    {
        auto target_fasta = open("Protospacers.fasta"), extension_fasta = open("ContiguousPaths.fasta");
        for (const auto &[id, s] : targets)
            target_fasta << '>' << id << '\n' << s << '\n';
        for (const auto &[id, s] : regions)
            extension_fasta << '>' << id << '\n' << s << '\n';
    }
    {
        auto input = open("sources.tsv"), queries = open("queries.tsv"), candidates = open("candidates.tsv"),
             rejected = open("unresolved_and_rejected.tsv");
        input << "source_id\tarray_id\toriginal_spacer\toriginal_repeat\tspacer\trepeat\ttrim_left\ttrim_"
                 "right\n";
        for (const auto &s : sources_)
            input << s.id << '\t' << s.array << '\t' << s.original_spacer << '\t' << s.original_repeat << '\t'
                  << s.spacer << '\t' << s.repeat << '\t' << s.trim_left << '\t' << s.trim_right << '\n';
        queries << "query_id\tspacer\tsource_count\tstatus\twork\tsearch_limited\n";
        auto metrics = open("query_metrics.tsv");
        metrics
            << "query_id\tmapping_work\textension_work\tmapping_seconds\textension_seconds\tfilter_seconds\n";
        candidates << "query_id\ttarget_id\tstatus\treason\tquery_begin\tquery_end\tidentity\tquery_"
                      "coverage\tsequence\tleft_bp\tright_bp\tleft_stop\tright_stop\tsearch_limited\n";
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
                           << c.target << '\t' << c.left.bases.size() << '\t' << c.right.bases.size() << '\t'
                           << c.left.stop << '\t' << c.right.stop << '\t' << c.limited << '\n';
                if (!c.extended)
                    rejected << id << '\t' << c.status << '\t' << c.reason << '\n';
            }
        }
    }
    {
        auto links = open("associations.tsv"), evidence = open("coverage_evidence.tsv");
        links << "source_id\tquery_id\ttarget_id\textension_id\tquery_begin\tquery_end\tidentity\tquery_"
                 "coverage\ttarget_strand\tseed_begin\tseed_end\tsupported_fragment_begin\tsupported_"
                 "fragment_end\tleft_bp\tright_bp\tleft_stop\tright_stop\tsearch_limited\n";
        evidence << "target_id\textension_id\tside\tadded_base_position\tcount\tpath_capacity\n";
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
                        seed_begin = offset + (reverse ? c.right.bases.size() : c.left.bases.size());
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
                        links << source.id << '\t' << numbered("Q", c.query + 1) << '\t' << c.target_id
                              << '\t' << c.extension_id << '\t' << begin + source.trim_left << '\t'
                              << end + source.trim_left << '\t' << c.identity << '\t' << c.coverage << '\t'
                              << ((reverse ^ source_flip) ? '-' : '+') << '\t';
                        if (c.extended)
                            links << seed_begin << '\t' << seed_end << '\t' << offset << '\t'
                                  << offset + c.sequence.size();
                        else
                            links << "\t\t\t";
                        links << '\t' << c.left.bases.size() << '\t' << c.right.bases.size() << '\t'
                              << c.left.stop << '\t' << c.right.stop << '\t' << c.limited << '\n';
                    }
                    for (bool back : {true, false}) {
                        const auto &w = back ? c.left : c.right;
                        for (size_t j = 0; j < w.counts.size(); ++j)
                            evidence << c.target_id << '\t' << c.extension_id << '\t'
                                     << (back ? "left" : "right") << '\t' << w.sample_positions[j] << '\t'
                                     << w.counts[j] << '\t' << w.capacities[j] << '\n';
                    }
                }
    }
    {
        auto contexts = open("target_context.tsv");
        contexts << "query_id\ttarget_id\tleft_outward\tright_outward\tleft_clear_at\tright_clear_at\tleft_"
                    "stop_before_filter\tright_stop_before_filter\n";
        auto paths = open("graph_paths.tsv");
        paths << "query_id\ttarget_id\textension_id\tanchor_nodes\tleft_nodes_outward\tright_nodes_outward\n";
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
                    paths << '\n';
                }
    }
    for (const auto &f : fs::directory_iterator(staging))
        fs::rename(f.path(), options_.output / f.path().filename());
    fs::remove(staging);
}

void ProtospacerWorkflow::write_summary(const std::string &status, const std::string &error) const {
    std::set<std::string> arrays, targets, extensions;
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
                targets.insert(c.target_id);
            if (c.extended)
                extensions.insert(c.extension_id);
        }
        detected += d;
        extended += e;
        limited += l;
    }
    limited = std::count(query_limited_.begin(), query_limited_.end(), uint8_t(1));
    fs::path temporary = options_.output / "run_summary.json.tmp";
    std::ofstream out(temporary);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << "{\n  \"implementation\": \"protospacer-workflow-v2\",\n  \"status\": " << json(status)
        << ",\n  \"error\": " << json(error) << ",\n  \"arrays\": " << json(options_.arrays.string())
        << ",\n  \"graph\": " << json(options_.graph.string())
        << ",\n  \"graph_k\": " << (graph_ ? graph_->k() : 0) << ",\n  \"array_count\": " << arrays.size()
        << ",\n  \"spacer_observations\": " << sources_.size()
        << ",\n  \"distinct_spacers\": " << queries_.size() << ",\n  \"detected_queries\": " << detected
        << ",\n  \"protospacer_regions\": " << targets.size() << ",\n  \"extended_queries\": " << extended
        << ",\n  \"unique_extensions\": " << extensions.size()
        << ",\n  \"queries_without_extension\": " << queries_.size() - extended
        << ",\n  \"limited_queries\": " << limited << ",\n  \"threads\": " << options_.threads
        << ",\n  \"beam_width\": " << options_.beam_width
        << ",\n  \"max_flank_length\": " << options_.max_flank
        << ",\n  \"min_extension_flank\": " << options_.min_flank
        << ",\n  \"min_identity\": " << options_.identity
        << ",\n  \"min_query_coverage\": " << options_.coverage
        << ",\n  \"array_repeat_identity\": " << options_.repeat_identity
        << ",\n  \"boundary_support\": " << options_.boundary_support
        << ",\n  \"correct_spacer_boundaries\": " << (options_.correct_boundaries ? "true" : "false")
        << ",\n  \"max_query_states\": " << options_.max_work
        << ",\n  \"max_placement_states\": " << options_.placement_work
        << ",\n  \"max_placements\": " << options_.max_placements
        << ",\n  \"branch_ratio\": " << options_.branch_ratio
        << ",\n  \"max_node_visits\": " << options_.max_visits
        << ",\n  \"extension_method\": \"unitig_widest_path\""
        << ",\n  \"probabilities_of_correctness\": false,\n  \"reference_used_for_selection\": false,\n  "
           "\"source_array_crossing\": false,\n  \"timings_seconds\": {";
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
    group_queries();
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
    std::cout << "Workflow complete. Final counts: " << (options_.output / "run_summary.json") << '\n';
    return interrupted_ ? 128 + interrupted_ : 0;
}

void ProtospacerWorkflow::self_test() {
    auto require = [](bool condition, const std::string &reason) {
        if (!condition)
            throw std::runtime_error("self-test failed: " + reason);
    };
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
    ExtensionWalk coverage;
    for (double count : {20., 25., 20.}) {
        size_t before = coverage.counts.size();
        coverage.counts.push_back(count);
        coverage.capacities.push_back(0);
        workflow.update_flow(coverage, before);
    }
    require(coverage.flow == 20., "shared repeat does not penalize widest-path capacity");
    auto walk = [](size_t length, double flow, uint64_t branch) {
        ExtensionWalk w;
        w.external = true;
        w.clear_at = 1;
        w.flow = flow;
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
    require(workflow.choose_walk({long_path, weak_path}).nodes == long_path.nodes,
            "widest path wins against weak alternative");
    auto tied_path = walk(250, 10, 4000);
    require(workflow.choose_walk({long_path, tied_path}).bases.empty(),
            "equal flow preserves branch ambiguity");
    std::cout << "Workflow arithmetic and context checks passed. No dataset was run.\n";
}
} // namespace mcaat::proto
