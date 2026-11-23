// FindProtospacers and ProtospacerExtend as defined in the thesis, with the
// multiplicity criterion of the beam search exchangeable.
//
// For every array the tangle T is the node set of its repeat and spacers. For
// every node v of a spacer, the neighbours outside T are collected: those that
// point into v give B_in, those v points to give B_out. A spacer is a candidate
// only when both sets are non-empty. Each node of B_out is extended forward and
// each node of B_in backward by a beam search of width 50, up to 1250 added
// bases, admitting a node w when mult(w) lies in [0.1*mult(b), 5*mult(b)] for
// the boundary node b the walk started from and w is outside T. The best path of
// each side is kept, and the output is left + spacer + right.
//
// Only graph topology and multiplicity are read. The criteria differ solely in
// how a partial path is scored.
//
// Two things differ from the thesis, both measured on SRR14017565 and both
// defaults here. The window is keyed to the running median of the walk instead
// of to the boundary node. And the nodes of a spacer are not required to form
// one connected path, which the thesis never asks for: requiring it discarded
// 621 of 1,962 spacer observations, and dropping it took the two-sided
// candidates from 245 to 575 and the strict matches from 22 to 93. Both sides
// are still required.
//
// A walk is kept when the harmonic mean of its node multiplicities, pooled over
// both flanks, is at least --min-harmonic-mean. At 3 on that sample the run
// gives 268 contigs with 93 strict matches at 56 loci, against 286 contigs with
// 54 matches at 5 loci from mcaat phage-curate. On the same run at 5, 52 of its
// 92 matches have the alignment covering the spacer and at least 50 bases on
// each side of it; the production run has none at all, every one of its 286
// contigs being extended on one side only. Those figures are for harmonic_mean, which beat the
// thesis's plain mean by 1.37x. One sample only, so far.
//
//   boundary_walk GRAPH_PREFIX SOURCES_TSV OUT_DIR [--flank N] [--beam N]
//
// Writes OUT_DIR/<criterion>.fasta and OUT_DIR/summary.tsv.
#include "sdbg/sdbg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {
constexpr uint64_t missing = std::numeric_limits<uint64_t>::max();
const std::string ALPHABET = "ACGT";

std::string reverse_complement(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char &c : s) {
        const auto i = ALPHABET.find(c);
        c = i == std::string::npos ? c : ALPHABET[3 - i];
    }
    return s;
}

bool dna(const std::string &s) { return !s.empty() && s.find_first_not_of("ACGT") == std::string::npos; }

struct Arc {
    uint64_t node = 0, edge = 0;
    double mult = 0;
};

// The six scoring criteria. Every one is maximised.
enum Criterion { Mean, LengthWeighted, Geometric, Harmonic, BoundaryStable, LocallySmooth, Bottleneck, Count };
const std::array<const char *, Count> NAMES{"mean_multiplicity", "length_weighted_mean", "geometric_mean",
                                            "harmonic_mean",     "boundary_stable",      "locally_smooth",
                                            "bottleneck"};

struct Walk {
    std::vector<uint64_t> nodes;
    std::string bases;
    double accumulated = 0., score = 0., first = 0., last = 0., lowest = 0.;
    std::vector<double> multiplicities;
};

class Graph {
  public:
    explicit Graph(const char *prefix) { sdbg_.LoadFromFile(prefix); }
    size_t k() const { return sdbg_.k(); }
    uint64_t size() const { return sdbg_.size(); }

    std::string label(uint64_t node) const {
        if (node >= sdbg_.size()) return {};
        std::vector<uint8_t> encoded(sdbg_.k());
        sdbg_.GetLabel(node, encoded.data());
        std::string s;
        for (auto c : encoded) {
            if (c < 1 || c > 4) return {};
            s += ALPHABET[c - 1];
        }
        return s;
    }

    uint64_t find(const std::string &word) const {
        if (word.size() != k() || !dna(word)) return missing;
        std::vector<uint8_t> encoded;
        for (char c : word) encoded.push_back(ALPHABET.find(c) + 1);
        const uint64_t n = sdbg_.IndexBinarySearch(encoded.data());
        return n != missing && label(n) == word ? n : missing;
    }

    // The library's own edge lists. A hand-written version of the forward scan
    // was wrong here: it collected the group around the node itself rather than
    // the group below its successor, which both invented steps the reads never
    // supported and read the multiplicity off a neighbouring edge.
    //
    // An edge id is the node it leads to. Going backwards the node an edge
    // comes from is the end of that edge's group, which is what GetLastIndex
    // gives.
    std::vector<Arc> neighbours(uint64_t node, bool incoming) const {
        std::vector<Arc> arcs;
        if (node >= sdbg_.size() || !sdbg_.IsValidEdge(node)) return arcs;
        std::array<uint64_t, 8> buffer{};
        const int found = incoming ? sdbg_.IncomingEdges(node, buffer.data())
                                   : sdbg_.OutgoingEdges(node, buffer.data());
        for (int i = 0; i < found; ++i) {
            const uint64_t edge = buffer[i];
            const uint64_t to = incoming ? uint64_t(sdbg_.GetLastIndex(edge)) : edge;
            if (to >= sdbg_.size()) continue;
            arcs.push_back({to, edge, double(sdbg_.EdgeMultiplicity(edge))});
        }
        return arcs;
    }

    // Every k-mer of the sequence that is in the graph, in whichever
    // orientation has more of them, and returns that orientation. The boundary
    // sets are collected over every node of a spacer; the thesis does not ask
    // those nodes to form one connected path, and one missing edge in the
    // middle should not discard the spacer. 183 of the 621 spacers that failed
    // to trace on SRR14017565 have both of their end k-mers in the graph.
    std::string locate(const std::string &seq, std::vector<uint64_t> &nodes) const {
        std::string chosen = seq;
        std::vector<uint64_t> best;
        for (const auto &oriented : {seq, reverse_complement(seq)}) {
            std::vector<uint64_t> found;
            for (size_t i = 0; i + k() <= oriented.size(); ++i) {
                const auto n = find(oriented.substr(i, k()));
                if (n != missing) found.push_back(n);
            }
            if (found.size() > best.size()) {
                best = std::move(found);
                chosen = oriented;
            }
        }
        nodes = std::move(best);
        return chosen;
    }

    // Nodes and entering multiplicities along a sequence, or empty when the
    // sequence is not a connected path.
    bool trace(const std::string &seq, std::vector<uint64_t> &nodes, std::vector<double> &mult) const {
        nodes.clear();
        mult.clear();
        if (seq.size() < k() || !dna(seq)) return false;
        auto node = find(seq.substr(0, k()));
        if (node == missing) return false;
        nodes.push_back(node);
        for (size_t at = k(); at < seq.size(); ++at) {
            bool stepped = false;
            for (const auto &arc : neighbours(node, false)) {
                if (base_of(arc, false) != seq[at]) continue;
                nodes.push_back(arc.node);
                mult.push_back(arc.mult);
                node = arc.node;
                stepped = true;
                break;
            }
            if (!stepped) return false;
        }
        return true;
    }

    // The base a step adds is the one the node it moves to ends with, going
    // forward, and the one it begins with, going backward. GetW on the edge is
    // not that base: it agreed with the node's own label in only 966 of 3,315
    // forward steps checked against the graph.
    char base_of(const Arc &arc, bool incoming) const {
        const auto s = label(arc.node);
        if (s.empty()) return 'N';
        return incoming ? s.front() : s.back();
    }

  private:
    SDBG sdbg_;
};

void score(Walk &w, int criterion) {
    const double n = double(w.nodes.size());
    switch (criterion) {
    case Mean: w.score = w.accumulated / n; break;
    case LengthWeighted: w.score = (w.accumulated / n) * std::sqrt(n + 1.); break;
    case Geometric: w.score = w.accumulated / n; break;          // accumulated holds sum of logs
    case Harmonic: w.score = -w.accumulated / n; break;          // accumulated holds sum of reciprocals
    case BoundaryStable: w.score = -w.accumulated / n; break;    // sum |log(c/c_first)|
    case LocallySmooth: w.score = n > 1 ? -w.accumulated / (n - 1.) : 0.; break;
    case Bottleneck: w.score = w.lowest; break;
    default: w.score = w.accumulated / n;
    }
}

void extend_one(Walk &w, const Arc &arc, int criterion, bool incoming, const Graph &g) {
    const double c = std::max(arc.mult, 1e-9);
    switch (criterion) {
    case Geometric: w.accumulated += std::log(c); break;
    case Harmonic: w.accumulated += 1. / c; break;
    case BoundaryStable: w.accumulated += std::abs(std::log(c / std::max(w.first, 1e-9))); break;
    case LocallySmooth:
        if (!w.nodes.empty()) w.accumulated += std::abs(std::log(c / std::max(w.last, 1e-9)));
        break;
    default: w.accumulated += c;
    }
    w.nodes.push_back(arc.node);
    w.multiplicities.push_back(c);
    w.bases += g.base_of(arc, incoming);
    w.last = c;
    w.lowest = w.nodes.empty() ? c : std::min(w.lowest, c);
    score(w, criterion);
}

struct Source {
    std::string id, array, spacer, repeat;
};

std::vector<Source> read_sources(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "cannot read " << path << "\n";
        std::exit(2);
    }
    std::string line;
    std::getline(in, line);
    std::vector<std::string> header;
    {
        std::istringstream h(line);
        std::string f;
        while (std::getline(h, f, '\t')) header.push_back(f);
    }
    auto column = [&](const std::string &name) {
        for (size_t i = 0; i < header.size(); ++i)
            if (header[i] == name) return i;
        std::cerr << "sources.tsv has no column " << name << "\n";
        std::exit(2);
    };
    const auto ci = column("source_id"), ca = column("array_id"), cs = column("spacer"), cr = column("repeat");
    std::vector<Source> out;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::istringstream row(line);
        std::string cell;
        while (std::getline(row, cell, '\t')) f.push_back(cell);
        if (f.size() <= std::max(std::max(ci, ca), std::max(cs, cr))) continue;
        out.push_back({f[ci], f[ca], f[cs], f[cr]});
    }
    return out;
}
} // namespace

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "usage: boundary_walk GRAPH_PREFIX SOURCES_TSV OUT_DIR [--flank N] "
                     "[--beam N] [--window-low F] [--window-high F] [--window-running true|false] "
                     "[--min-harmonic-mean F] [--require-connected-spacer true|false] "
                     "[--threads N]\n"
                     "defaults: --flank 1250 --beam 50 --window-low 0.1 --window-high 5 "
                     "--window-running true --min-harmonic-mean 3 "
                     "--require-connected-spacer false, every core for --threads\n";
        return 2;
    }
    size_t flank = 1250, beam = 50;
    double window_low = .1, window_high = 5.;
    // Keyed to the running median of the walk rather than to the boundary
    // node: a boundary node sitting on a low-coverage k-mer otherwise sets
    // a window the rest of the flank cannot stay inside.
    bool window_running = true;
    // The score picks the path out of the spacer; this decides whether the path
    // is kept. A walk whose node multiplicities have harmonic mean 1 was built
    // entirely out of k-mers seen once, which a sequencing error also produces.
    // On SRR14017565 every strict match survives up to 3 and all but one up to
    // 5, so 3 sits inside the plateau rather than on its edge.
    double min_harmonic = 3.;
    // The thesis's own wording asks for every node of the spacer, not for a
    // connected path through it. Left on, a single missing edge in the middle
    // discards the spacer; turned off, the boundary sets are collected over
    // whichever of its k-mers are in the graph. Off by default: requiring the
    // path discarded 621 of 1,962 spacer observations on SRR14017565, and 178
    // of the recovered candidates have every k-mer present with one edge
    // missing. Both boundary sets are still required to be non-empty.
    bool require_connected = false;
    // 0 leaves OpenMP its own default, which is every core.
    int threads = 0;
    for (int i = 4; i + 1 < argc; i += 2) {
        const std::string flag = argv[i];
        if (flag == "--flank") flank = std::stoul(argv[i + 1]);
        else if (flag == "--beam") beam = std::stoul(argv[i + 1]);
        else if (flag == "--window-low") window_low = std::stod(argv[i + 1]);
        else if (flag == "--window-high") window_high = std::stod(argv[i + 1]);
        else if (flag == "--window-running") window_running = std::string(argv[i + 1]) == "true";
        else if (flag == "--min-harmonic-mean") min_harmonic = std::stod(argv[i + 1]);
        else if (flag == "--require-connected-spacer")
            require_connected = std::string(argv[i + 1]) == "true";
        else if (flag == "--threads") threads = std::stoi(argv[i + 1]);
        else { std::cerr << "unknown option " << flag << "\n"; return 2; }
    }
    const auto sources = read_sources(argv[2]);
    std::filesystem::create_directories(argv[3]);
    Graph graph(argv[1]);
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    const int in_use = omp_get_max_threads();
#else
    (void)threads;
    const int in_use = 1;
#endif
    std::cerr << "graph     k=" << graph.k() << "  edges=" << graph.size() << "  threads=" << in_use
              << "\nsettings  flank " << flank << ", beam " << beam << ", window [" << window_low
              << ", " << window_high << "] on " << (window_running ? "the running median" : "the boundary node")
              << ", rule " << min_harmonic << "\n";

    // Tangle of each array: the nodes of its repeat and its spacers, both strands.
    std::map<std::string, std::unordered_set<uint64_t>> tangle;
    std::map<std::string, std::vector<const Source *>> by_array;
    for (const auto &s : sources) by_array[s.array].push_back(&s);
    std::vector<uint64_t> nodes;
    std::vector<double> mult;
    for (const auto &[array, members] : by_array) {
        auto &T = tangle[array];
        // repeat + spacer + repeat, so that the k-mers spanning the junctions
        // belong to the tangle as well. Built from pieces they would not, and a
        // walk leaving the spacer towards the array would run k-1 bases into the
        // junction before being stopped.
        for (const auto *s : members) {
            const auto context = s->repeat + s->spacer + s->repeat;
            for (const auto &oriented : {context, reverse_complement(context)})
                for (size_t i = 0; i + graph.k() <= oriented.size(); ++i) {
                    const auto n = graph.find(oriented.substr(i, graph.k()));
                    if (n != missing) T.insert(n);
                }
        }
    }
    std::cerr << "tangles built for " << tangle.size() << " arrays\n";

    std::array<std::ofstream, Count> out;
    for (int c = 0; c < Count; ++c) {
        out[c].open(std::string(argv[3]) + "/" + NAMES[c] + ".fasta");
        out[c].exceptions(std::ios::badbit | std::ios::failbit);
    }
    std::ofstream summary(std::string(argv[3]) + "/summary.tsv");
    summary << "criterion\tsource_id\tarray_id\tspacer_nodes\tb_in\tb_out\tleft_bp\tright_bp\tscore"
               "\tharmonic_mean\n";

    // One line per spacer, whether or not it became a candidate, so the places
    // where spacers are lost can be counted rather than guessed at.
    std::ofstream stages(std::string(argv[3]) + "/spacers.tsv");
    stages << "source_id\tarray_id\tspacer_bp\tspacer_nodes\tb_in\tb_out\tends_in_graph\tconnected"
              "\toutcome\n";
    // Whether the spacer's first and last k-mer are in the graph, on either
    // strand. Those two are what the boundary sets are read from, so a spacer
    // can have both of them and still fail to trace as one connected path.
    auto ends_present = [&](const std::string &seq) {
        if (seq.size() < graph.k()) return 0;
        int best = 0;
        for (const auto &oriented : {seq, reverse_complement(seq)}) {
            const int here = (graph.find(oriented.substr(0, graph.k())) != missing ? 1 : 0) +
                             (graph.find(oriented.substr(oriented.size() - graph.k())) != missing ? 1 : 0);
            best = std::max(best, here);
        }
        return best;
    };
    auto record = [&](const Source &s, size_t nodes_seen, size_t in, size_t out, int connected,
                      const char *outcome) {
        std::ostringstream line;
        line << s.id << '\t' << s.array << '\t' << s.spacer.size() << '\t' << nodes_seen << '\t'
             << in << '\t' << out << '\t' << ends_present(s.spacer) << '\t' << connected << '\t'
             << outcome << '\n';
        return line.str();
    };

    // Each spacer is walked on its own: the graph is only read, and nothing one
    // spacer computes is used by another. The loop therefore runs in parallel,
    // and every line it would have written goes into that spacer's own slot.
    // The files are written afterwards in source order, so the output does not
    // depend on the order the threads happen to finish in.
    struct Result {
        std::string stage, summary;
        std::array<std::string, Count> fasta;
        bool traced = false, candidate = false;
        std::array<size_t, Count> saved{}, thin{};
        std::array<size_t, 5> blocked{};
    };
    std::vector<Result> results(sources.size());
    const long long count = static_cast<long long>(sources.size());
#pragma omp parallel for schedule(dynamic)
    for (long long index = 0; index < count; ++index) {
        const Source &s = sources[static_cast<size_t>(index)];
        Result &r = results[static_cast<size_t>(index)];
        // Per spacer, where the serial version shared one buffer across the loop.
        std::vector<uint64_t> nodes;
        std::vector<double> mult;
        if (!dna(s.spacer) || s.spacer.size() < graph.k()) {
            r.stage = record(s, 0, 0, 0, 0, "not_usable");
            continue;
        }
        std::string spacer = s.spacer;
        bool connected = graph.trace(spacer, nodes, mult);
        if (!connected) {
            spacer = reverse_complement(s.spacer);
            connected = graph.trace(spacer, nodes, mult);
        }
        if (!connected) {
            if (require_connected) {
                r.stage = record(s, 0, 0, 0, 0, "no_connected_path");
                continue;
            }
            spacer = graph.locate(s.spacer, nodes);
            if (nodes.empty()) {
                r.stage = record(s, 0, 0, 0, 0, "no_nodes_in_graph");
                continue;
            }
        }
        r.traced = true;
        const auto &T = tangle.at(s.array);
        // FindProtospacers: outside neighbours over every node of the spacer.
        std::vector<Arc> b_in, b_out;
        const std::unordered_set<uint64_t> own(nodes.begin(), nodes.end());
        for (auto v : nodes)
            for (bool incoming : {true, false})
                for (const auto &arc : graph.neighbours(v, incoming)) {
                    if (T.count(arc.node) || own.count(arc.node)) continue;
                    (incoming ? b_in : b_out).push_back(arc);
                }
        if (b_in.empty() || b_out.empty()) {
            r.stage = record(s, nodes.size(), b_in.size(), b_out.size(), connected,
                             b_in.empty() && b_out.empty()
                                 ? "no_boundary_either_side"
                                 : b_in.empty() ? "no_boundary_in" : "no_boundary_out");
            continue;
        }
        r.stage = record(s, nodes.size(), b_in.size(), b_out.size(), connected, "candidate");
        r.candidate = true;
        for (int criterion = 0; criterion < Count; ++criterion) {
            std::array<Walk, 2> best{};
            std::array<bool, 2> have{false, false};
            for (int side = 0; side < 2; ++side) {
                const bool incoming = side == 0;
                for (const auto &boundary : (incoming ? b_in : b_out)) {
                    const double mu = std::max(boundary.mult, 1e-9);
                    Walk start;
                    start.first = mu;
                    start.lowest = mu;
                    extend_one(start, boundary, criterion, incoming, graph);
                    // The beam returns the best path among those that got
                    // furthest. A path is only compared with another of the same
                    // length, because none of these criteria grows with length
                    // and the shortest path would otherwise always win.
                    auto further = [](const Walk &a, const Walk &b) {
                        return a.bases.size() != b.bases.size() ? a.bases.size() > b.bases.size()
                                                                : a.score > b.score;
                    };
                    std::vector<Walk> frontier{start};
                    Walk terminal = start;
                    while (!frontier.empty()) {
                        std::vector<Walk> next;
                        for (const auto &w : frontier) {
                            if (w.bases.size() >= flank) {
                                if (further(w, terminal)) terminal = w;
                                continue;
                            }
                            bool grew = false;
                            // Keyed either to the boundary node, as the thesis
                            // specifies, or to the median of the walk so far.
                            double centre = mu;
                            if (window_running && !w.multiplicities.empty()) {
                                auto v = w.multiplicities;
                                std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
                                centre = std::max(v[v.size() / 2], 1e-9);
                            }
                            const double low = window_low * centre, high = window_high * centre;
                            const auto around = graph.neighbours(w.nodes.back(), incoming);
                            if (criterion == Harmonic && around.empty()) ++r.blocked[4];
                            for (const auto &arc : around) {
                                if (T.count(arc.node) || own.count(arc.node)) {
                                    if (criterion == Harmonic) ++r.blocked[0];
                                    continue;
                                }
                                if (arc.mult < low) {
                                    if (criterion == Harmonic) ++r.blocked[1];
                                    continue;
                                }
                                if (arc.mult > high) {
                                    if (criterion == Harmonic) ++r.blocked[2];
                                    continue;
                                }
                                if (std::find(w.nodes.begin(), w.nodes.end(), arc.node) != w.nodes.end()) {
                                    if (criterion == Harmonic) ++r.blocked[3];
                                    continue;
                                }
                                Walk child = w;
                                extend_one(child, arc, criterion, incoming, graph);
                                next.push_back(std::move(child));
                                grew = true;
                            }
                            if (!grew && further(w, terminal)) terminal = w;
                        }
                        std::sort(next.begin(), next.end(),
                                  [](const Walk &a, const Walk &b) { return a.score > b.score; });
                        if (next.size() > beam) next.resize(beam);
                        frontier = std::move(next);
                    }
                    if (!have[side] || further(terminal, best[side])) {
                        best[side] = terminal;
                        have[side] = true;
                    }
                }
            }
            if (!have[0] || !have[1]) continue;
            auto left = best[0].bases;
            std::reverse(left.begin(), left.end());
            const auto sequence = left + spacer + best[1].bases;
            // Over the node multiplicities of both sides, so the rule does not
            // depend on which criterion chose the path.
            double reciprocal = 0.;
            size_t counted = 0;
            for (const auto &side : {best[0].multiplicities, best[1].multiplicities})
                for (const double c : side) {
                    reciprocal += 1. / std::max(c, 1e-9);
                    ++counted;
                }
            const double harmonic = counted ? double(counted) / reciprocal : 0.;
            if (harmonic < min_harmonic) {
                ++r.thin[criterion];
                continue;
            }
            ++r.saved[criterion];
            std::ostringstream fasta;
            fasta << '>' << NAMES[criterion] << '_' << s.id << " left=" << left.size()
                  << " right=" << best[1].bases.size() << " spacer=" << spacer.size() << '\n'
                  << sequence << '\n';
            r.fasta[criterion] = fasta.str();
            std::ostringstream row;
            row << NAMES[criterion] << '\t' << s.id << '\t' << s.array << '\t' << nodes.size() << '\t'
                << b_in.size() << '\t' << b_out.size() << '\t' << left.size() << '\t'
                << best[1].bases.size() << '\t' << (best[0].score + best[1].score) / 2. << '\t'
                << harmonic << '\n';
            r.summary += row.str();
        }
    }

    size_t traced = 0, candidates = 0;
    std::array<size_t, Count> saved{}, thin{};
    // Why a step was refused, counted for one criterion so the numbers mean
    // something: in the array, out of the window below, above, already on the
    // walk, and no neighbour in the graph at all.
    std::array<size_t, 5> blocked{};
    for (const auto &r : results) {
        stages << r.stage;
        summary << r.summary;
        for (int c = 0; c < Count; ++c) out[c] << r.fasta[c];
        traced += r.traced ? 1 : 0;
        candidates += r.candidate ? 1 : 0;
        for (int c = 0; c < Count; ++c) {
            saved[c] += r.saved[c];
            thin[c] += r.thin[c];
        }
        for (int i = 0; i < 5; ++i) blocked[i] += r.blocked[i];
    }
    std::cerr << "spacers   " << sources.size() << " given, " << traced << " traced, " << candidates
              << " two-sided\n\n"
              << std::left << std::setw(22) << "criterion" << std::right << std::setw(8) << "kept"
              << std::setw(12) << "below rule" << "\n";
    for (int c = 0; c < Count; ++c)
        std::cerr << std::left << std::setw(22) << NAMES[c] << std::right << std::setw(8) << saved[c]
                  << std::setw(12) << thin[c] << "\n";
    // Ordered by how often each one fires, so the step that is actually
    // limiting the walks is the first line rather than the fourth.
    const char *why[5] = {"in the array or spacer", "below the window", "above the window",
                          "already on the walk", "no neighbour in graph"};
    std::array<int, 5> order{1, 2, 3, 4, 0};
    std::cerr << "\nsteps refused, harmonic_mean only\n";
    for (const int i : order)
        std::cerr << std::left << std::setw(22) << why[i] << std::right << std::setw(20) << blocked[i]
                  << "\n";
    return 0;
}
