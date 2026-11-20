// Measure, for every reported spacer, how much of its graph depth the source
// array alone can explain.
//
// A spacer's k-mers occur once per array copy. The array's repeat occurs
// (spacers + 1) times per copy. So the repeat depth divided by (spacers + 1)
// estimates the depth one array copy contributes to each of its spacers. Depth
// above that estimate must come from another molecule carrying the same
// sequence, which for a real protospacer is the phage the spacer was acquired
// from. A spacer with no excess has no second copy in the sample, so there is
// no protospacer region to reconstruct.
//
// This uses only graph topology and multiplicity. No reference is read.
//
//   spacer_depth GRAPH_PREFIX SOURCES_TSV > depth.tsv
//
// SOURCES_TSV is the sources.tsv written by phage-curate.
#include "sdbg/sdbg.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
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

std::string label(const SDBG &g, uint64_t node) {
    if (node >= g.size()) return {};
    std::vector<uint8_t> encoded(g.k());
    g.GetLabel(node, encoded.data());
    std::string s;
    for (auto c : encoded) {
        if (c < 1 || c > 4) return {};
        s += ALPHABET[c - 1];
    }
    return s;
}

uint64_t find_node(const SDBG &g, const std::string &word) {
    if (word.size() != size_t(g.k()) || !dna(word)) return missing;
    std::vector<uint8_t> encoded;
    for (char c : word) encoded.push_back(ALPHABET.find(c) + 1);
    const uint64_t n = g.IndexBinarySearch(encoded.data());
    return n != missing && label(g, n) == word ? n : missing;
}

// Multiplicity of every (k+1)-mer of `seq`, following the actual graph path.
// Returns an empty vector when the sequence is not a connected path.
std::vector<double> path_counts(const SDBG &g, const std::string &seq) {
    std::vector<double> counts;
    const size_t k = g.k();
    if (seq.size() <= k || !dna(seq)) return counts;
    uint64_t node = find_node(g, seq.substr(0, k));
    if (node == missing) return counts;
    for (size_t at = k; at < seq.size(); ++at) {
        uint64_t chosen = missing, next = missing;
        uint64_t e = node;
        do {
            auto c = g.GetW(e);
            if (c > 4) c -= 4;
            if (c >= 1 && c <= 4 && g.IsValidEdge(e) && ALPHABET[c - 1] == seq[at]) {
                chosen = e;
                next = g.Forward(e);
                break;
            }
            --e;
        } while (e != missing && !g.IsLastOrTip(e));
        if (chosen == missing) return {};
        counts.push_back(double(g.EdgeMultiplicity(chosen)));
        node = next;
    }
    return counts;
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.;
    const size_t h = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + h, v.end());
    double m = v[h];
    if (v.size() % 2 == 0) m = (m + *std::max_element(v.begin(), v.begin() + h)) / 2.;
    return m;
}

// Best of the two orientations: a reported sequence can sit on either strand.
std::vector<double> oriented_counts(const SDBG &g, const std::string &seq) {
    auto forward = path_counts(g, seq);
    auto reverse = path_counts(g, reverse_complement(seq));
    return forward.size() >= reverse.size() ? forward : reverse;
}

struct Source {
    std::string id, array, spacer, repeat;
};
} // namespace

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "usage: spacer_depth GRAPH_PREFIX SOURCES_TSV\n";
        return 2;
    }
    std::ifstream in(argv[2]);
    if (!in) {
        std::cerr << "cannot read " << argv[2] << "\n";
        return 2;
    }
    std::string line;
    if (!std::getline(in, line)) return 2;
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
    const auto c_id = column("source_id"), c_array = column("array_id");
    const auto c_spacer = column("spacer"), c_repeat = column("repeat");
    std::vector<Source> sources;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::istringstream row(line);
        std::string cell;
        while (std::getline(row, cell, '\t')) f.push_back(cell);
        if (f.size() <= std::max(std::max(c_id, c_array), std::max(c_spacer, c_repeat))) continue;
        sources.push_back({f[c_id], f[c_array], f[c_spacer], f[c_repeat]});
    }
    std::cerr << "read " << sources.size() << " source observations\n";

    SDBG graph;
    graph.LoadFromFile(argv[1]);
    std::cerr << "graph loaded, k=" << graph.k() << ", edges=" << graph.size() << "\n";

    // Spacers per array, and the array's repeat.
    std::map<std::string, std::set<std::string>> array_spacers;
    std::map<std::string, std::string> array_repeat;
    for (const auto &s : sources) {
        array_spacers[s.array].insert(s.spacer);
        if (!s.repeat.empty()) array_repeat[s.array] = s.repeat;
    }
    std::map<std::string, double> repeat_depth;
    for (const auto &[array, repeat] : array_repeat)
        repeat_depth[array] = median(oriented_counts(graph, repeat));

    std::map<std::string, double> spacer_depth;
    for (const auto &s : sources)
        if (!spacer_depth.count(s.spacer)) spacer_depth[s.spacer] = median(oriented_counts(graph, s.spacer));

    std::cout << "source_id\tarray_id\tspacers_in_array\trepeat_depth\tspacer_depth\t"
                 "array_share\texcess\tratio\n";
    for (const auto &s : sources) {
        const double n = double(array_spacers[s.array].size());
        const double rd = repeat_depth.count(s.array) ? repeat_depth[s.array] : 0.;
        const double sd = spacer_depth[s.spacer];
        const double share = rd / (n + 1.);
        std::cout << s.id << '\t' << s.array << '\t' << size_t(n) << '\t' << rd << '\t' << sd << '\t'
                  << share << '\t' << (sd - share) << '\t' << (share > 0. ? sd / share : 0.) << '\n';
    }
    return 0;
}
