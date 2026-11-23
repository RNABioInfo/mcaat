// Checks boundary_walk's hand-written neighbour lookup against the library's
// own IncomingEdges and OutgoingEdges, on nodes taken from the spacers of a
// real array report.
//
//   neighbour_check GRAPH_PREFIX SOURCES_TSV [MAX_NODES]
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

bool dna(const std::string &s) {
    return !s.empty() && s.find_first_not_of("ACGT") == std::string::npos;
}

std::string reverse_complement(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char &c : s) {
        const auto i = ALPHABET.find(c);
        c = i == std::string::npos ? c : ALPHABET[3 - i];
    }
    return s;
}
} // namespace

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "usage: neighbour_check GRAPH_PREFIX SOURCES_TSV [MAX_NODES]\n";
        return 2;
    }
    const size_t limit = argc > 3 ? std::stoul(argv[3]) : 4000;

    SDBG g;
    g.LoadFromFile(argv[1]);
    const size_t k = g.k();
    std::cerr << "graph k=" << k << " edges=" << g.size() << "\n";

    auto label = [&](uint64_t node) {
        std::vector<uint8_t> e(k);
        g.GetLabel(node, e.data());
        std::string s;
        for (auto c : e) {
            if (c < 1 || c > 4) return std::string();
            s += ALPHABET[c - 1];
        }
        return s;
    };
    auto find = [&](const std::string &word) -> uint64_t {
        if (word.size() != k || !dna(word)) return missing;
        std::vector<uint8_t> e;
        for (char c : word) e.push_back(ALPHABET.find(c) + 1);
        const uint64_t n = g.IndexBinarySearch(e.data());
        return n != missing && label(n) == word ? n : missing;
    };

    // boundary_walk's version, copied exactly. What it hands the walk is the
    // destination node and the multiplicity, so that is what is compared: the
    // edge id it stores alongside is only used to read the next base.
    auto mine = [&](uint64_t node, bool incoming) {
        std::vector<std::pair<uint64_t, uint64_t>> arcs;
        if (node >= g.size()) return arcs;
        auto add = [&](uint64_t edge, uint64_t to) {
            auto c = g.GetW(edge);
            if (c > 4) c -= 4;
            if (c >= 1 && c <= 4 && g.IsValidEdge(edge))
                arcs.push_back({to, g.EdgeMultiplicity(edge)});
        };
        if (!incoming) {
            uint64_t e = node;
            do {
                auto c = g.GetW(e);
                if (c > 4) c -= 4;
                if (c >= 1 && c <= 4 && g.IsValidEdge(e)) add(e, g.Forward(e));
                --e;
            } while (e != missing && !g.IsLastOrTip(e));
        } else {
            uint64_t first = g.Backward(node);
            if (first >= g.size()) return arcs;
            const auto c = g.GetW(first);
            unsigned ends = g.IsLastOrTip(first);
            add(first, g.GetLastIndex(first));
            for (uint64_t e = first + 1; ends < 5 && e < g.size(); ++e) {
                ends += g.IsLastOrTip(e);
                const auto w = g.GetW(e);
                if (w == c) break;
                if (w == c + 4) add(e, g.GetLastIndex(e));
            }
        }
        std::sort(arcs.begin(), arcs.end());
        return arcs;
    };

    auto theirs = [&](uint64_t node, bool incoming) {
        std::vector<uint64_t> buf(8, 0);
        const int n = incoming ? g.IncomingEdges(node, buf.data())
                               : g.OutgoingEdges(node, buf.data());
        std::vector<std::pair<uint64_t, uint64_t>> arcs;
        for (int i = 0; i < n; ++i) {
            // An edge id is the node it leads to, forward; backward, the node it
            // comes from is the end of that edge's group.
            arcs.push_back({incoming ? uint64_t(g.GetLastIndex(buf[i])) : buf[i],
                            g.EdgeMultiplicity(buf[i])});
        }
        std::sort(arcs.begin(), arcs.end());
        return arcs;
    };

    // Nodes to test: every k-mer of every spacer in the report, both strands.
    std::ifstream in(argv[2]);
    std::string line;
    std::getline(in, line);
    size_t cs = 0;
    {
        std::istringstream h(line);
        std::string f;
        for (size_t i = 0; std::getline(h, f, '\t'); ++i)
            if (f == "spacer") cs = i;
    }
    std::set<uint64_t> nodes;
    while (std::getline(in, line) && nodes.size() < limit) {
        std::istringstream row(line);
        std::vector<std::string> f;
        std::string cell;
        while (std::getline(row, cell, '\t')) f.push_back(cell);
        if (f.size() <= cs) continue;
        for (const auto &oriented : {f[cs], reverse_complement(f[cs])})
            for (size_t i = 0; i + k <= oriented.size(); ++i) {
                const auto n = find(oriented.substr(i, k));
                if (n != missing) nodes.insert(n);
            }
    }
    std::cerr << "comparing " << nodes.size() << " nodes\n";

    size_t same_out = 0, diff_out = 0, same_in = 0, diff_in = 0;
    size_t mine_out = 0, their_out = 0, mine_in = 0, their_in = 0;
    int shown = 0;
    for (const uint64_t n : nodes) {
        const auto mo = mine(n, false), to = theirs(n, false);
        const auto mi = mine(n, true), ti = theirs(n, true);
        mine_out += mo.size();
        their_out += to.size();
        mine_in += mi.size();
        their_in += ti.size();
        if (mo == to) ++same_out; else ++diff_out;
        if (mi == ti) ++same_in; else ++diff_in;
        if ((mo != to || mi != ti) && shown < 5) {
            ++shown;
            std::cerr << "node " << n << " " << label(n) << "\n";
            auto show = [](const char *w, const std::vector<std::pair<uint64_t, uint64_t>> &v) {
                std::cerr << w;
                for (const auto &a : v) std::cerr << "  " << a.first << " x" << a.second;
                std::cerr << "\n";
            };
            show("   out mine:", mo);
            show("   out lib :", to);
            show("   in  mine:", mi);
            show("   in  lib :", ti);
        }
    }
    std::cout << "outgoing: " << same_out << " agree, " << diff_out << " differ"
              << "   (edges: mine " << mine_out << ", library " << their_out << ")\n";
    std::cout << "incoming: " << same_in << " agree, " << diff_in << " differ"
              << "   (edges: mine " << mine_in << ", library " << their_in << ")\n";

    // A step is only a step if the two labels overlap in k-1 bases. This says
    // which edge list describes a real path and which base the step adds.
    size_t ok_out = 0, bad_out = 0, ok_in = 0, bad_in = 0;
    size_t mine_ok_out = 0, mine_bad_out = 0;
    for (const uint64_t n : nodes) {
        const auto here = label(n);
        if (here.empty()) continue;
        for (const auto &a : theirs(n, false)) {
            const auto there = label(a.first);
            (!there.empty() && here.substr(1) == there.substr(0, k - 1) ? ok_out : bad_out)++;
        }
        for (const auto &a : theirs(n, true)) {
            const auto there = label(a.first);
            (!there.empty() && there.substr(1) == here.substr(0, k - 1) ? ok_in : bad_in)++;
        }
        for (const auto &a : mine(n, false)) {
            const auto there = label(a.first);
            (!there.empty() && here.substr(1) == there.substr(0, k - 1) ? mine_ok_out : mine_bad_out)++;
        }
    }
    std::cout << "library outgoing steps overlapping in k-1 bases: " << ok_out << " of "
              << ok_out + bad_out << "\n";
    std::cout << "library incoming steps overlapping in k-1 bases: " << ok_in << " of "
              << ok_in + bad_in << "\n";
    std::cout << "mine    outgoing steps overlapping in k-1 bases: " << mine_ok_out << " of "
              << mine_ok_out + mine_bad_out << "\n";

    // How the outgoing lists differ, split by kind.
    size_t extra = 0, absent = 0, wrong_mult = 0, nodes_extra = 0, nodes_mult = 0;
    double worst = 1.;
    for (const uint64_t n : nodes) {
        const auto mo = mine(n, false), to = theirs(n, false);
        std::set<uint64_t> md, td;
        std::map<uint64_t, uint64_t> mm, tm;
        for (const auto &a : mo) { md.insert(a.first); mm[a.first] = a.second; }
        for (const auto &a : to) { td.insert(a.first); tm[a.first] = a.second; }
        bool had_extra = false, had_mult = false;
        for (const auto d : md)
            if (!td.count(d)) { ++extra; had_extra = true; }
        for (const auto d : td)
            if (!md.count(d)) ++absent;
        for (const auto d : md)
            if (td.count(d) && mm[d] != tm[d]) {
                ++wrong_mult;
                had_mult = true;
                const double r = double(mm[d]) / double(std::max<uint64_t>(tm[d], 1));
                if (std::max(r, 1. / r) > worst) worst = std::max(r, 1. / r);
            }
        nodes_extra += had_extra;
        nodes_mult += had_mult;
    }
    std::cout << "\nhow the outgoing lists differ, over " << nodes.size() << " nodes\n";
    std::cout << "   destinations mine has and the library does not: " << extra
              << "  (on " << nodes_extra << " nodes)\n";
    std::cout << "   destinations the library has and mine does not: " << absent << "\n";
    std::cout << "   same destination, different multiplicity:       " << wrong_mult
              << "  (on " << nodes_mult << " nodes, worst ratio " << worst << ")\n";

    // Is GetW enough to name the base a step adds, or is the label needed?
    size_t getw_ok = 0, getw_bad = 0;
    for (const uint64_t n : nodes)
        for (const auto &a : theirs(n, false)) {
            const auto there = label(a.first);
            auto c = g.GetW(a.first);
            if (c > 4) c -= 4;
            const char from_w = c >= 1 && c <= 4 ? ALPHABET[c - 1] : 'N';
            (!there.empty() && from_w == there.back() ? getw_ok : getw_bad)++;
        }
    std::cout << "   GetW(edge) names the added base correctly: " << getw_ok << " of "
              << getw_ok + getw_bad << "\n";

    // boundary_walk writes the base as ALPHABET[GetW(sibling edge)] and the
    // destination as Forward(sibling edge). Does that pair agree with the
    // destination's own label, which is what the sequence has to spell?
    size_t base_ok = 0, base_bad = 0;
    for (const uint64_t n : nodes) {
        if (n >= g.size()) continue;
        uint64_t e = n;
        do {
            auto c = g.GetW(e);
            if (c > 4) c -= 4;
            if (c >= 1 && c <= 4 && g.IsValidEdge(e)) {
                const auto there = label(g.Forward(e));
                const char written = ALPHABET[c - 1];
                (!there.empty() && written == there.back() ? base_ok : base_bad)++;
            }
            --e;
        } while (e != missing && !g.IsLastOrTip(e));
    }
    std::cout << "   the base boundary_walk writes matches the node it moves to: " << base_ok
              << " of " << base_ok + base_bad << "\n";

    // And the library's edge id, used as both destination and base source.
    size_t lib_ok = 0, lib_bad = 0;
    for (const uint64_t n : nodes) {
        const auto here = label(n);
        if (here.empty()) continue;
        for (const auto &a : theirs(n, false)) {
            const auto there = label(a.first);
            (!there.empty() && here.substr(1) + there.back() == there ? lib_ok : lib_bad)++;
        }
    }
    std::cout << "   library: label(node).substr(1) + label(dest).back() == label(dest): "
              << lib_ok << " of " << lib_ok + lib_bad << "\n";
    return 0;
}
