#include "crispr_array/unitig_cycle_search.h"
#include <cstdlib>
#include <iostream>
#include <queue>
#include <random>
#include <stdexcept>

struct Graph {
    std::vector<std::vector<uint64_t>> out, in;
    std::vector<bool> valid;
    std::vector<uint64_t> multiplicity;
    explicit Graph(size_t n) : out(n), in(n), valid(n, true), multiplicity(n, 100) {}
    void add(uint64_t u, uint64_t v) { out[u].push_back(v); in[v].push_back(u); }
    bool IsValidEdge(uint64_t v) const { return v < valid.size() && valid[v]; }
    uint64_t EdgeMultiplicity(uint64_t v) const { return multiplicity.at(v); }
    int OutgoingEdges(uint64_t v, uint64_t* a) const {
        if (!IsValidEdge(v)) return -1;
        int n = 0;
        for (auto u : out[v]) if (IsValidEdge(u)) a[n++] = u;
        return n;
    }
    int EdgeOutdegree(uint64_t v) const { uint64_t a[4]; return OutgoingEdges(v, a); }
    int IncomingEdges(uint64_t v, uint64_t* a) const {
        if (!IsValidEdge(v)) return -1;
        int n = 0;
        for (auto u : in[v]) if (IsValidEdge(u)) a[n++] = u;
        return n;
    }
    int EdgeIndegree(uint64_t v) const {
        int n = 0;
        for (auto u : in[v]) n += IsValidEdge(u);
        return n;
    }
};

struct OffsetGraph {
    const Graph& g;
    static constexpr uint64_t base = (uint64_t(1) << 40);
    bool IsValidEdge(uint64_t v) const { return v >= base && g.IsValidEdge(v-base); }
    uint64_t EdgeMultiplicity(uint64_t v) const { return g.EdgeMultiplicity(v-base); }
    int EdgeIndegree(uint64_t v) const { return g.EdgeIndegree(v-base); }
    int EdgeOutdegree(uint64_t v) const { return g.EdgeOutdegree(v-base); }
    int OutgoingEdges(uint64_t v, uint64_t* a) const {
        int n = g.OutgoingEdges(v-base,a);
        for (int i=0;i<n;++i) a[i] += base;
        return n;
    }
    int IncomingEdges(uint64_t v, uint64_t* a) const {
        int n = g.IncomingEdges(v-base,a);
        for (int i=0;i<n;++i) a[i] += base;
        return n;
    }
};

// Independent unit-edge BFS oracle, reading adjacency directly. A shortest
// positive return exists iff any cycle through the start fits the upper bound.
bool oracle(const Graph& g, uint64_t start, int limit) {
    if (limit < 1 || !g.IsValidEdge(start)) return false;
    std::vector<int> distance(g.out.size(), -1);
    std::queue<uint64_t> q;
    q.push(start); distance[start] = 0;
    while (!q.empty()) {
        auto v = q.front(); q.pop();
        if (distance[v] >= limit) continue;
        for (auto u : g.out[v]) {
            if (!g.IsValidEdge(u) || g.multiplicity[u] == 0 ||
                g.multiplicity[start] / g.multiplicity[u] > 500) continue;
            if (u == start) return true;
            if (distance[u] == -1) { distance[u] = distance[v] + 1; q.push(u); }
        }
    }
    return false;
}
uint64_t checks = 0;
void check(bool condition, const char* message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
void compare(Graph& g, int max_limit, size_t cache_capacity = 16384) {
    mcaat::UnitigCycleSearch<Graph> search(g, cache_capacity);
    for (int limit = 0; limit <= max_limit; ++limit)
        for (int round = 0; round < 2; ++round)
            for (uint64_t v = 0; v < g.out.size(); ++v) {
                const bool actual = search.has_cycle(v, limit), expected = oracle(g, v, limit);
                if (actual != expected) {
                    std::cerr << "n="<<g.out.size()<<" start="<<v<<" limit="<<limit
                              <<" actual="<<actual<<" expected="<<expected<<'\n';
                    for (size_t u=0;u<g.out.size();++u) {
                        std::cerr<<u<<" valid="<<g.valid[u]<<" mult="<<g.multiplicity[u]<<":";
                        for (auto w:g.out[u]) std::cerr<<' '<<w;
                        std::cerr<<'\n';
                    }
                }
                check(actual == expected, "cached search disagrees with BFS oracle");
                check(search.cache_size() <= cache_capacity, "unbounded cache");
            }
}

int main() {
    try {
        // All directed four-node graphs, including self-loops and disconnected
        // components. The DNA alphabet admits all degrees in this test.
        for (unsigned mask = 0; mask < (1u << 16); ++mask) {
            Graph g(4);
            for (int u = 0; u < 4; ++u)
                for (int v = 0; v < 4; ++v)
                    if (mask & (1u << (4*u + v))) g.add(u, v);
            compare(g, 5, 3); // also exercises cache eviction
        }
        // Original bounded-DFS counterexample: DFS reaches a shared node by the
        // longer branch first and must reopen it after the shorter arrival.
        Graph revisit(7);
        for (auto e : std::vector<std::pair<int,int>>{
                {0,1},{0,2},{1,4},{2,3},{3,4},{4,5},{5,0},{5,6},{6,0}})
            revisit.add(e.first,e.second);
        compare(revisit, 8);
        // Strict length bounds, long unitigs, pure rings and extra incoming arms.
        for (int length : {1,2,27,76,77,78,79,150}) {
            Graph ring(length + 1);
            for (int i=0; i<length; ++i) ring.add(i,(i+1)%length);
            compare(ring, 80);
            ring.add(length,0);
            compare(ring, 80);
        }
        // Shared repeat + two spacer paths + escape into an unrelated cycle.
        Graph tangle(30);
        for (int i=0;i<12;++i) tangle.add(i,i+1);
        tangle.add(12,0); tangle.add(0,13);
        for (int i=13;i<24;++i) tangle.add(i,i+1);
        tangle.add(24,0); tangle.add(7,25);
        for (int i=25;i<29;++i) tangle.add(i,i+1);
        tangle.add(29,25);
        compare(tangle, 35);
        OffsetGraph large_ids{tangle};
        mcaat::UnitigCycleSearch<OffsetGraph> wide_search(large_ids);
        for (uint64_t v=0;v<tangle.out.size();++v)
            check(wide_search.has_cycle(v+OffsetGraph::base,35) == oracle(tangle,v,35),
                  "64-bit node identifiers");
        check(mcaat::upstream_precheck_bound(77)==78, "odd legacy bound");
        check(mcaat::upstream_precheck_bound(78)==80, "even legacy bound");
        // The cached minimum must be checked against each new start's coverage.
        tangle.multiplicity[0] = 1002; tangle.multiplicity[3] = 2;
        tangle.multiplicity[17] = 0;
        compare(tangle, 35);
        tangle.multiplicity[0] = 1001; // integer quotient 500: permitted
        compare(tangle, 35);
        mcaat::UnitigCycleSearch<Graph> persistent(tangle);
        persistent.has_cycle(0,35);
        tangle.valid[5] = false;
        persistent.clear_cache();
        check(persistent.has_cycle(0,35) == oracle(tangle,0,35), "pruning invalidation");
        compare(tangle, 35, 0); // no retained cache
        // Randomly subdivided DNA-degree graphs make long paths between junctions.
        std::mt19937 rng(1729);
        for (int trial=0; trial<300; ++trial) {
            Graph g(120); int next=6;
            for (int u=0;u<6;++u) for (int v=0;v<6;++v) {
                if (rng()%5 || g.out[u].size()>=4 || g.in[v].size()>=4) continue;
                int prev=u, length=1+rng()%8;
                while (--length && next<120) {g.add(prev,next);prev=next++;}
                g.add(prev,v);
            }
            for (int i=0;i<next;++i) {
                g.multiplicity[i] = rng()%13 == 0 ? 0 : (rng()%10000 + 1);
                if (rng()%23 == 0) g.valid[i] = false;
            }
            compare(g, 25, 7);
        }
        std::cout << "PASS: " << checks << " checks; exhaustive 4-node graphs, "
                     "subdivided graphs, length/coverage boundaries and cache lifecycle\n";
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n'; return 1;
    }
}
