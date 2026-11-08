#ifndef MCAAT_UNITIG_CYCLE_SEARCH_H
#define MCAAT_UNITIG_CYCLE_SEARCH_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>
#include <parallel_hashmap/phmap.h>

namespace mcaat {

// The existing bidirectional precheck expands floor(limit/2)+1 steps in EACH
// direction. Keep that effective bound at the integration boundary so enabling
// the optimization does not also change which candidates MCAAT accepts.
// has_cycle() itself always obeys the exact bound supplied to it.
inline int upstream_precheck_bound(int limit) {
    const int64_t bound = 2 * (int64_t(limit) / 2 + 1);
    return static_cast<int>(std::min<int64_t>(bound, std::numeric_limits<int>::max()));
}

// Exact bounded cycle-existence precheck on the SDBG edge-as-node interface.
// Own one instance per worker and graph-pruning epoch. The graph must not change
// while the cache is live; clear_cache() is required after any graph mutation.
// This is a predicate only: no edges, sequences, parallel paths or cycle output
// are rewritten. Lengths count original directed transitions, not cached jumps.
template<class Graph>
class UnitigCycleSearch {
public:
    struct Statistics {
        uint64_t queries = 0, expansions = 0, cache_hits = 0, cache_misses = 0;
        uint64_t built_steps = 0, cached_steps = 0, cache_resets = 0;
    };

    explicit UnitigCycleSearch(const Graph& graph, size_t cache_capacity = 16384)
        : graph_(graph), capacity_(cache_capacity) {}

    bool has_cycle(uint64_t start, int limit, bool use_jumps = true) {
        ++stats_.queries;
        if (limit <= 0 || !graph_.IsValidEdge(start)) return false;
        if (limit != cache_limit_) {
            clear_cache();
            cache_limit_ = limit;
        }
        const uint64_t start_multiplicity = graph_.EdgeMultiplicity(start);
        if (start_multiplicity == 0) return false;
        // A cached chain must not skip the query start. MCAAT candidates have
        // indegree >= 2, so are always boundaries. For general callers querying
        // inside a unitig, use the same search with single-step edges instead.
        use_jumps = use_jumps &&
            (graph_.EdgeIndegree(start) != 1 || graph_.EdgeOutdegree(start) != 1);
        for (int direction = 0; direction < 2; ++direction) {
            best_depth_[direction].clear();
            frontier_[direction].clear();
            best_depth_[direction][start] = 0;
            frontier_[direction].push_back({start, 0});
        }
        const auto max_depth = static_cast<uint32_t>(limit);
        // Weighted bidirectional search: a jump costs its original hop count.
        // An intersection is accepted only if the two distances fit the bound.
        // The initial zero-length intersection at start is deliberately ignored.
        while (true) {
            discard_stale(0);
            discard_stale(1);
            if (frontier_[0].empty() || frontier_[1].empty()) return false;
            const auto forward_depth = frontier_[0].front().depth;
            const auto backward_depth = frontier_[1].front().depth;
            if (uint64_t(forward_depth) + backward_depth > max_depth) return false;
            const int direction = forward_depth <= backward_depth ? 0 : 1;
            auto& frontier = frontier_[direction];
            auto& best = best_depth_[direction];
            auto& other = best_depth_[1-direction];
            const State state = frontier.front();
            std::pop_heap(frontier.begin(), frontier.end(), Later{});
            frontier.pop_back();
            ++stats_.expansions;
            uint64_t neighbors[4];
            const int degree = direction ? graph_.IncomingEdges(state.node, neighbors)
                                         : graph_.OutgoingEdges(state.node, neighbors);
            for (int i = 0; i < degree; ++i) {
                const auto next = neighbors[i];
                if (!graph_.IsValidEdge(next)) continue;
                const Jump jump = use_jumps ? advance(next, direction) :
                    Jump{next, graph_.EdgeMultiplicity(next), 1};
                if (jump.length > max_depth - state.depth) continue;
                // Preserve the existing integer-division multiplicity filter,
                // including checks on every skipped internal node, not only ends.
                if (jump.min_multiplicity == 0 ||
                    start_multiplicity / jump.min_multiplicity > 500) continue;
                if (jump.end == start) return true;
                const uint32_t depth = state.depth + jump.length;
                const auto meeting = other.find(jump.end);
                if (meeting != other.end() && uint64_t(depth) + meeting->second <= max_depth)
                    return true;
                if (depth >= max_depth) continue;
                const auto old = best.find(jump.end);
                if (old == best.end() || depth < old->second) {
                    best[jump.end] = depth;
                    frontier.push_back({jump.end, depth});
                    std::push_heap(frontier.begin(), frontier.end(), Later{});
                }
            }
        }
        return false;
    }

    void clear_cache() { cache_[0].clear(); cache_[1].clear(); }
    void reset_statistics() { stats_ = {}; }
    const Statistics& statistics() const { return stats_; }
    size_t cache_size() const { return cache_[0].size() + cache_[1].size(); }

private:
    struct State { uint64_t node; uint32_t depth; };
    struct Later {
        bool operator()(const State& a, const State& b) const {
            return a.depth > b.depth || (a.depth == b.depth && a.node > b.node);
        }
    };
    void discard_stale(int direction) {
        auto& q = frontier_[direction];
        while (!q.empty() && best_depth_[direction].at(q.front().node) != q.front().depth) {
            std::pop_heap(q.begin(), q.end(), Later{});
            q.pop_back();
        }
    }
    // The key is the first node after a boundary. Length includes the edge into
    // that node and all following forced steps. Min includes the final endpoint.
    struct Jump { uint64_t end, min_multiplicity; uint32_t length; };

    Jump advance(uint64_t first, int direction) {
        auto& cache = cache_[direction];
        const auto found = cache.find(first);
        if (found != cache.end()) {
            ++stats_.cache_hits;
            stats_.cached_steps += found->second.length - 1;
            return found->second;
        }
        ++stats_.cache_misses;
        Jump jump{first, graph_.EdgeMultiplicity(first), 1};
        // Do not walk a long genomic unitig just to reject a short cycle query.
        // A prefix of limit+1 steps already cannot close an eligible cycle.
        const uint32_t cap = static_cast<uint32_t>(cache_limit_) + 1;
        while (jump.length < cap) {
            const int opposite_degree = direction ? graph_.EdgeOutdegree(jump.end)
                                                  : graph_.EdgeIndegree(jump.end);
            if (opposite_degree != 1) break;
            uint64_t next[4];
            const int degree = direction ? graph_.IncomingEdges(jump.end, next)
                                         : graph_.OutgoingEdges(jump.end, next);
            if (degree != 1) break;
            if (!graph_.IsValidEdge(next[0])) break;
            jump.end = next[0];
            jump.min_multiplicity = std::min(
                jump.min_multiplicity, uint64_t(graph_.EdgeMultiplicity(jump.end)));
            ++jump.length;
            ++stats_.built_steps;
        }
        if (capacity_ != 0) {
            if (cache_size() >= capacity_) {
                clear_cache();
                ++stats_.cache_resets;
            }
            cache.emplace(first, jump);
        }
        return jump;
    }

    const Graph& graph_;
    const size_t capacity_;
    int cache_limit_ = -1;
    phmap::flat_hash_map<uint64_t, Jump> cache_[2];
    phmap::flat_hash_map<uint64_t, uint32_t> best_depth_[2];
    std::vector<State> frontier_[2];
    Statistics stats_;
};
} // namespace mcaat
#endif
