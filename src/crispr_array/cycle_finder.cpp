#include "crispr_array/cycle_finder.h"
#include "crispr_array/filters.h"
#include "core/settings.h"
#include <stdexcept>
#include <fstream>
#include <stack>

// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

/**
 * @file cycle_finder.cpp
 * @brief Implementation of functions for cycle detection and analysis in a sequence graph.
 * - Checking if any incoming edge of a node is not equal to the node itself.
 * - Performing a background check on a neighbor node to determine if it meets certain criteria.
 * - Getting the outgoing edges of a node that pass the background check.
 * - Retrieving the incoming edges of a node that pass the background check.
 * - Finding cycles in the graph.
 * - Performing a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * - Collecting tips (nodes with no outgoing edges) in the graph.
 * - Performing recursive reduction to invalidate edges.
 * - Invalidating nodes with edge multiplicity of one.
 * - Chunking start nodes based on their multiplicity for parallel processing.
 * - Finding all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */


/**
 * @brief Checks if any incoming edge of a node is not equal to the node itself.
 */
bool CycleFinder::_IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree) {
    uint64_t incomings[edge_indegree];
    this->settings.sdbg->IncomingEdges(node, incomings);
    for (const auto& incoming : incomings)
        if (node==incoming)
            return true;
    return false;
}
/**
 * @brief Performs a background check on a neighbor node to determine if it meets certain criteria.
 */
bool CycleFinder::_BackgroundCheck(uint64_t original_node, size_t repeat_multiplicity, uint64_t neighbor_node) {
    auto neighbor_node_multiplicity = this->settings.sdbg->EdgeMultiplicity(neighbor_node);
    // Use lock-free visited test (atomic bitset) to avoid contention
    size_t idx = neighbor_node >> 6;
    uint64_t mask = 1ULL << (neighbor_node & 63);
    int tid = omp_get_thread_num();
    if ((per_thread_visited[tid][idx] & mask) != 0) {
        return false;
    }
    if (repeat_multiplicity / neighbor_node_multiplicity > 500) {
        return false;
    }
    if (original_node == neighbor_node) {
        return false;
    }
    return true; 
}


/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set, size_t repeat_multiplicity) {
   
    int edge_outdegree = this->settings.sdbg->EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
     uint64_t outgoings[edge_outdegree];
    int flag = this->settings.sdbg->OutgoingEdges(node, outgoings);
    if(flag!=-1)    
        for (const auto& outgoing : outgoings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, outgoing) && this->settings.sdbg->IsValidEdge(outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& incomings_set, size_t repeat_multiplicity) {
  
    int edge_indegree = this->settings.sdbg->EdgeIndegree(node);
    if (edge_indegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = this->settings.sdbg->IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, incoming) && this->settings.sdbg->IsValidEdge(incoming))
                incomings_set.insert(incoming);
}
// ## START: HELPER FUNCTIONS FOR DLS ##
/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set) {
   
    int edge_outdegree = this->settings.sdbg->EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t outgoings[edge_outdegree];
    int flag = this->settings.sdbg->OutgoingEdges(node, outgoings);
    if(flag!=-1)
        for (const auto& outgoing : outgoings)
            if (this->settings.sdbg->IsValidEdge(outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& incomings_set) {
  
    int edge_indegree = this->settings.sdbg->EdgeIndegree(node);
    if (edge_indegree == 0 || !this->settings.sdbg->IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = this->settings.sdbg->IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (this->settings.sdbg->IsValidEdge(incoming))
                incomings_set.insert(incoming);
}
// ## END: HELPER FUNCTIONS FOR DLS ##

/**
 * @brief Finds the cycles in the graph using settings configuration.
 * 
 * @param settings The application Settings containing SDBG pointer and cycle finder configuration.
 */
CycleFinder::CycleFinder(Settings& settings)
    : settings(settings), cluster_bounds(500) {
    // Ensure settings.sdbg is valid; this should be set by caller before constructing CycleFinder
    if (this->settings.sdbg == nullptr) {
        throw std::runtime_error("CycleFinder requires settings.sdbg to be set to a valid SDBG instance");
    }
    this->FindApproximateCRISPRArrays();
}

vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t>& path, phmap::flat_hash_map<uint64_t,int>& lock, vector<SmallNodeSet>& stack, 
                                        vector<int>& backtrack_lengths) {
    int counter = 0;
    uint64_t current_node = start_node;
    vector<vector<uint64_t>> cycles;
    int steps_counter = 0;

    while (!stack.empty()) {
        steps_counter += 1;
        if (steps_counter > 10000000) {
            break;
        }
        
        SmallNodeSet neighbors = stack.back();
        bool flag = true;
        for (auto neighbor : neighbors) {
            current_node = neighbor;
            if (current_node == start_node ) {
                backtrack_lengths.back() = 1;
                if (path.size() > static_cast<size_t>(this->settings.cycle_finder_settings.cycle_min_length)) {
                    cycles.push_back(path);
                    counter += 1;
                    if (counter >= this->cluster_bounds) {
                        cycles.clear();
                        flag=false;
                    }
                }
            } 
            else if (static_cast<int>(path.size()) < lock.try_emplace(neighbor, this->settings.cycle_finder_settings.cycle_max_length).first->second) {
                neighbors.erase(neighbor);
                path.push_back(neighbor);
                backtrack_lengths.push_back(this->settings.cycle_finder_settings.cycle_max_length);
                lock[neighbor] = path.size();
                stack.back().erase(neighbor);
                SmallNodeSet outgoings;
                this->_GetOutgoings(neighbor, outgoings, this->settings.sdbg->EdgeMultiplicity(start_node));
                stack.push_back(outgoings);
                flag = false;
                break;
            }
        }
        if (flag) {
            stack.pop_back();
            uint64_t v = path.back();
            path.pop_back();
            int backtrack_length = backtrack_lengths.back();
            backtrack_lengths.pop_back();

            if (!backtrack_lengths.empty()) {
                backtrack_lengths.back() = min(backtrack_lengths.back(), backtrack_length);
            }
            if (backtrack_length < this->settings.cycle_finder_settings.cycle_max_length) {
                vector<pair<int, int>> relax_stack;
                relax_stack.push_back(make_pair(backtrack_length, v));

                phmap::flat_hash_set<uint64_t> path_set(path.begin(), path.end());

                while (!relax_stack.empty()) {
                    int bl = relax_stack.back().first;
                    int u = relax_stack.back().second;
                    relax_stack.pop_back();
                    if (lock.try_emplace(u, this->settings.cycle_finder_settings.cycle_max_length).first->second < this->settings.cycle_finder_settings.cycle_max_length - bl + 1) {
                        lock[u] = this->settings.cycle_finder_settings.cycle_max_length - bl + 1;
                        SmallNodeSet incomings;
                        this->_GetIncomings(u, incomings, this->settings.sdbg->EdgeMultiplicity(start_node));
                        for (auto w : incomings)
                            if (path_set.find(w) == path_set.end())
                                relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            }
        }
    }
    
    
    if (cycles.empty()) return {};

    // Mark visited nodes atomically (no global critical section)
    for (const auto& cycle : cycles)
        for (const auto& node : cycle) {
            size_t idx = node >> 6;
            uint64_t mask = 1ULL << (node & 63);
            int tid = omp_get_thread_num();
            per_thread_visited[tid][idx] |= mask;
        }

    return cycles;
}


/**
 * @brief Utility function to initialize and start the cycle finding process from a given node.
 */
vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node) {
    vector<uint64_t> path;
    phmap::flat_hash_map<uint64_t, int> lock;
    vector<SmallNodeSet> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    SmallNodeSet outgoings;
    this->_GetOutgoings(start_node, outgoings, this->settings.sdbg->EdgeMultiplicity(start_node));
    stack.push_back(outgoings);
    backtrack_lengths.push_back(this->settings.cycle_finder_settings.cycle_max_length);
    return FindCycle(start_node, path, lock, stack, backtrack_lengths);
}
/**
 * @brief Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * Optimized using megahit's graph traversal patterns.
 */
bool CycleFinder::DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) {
    // Megahit-style memory pool: Dynamic thread-local pools for reuse (no heuristic sizing)
    struct StackEntry {
        uint64_t node;
        int depth;
    };
    
    static thread_local std::vector<StackEntry> dls_stack_pool;
    static thread_local phmap::flat_hash_set<uint64_t> dls_visited_pool;
    
    // Clear but keep capacity (megahit memory reuse pattern) - no fixed reserve
    dls_stack_pool.clear();
    dls_visited_pool.clear();
    
    // Reserve reasonable initial capacity to avoid reallocations
    if (dls_stack_pool.capacity() == 0) {
        dls_stack_pool.reserve(64);  // Small initial reserve
    }

    // Use pooled structures
    auto& dls_stack = dls_stack_pool;
    auto& dls_visited = dls_visited_pool;

    // Cache values for faster comparison
    const uint64_t target_node = target;
    const uint64_t start_node = start;

    dls_stack.push_back({start_node, 0});
    reached_depth = 0;

    while (!dls_stack.empty()) {
        StackEntry current = dls_stack.back();
        dls_stack.pop_back();
        uint64_t v = current.node;
        int depth = current.depth;

        // Check if the current node is valid
        if (!this->settings.sdbg->IsValidEdge(v)) {
            continue;
        }

        // Update reached depth
        reached_depth = depth;

        // Get neighbors using SDBG API directly (megahit's pattern) - moved up for better branch prediction
        // Use megahit's efficient zero-degree check for early termination
        if (__builtin_expect(this->settings.sdbg->EdgeOutdegreeZero(v), 0)) {
            continue;
        }
        
        int outdegree = this->settings.sdbg->EdgeOutdegree(v);

        // Use fixed-size array for neighbors (megahit's pattern)
        uint64_t neighbors[MAX_EDGE_COUNT];
        int flag = this->settings.sdbg->OutgoingEdges(v, neighbors);

        if (__builtin_expect(flag == -1, 0)) {
            continue;
        }

        // Prefetch neighbors for better cache performance
        __builtin_prefetch(&neighbors[0], 0, 1);

        // Exceeded depth limit - check after we know we have neighbors
        if (__builtin_expect(depth >= limit, 0)) {
            continue;
        }

        // Process all neighbors to maintain correctness (removed faulty simple path optimization)
        // Process neighbors in forward order (SIMD-friendly pattern from megahit)
        // Unroll loop for small outdegrees to reduce overhead
            // Unrolled loop for common case (de Bruijn graph max degree = 4)
            for (int i = 0; i < outdegree; ++i) {
                uint64_t neighbor = neighbors[i];
                // Check if the neighbor is valid
                if (!this->settings.sdbg->IsValidEdge(neighbor)) {
                    continue;
                }
                auto visited_it = dls_visited.find(neighbor);
                bool not_visited = (visited_it == dls_visited.end());
                bool is_start_revisit = (neighbor == start_node && depth > 0);
                
                if (__builtin_expect(not_visited || is_start_revisit, 1)) {
                    dls_visited.insert(neighbor);
                    dls_stack.push_back({neighbor, depth + 1});
                }
            }
        
        // Megahit-style aggressive early cycle detection
        if (__builtin_expect(v == target_node && depth > 1, 0)) {
            return true;  // Found cycle - exit immediately
        }
    }

    return false;
}


vector<uint64_t> CycleFinder::CollectTips() {
    size_t n = this->settings.sdbg->size();
    int threads = static_cast<int>(this->settings.threads);
    vector<vector<uint64_t>> local_tips(threads);

    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(static)
        for (uint64_t node = 0; node < n; node++) {
            if (this->settings.sdbg->EdgeOutdegree(node) == 0 && this->settings.sdbg->IsValidEdge(node)) {
                local_tips[tid].push_back(node);
            }
        }
    }

    // Merge results while keeping memory usage minimal
    size_t total = 0;
    for (auto &v : local_tips) total += v.size();
    vector<uint64_t> tips;
    tips.reserve(total);
    for (auto &v : local_tips) {
        tips.insert(tips.end(), v.begin(), v.end());
    }
    return tips;
} 

void CycleFinder::RecursiveReduction(uint64_t tip) {
    std::stack<uint64_t> work;
    work.push(tip);
    while (!work.empty()) {
        uint64_t node = work.top();
        work.pop();
        if (this->settings.sdbg->EdgeOutdegree(node) > 0)
            continue;
        if (!this->settings.sdbg->IsValidEdge(node))
            continue;
        SmallNodeSet parents;
        this->_GetIncomings(node, parents);
        this->settings.sdbg->SetInvalidEdge(node);
        for (uint64_t parent : parents)
            if (this->settings.sdbg->IsValidEdge(parent))
                work.push(parent);
    }
}
void CycleFinder::InvalidateMultiplicityOneNodes() {
    uint64_t invalidated = 0;
    #pragma omp parallel for reduction(+:invalidated)
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
        if (this->settings.sdbg->EdgeMultiplicity(node) <=1) {
            this->settings.sdbg->SetInvalidEdge(node);
            invalidated += 1;
        }
    }
    std::cout << "  ▸ Pre-filter: removed " << invalidated << " low-multiplicity edges" << std::endl;
}

/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    const int chunk_size = 20000;
    if(!this->settings.cycle_finder_settings.low_abundance){
        this->InvalidateMultiplicityOneNodes();
    }
    size_t n = this->settings.sdbg->size();
    int threads = static_cast<int>(this->settings.threads);

    // Per-thread chunk collectors to avoid global critical regions
    vector<unordered_map<int, vector<uint64_t>>> local_chunks(threads);
    std::atomic<uint64_t> loaded(0);

    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t node = 0; node < n; node++) {
            if(!this->settings.sdbg->IsValidEdge(node)) continue;
            size_t edge_indegree = this->settings.sdbg->EdgeIndegree(node);
            uint64_t my_loaded = loaded.fetch_add(1, std::memory_order_relaxed) + 1;
            if (my_loaded % 1000000 == 0) std::cout << "  ▸ Scanning: " << (my_loaded / 1000000) << "M nodes..." << std::endl;
            if (edge_indegree >= 2 && this->settings.sdbg->EdgeMultiplicity(node) > this->settings.cycle_finder_settings.threshold_multiplicity)
            {
                if(this->_IncomingNotEqualToCurrentNode(node, edge_indegree)) continue;
                int reached_depth = 0;

                bool dls = this->DepthLevelSearch(node, node, this->settings.cycle_finder_settings.cycle_max_length, reached_depth);
                if(!dls) continue;
                int log2_mult = static_cast<int>(ceil(log2(double(this->settings.sdbg->EdgeMultiplicity(node)))));
                local_chunks[tid][log2_mult].push_back(node);
            }
        }
    }

    // Merge local chunks into the shared map (serial merge to avoid contention)
    for (int t = 0; t < threads; ++t) {
        for (auto &entry : local_chunks[t]) {
            auto &vec = start_nodes_chunked[entry.first];
            vec.insert(vec.end(), entry.second.begin(), entry.second.end());
        }
    }

   //writeStartNodesToFile(start_nodes_chunked, "start_nodes.txt");
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "  ▸ Bucket mult\u00d7" << (1 << key) << " (log2=" << key << "): " << value.size() << " nodes" << std::endl;
        sum_of_all_quantities_in_all_chunks += value.size();
    }
    return sum_of_all_quantities_in_all_chunks;
} 


/**
 * @brief Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */
int CycleFinder::FindApproximateCRISPRArrays()
 {
    
    vector<uint64_t> tips = this->CollectTips();
    std::cout << "  ▸ Graph: " << this->settings.sdbg->size() << " nodes, tips: " << tips.size() << std::endl;
    
    this->InvalidateMultiplicityOneNodes();
    for ( uint64_t tip : tips) {
        this->RecursiveReduction(tip);
    }
    int valid_edges = 0;
    #pragma omp parallel for reduction(+:valid_edges)
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
        if (this->settings.sdbg->IsValidEdge(node)) {
            valid_edges += 1;
        }
    }

    tips = this->CollectTips();
    std::cout << "  ▸ After tip pruning: " << valid_edges << " valid edges, " << tips.size() << " tips remain" << std::endl;
    // struct mallinfo mem_info = mallinfo(); // deprecated
    // size_t graph_mem_info = mem_info.uordblks; // unused
    int cumulative = 0;
    std::cout << "  ▸ Params: max-len=" << this->settings.cycle_finder_settings.cycle_max_length
              << "  min-len=" << this->settings.cycle_finder_settings.cycle_min_length
              << "  threads=" << this->settings.threads << std::endl;
        std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> all_cycles;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "  ▸ Start nodes: " << start_nodes_amount << std::endl;
    size_t counter = 0;
    int max_threads = static_cast<int>(this->settings.threads);
    size_t words = (this->settings.sdbg->size() + 63) / 64;
    per_thread_visited.resize(max_threads);
    for (auto &v : per_thread_visited) {
        v.resize(words, 0);
    }
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        // Reset visited bitmaps so nodes from a higher-multiplicity bucket don't
        // block traversal in lower-multiplicity buckets.
        for (auto &v : per_thread_visited)
            std::fill(v.begin(), v.end(), 0);
        size_t cumulative_at_bucket_start = cumulative;
        auto thread_count = static_cast<int>(this->settings.threads);
        if (static_cast<int>(nodes_iterator->second.size()) < thread_count)
            thread_count = nodes_iterator->second.size();

        // Per-thread results to avoid concurrent modification of shared unordered_map
        vector<unordered_map<uint64_t, vector<vector<uint64_t>>>> local_results(thread_count);
        vector<size_t> local_cycle_counts(thread_count, 0);
        vector<size_t> local_processed(thread_count, 0);

        #pragma omp parallel num_threads(thread_count)
        {
            int tid = omp_get_thread_num();
            #pragma omp for schedule(static)
            for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {
                uint64_t start_node = nodes_iterator->second[start_node_index];
                vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
                local_cycle_counts[tid] += cycles.size();
                local_processed[tid] += 1;
                if (!cycles.empty()) {
                    local_results[tid][start_node] = std::move(cycles);
                }
            }
        }

        // Merge local results into shared structures (serial merge)
        for (int t = 0; t < thread_count; ++t) {
            for (auto &entry : local_results[t]) {
                this->results[entry.first] = std::move(entry.second);
            }
            cumulative += local_cycle_counts[t];
            counter += local_processed[t];
        }

        malloc_trim(0);
        // Summarize work performed for this multiplicity bucket and avoid printing too often
        size_t cycles_in_bucket = cumulative - cumulative_at_bucket_start;
        std::cout << "  ▸ Bucket log2=" << nodes_iterator->first << ": "
                  << nodes_iterator->second.size() << " nodes, "
                  << cycles_in_bucket << " cycles (total " << cumulative << ")" << std::endl;
    }
        // Completed cycle enumeration
        std::cout << "  ▸ Total: " << cumulative << " cycles across "
              << this->results.size() << " start nodes" << std::endl;

    const std::string cycles_path = this->settings.output_folder + "/cycles.txt";
    this->writeMapToFile(this->results, cycles_path);

    return cumulative;
}

void CycleFinder::writeMapToFile(
    const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot write cycles to: " << filename << std::endl;
        return;
    }
    const size_t k = this->settings.sdbg->k();
    std::vector<uint8_t> seq(k);
    for (const auto& [start_node, cycle_vec] : cycles) {
        for (const auto& cycle : cycle_vec) {
            for (size_t i = 0; i < cycle.size(); ++i) {
                seq.assign(k, 0);
                this->settings.sdbg->GetLabel(cycle[i], seq.data());
                std::string label(k, 'N');
                for (size_t j = 0; j < k; ++j) {
                    uint8_t c = seq[k - 1 - j];
                    label[j] = (c >= 1 && c <= 4) ? "ACGT"[c - 1] : 'N';
                }
                if (i > 0) out << ' ';
                out << label;
            }
            out << '\n';
        }
    }
    std::cout << "  ▸ Cycles saved: " << filename << std::endl;
}

CycleFinder::~CycleFinder() {}