#include "cycle_finder.h"
#include "filters.h"
#include "settings.h"
#include <stdexcept>

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
    if(this->visited[neighbor_node]) {
        return false;
    }
    if (repeat_multiplicity / neighbor_node_multiplicity > 500) {
        return false;
    }
    if(original_node == neighbor_node) {
        return false;
    }
    return true;
}


/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity) {
   
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
void CycleFinder::_GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity) {
  
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
void CycleFinder::_GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set) {
   
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
void CycleFinder::_GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set) {
  
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

vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, 
                                        vector<int> backtrack_lengths) {
    int counter = 0;
    uint64_t current_node = start_node;
    vector<vector<uint64_t>> cycles;
    int steps_counter = 0;

    while (!stack.empty()) {
        steps_counter += 1;
        if (steps_counter > 10000000) {
            break;
        }
        
        unordered_set<uint64_t> neighbors = stack.back();
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
                unordered_set<uint64_t> outgoings;
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

                unordered_set<uint64_t> path_set(path.begin(), path.end());

                while (!relax_stack.empty()) {
                    int bl = relax_stack.back().first;
                    int u = relax_stack.back().second;
                    relax_stack.pop_back();
                    if (lock.try_emplace(u, this->settings.cycle_finder_settings.cycle_max_length).first->second < this->settings.cycle_finder_settings.cycle_max_length - bl + 1) {
                        lock[u] = this->settings.cycle_finder_settings.cycle_max_length - bl + 1;
                        unordered_set<uint64_t> incomings;
                        this->_GetIncomings(u, incomings, this->settings.sdbg->EdgeMultiplicity(start_node));
                        for (auto w : incomings)
                            if (path_set.find(w) == path_set.end())
                                relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            }
        }
    }
    
    
    if(cycles.empty()) return {};
    
    #pragma omp critical
    {
        for (const auto& cycle : cycles) 
            for (const auto& node : cycle) 
                this->visited[node] = true;
    }
    
    return cycles;
}


/**
 * @brief Utility function to initialize and start the cycle finding process from a given node.
 */
vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node) {
    vector<uint64_t> path;
    map<uint64_t, int> lock;
    vector<unordered_set<uint64_t>> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    unordered_set<uint64_t> outgoings;
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
    
    unordered_set<uint64_t> tips;

    #pragma omp parallel for
    for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) 
        if(this->settings.sdbg->EdgeOutdegree(node) == 0 && this->settings.sdbg->IsValidEdge(node)) {
            #pragma omp critical
            tips.insert(node);
        }
    return vector<uint64_t>(tips.begin(), tips.end());
}

void CycleFinder::RecursiveReduction(uint64_t tip) {
    if (this->settings.sdbg->EdgeOutdegree(tip)> 0) 
        return;
    unordered_set<uint64_t> parents;
    this->_GetIncomings(tip, parents);
    this->settings.sdbg->SetInvalidEdge(tip);
    for (uint64_t parent : parents) 
        if(this->settings.sdbg->IsValidEdge(parent)) 
            this->RecursiveReduction(parent);
        else 
            continue;
    return;
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
    std::cout << "Pre-filter: invalidated " << invalidated << " node(s) with multiplicity <= 1." << std::endl;
}

/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    uint64_t loaded = 0;
    const int chunk_size = 20000;
    (void)0; // jump_stride removed; placeholder to mark consideration for future optimisation
    if(!this->settings.cycle_finder_settings.low_abundance){
        this->InvalidateMultiplicityOneNodes();
    }
    #pragma omp parallel num_threads(static_cast<int>(this->settings.threads))
    {
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t node = 0; node < this->settings.sdbg->size(); node++) {
            if(!this->settings.sdbg->IsValidEdge(node)) continue;
            size_t edge_indegree = this->settings.sdbg->EdgeIndegree(node);
                // size_t edge_outdegree = this->settings.sdbg->EdgeOutdegree(node); // unused
                loaded+=1; 
                // Provide occasional progress updates during the expensive chunk scan (every 50M nodes)
                if(loaded % 50000000 == 0) std::cout << "ChunkStartNodes: scanned " << (loaded / 1000000) << "M nodes" << std::endl;
                if (edge_indegree >= 2 && this->settings.sdbg->EdgeMultiplicity(node) > this->settings.cycle_finder_settings.threshold_multiplicity)
                {
                    
                    
                    if(this->_IncomingNotEqualToCurrentNode(node,edge_indegree)) continue;
                    int reached_depth = 0;

                    bool dls = this->DepthLevelSearch(node, node, this->settings.cycle_finder_settings.cycle_max_length, reached_depth);
                    if(!dls) continue; //|-> the last version!
            
                    double log2_mult = ceil(log2(double(this->settings.sdbg->EdgeMultiplicity(node))));
                    #pragma omp critical
                    start_nodes_chunked[log2_mult].push_back(node);
            }
        }
    }
   //writeStartNodesToFile(start_nodes_chunked, "start_nodes.txt");
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "Chunked start nodes: multiplicity bucket (log2)=" << key << ", nodes=" << value.size() << std::endl;
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
    std::cout << "Graph size: " << this->settings.sdbg->size() << " nodes; gathered tips: " << tips.size() << std::endl;
    
    this->InvalidateMultiplicityOneNodes();
    for (uint64_t tip : tips) {
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
    std::cout << "After pruning, tips: " << tips.size() << ", valid edges: " << valid_edges << std::endl;
    // struct mallinfo mem_info = mallinfo(); // deprecated
    // size_t graph_mem_info = mem_info.uordblks; // unused
    int cumulative = 0;
    std::cout << "Total nodes in graph: " << this->settings.sdbg->size() << std::endl;
    string mode = "fastq";
    this->visited.resize(this->settings.sdbg->size(), false);
    std::cout << "Starting cycle enumeration: max_len=" << this->settings.cycle_finder_settings.cycle_max_length
              << " min_len=" << this->settings.cycle_finder_settings.cycle_min_length
              << " threads=" << this->settings.threads << std::endl;
        std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> all_cycles;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "Start nodes found in chunks: " << start_nodes_amount << std::endl;
    size_t counter = 0;
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        size_t cumulative_at_bucket_start = cumulative;
        auto thread_count = static_cast<int>(this->settings.threads);
        if (static_cast<int>(nodes_iterator->second.size()) < thread_count)
            thread_count = nodes_iterator->second.size();
        #pragma omp parallel for num_threads(thread_count) reduction(+:cumulative) shared(nodes_iterator, visited)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {
            uint64_t start_node = nodes_iterator->second[start_node_index];
            if (this->visited[start_node]) continue;
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            cumulative += cycles.size();
            this->results[start_node] = cycles;
            ++counter;
        }
        malloc_trim(0);
          // Summarize work performed for this multiplicity bucket and avoid printing too often
        size_t cycles_in_bucket = cumulative - cumulative_at_bucket_start;
          std::cout << "Bucket log2_mult=" << nodes_iterator->first << ": processed " << nodes_iterator->second.size()
                << " nodes, found " << cycles_in_bucket << " cycles (cumulative " << cumulative << ")" << std::endl;
    }
        // Completed cycle enumeration
        std::cout << "Cycle enumeration completed: total cycles=" << cumulative
              << ", result nodes=" << this->results.size() << std::endl;
    return cumulative;
}

CycleFinder::~CycleFinder() {}