#include "cycle_finder.h"
/**
 * @brief Checks if any incoming edge of a node is not equal to the node itself.
 */
bool CycleFinder::_IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree) {
    uint64_t incomings[edge_indegree];
    this->sdbg.IncomingEdges(node, incomings);
    for (const auto& incoming : incomings)
        if (node==incoming)
            return true;
    return false;
}
/**
 * @brief Performs a background check on a neighbor node to determine if it meets certain criteria.
 */
bool CycleFinder::_BackgroundCheck(auto original_node, size_t repeat_multiplicity, auto neighbor_node, bool& stop) {
    size_t neighbor_node_multiplicity = this->sdbg.EdgeMultiplicity(neighbor_node);
    bool base_case = neighbor_node_multiplicity > 1;

    return base_case;
}


/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity,bool& stop) {
   
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || node == -1) {
        outgoings_set.clear();
        return;
    }
    uint64_t outgoings[edge_outdegree];
    sdbg.OutgoingEdges(node, outgoings);
    for (const auto& outgoing : outgoings)
        if (this->_BackgroundCheck(node, repeat_multiplicity, outgoing,stop))
            outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(auto node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity,bool& stop) {
  
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || node == -1) {
        incomings_set.clear();
        return;
    }
    uint64_t incomings[edge_indegree];
    sdbg.IncomingEdges(node, incomings);
    for (const auto& incoming : incomings)
        if (this->_BackgroundCheck(node, repeat_multiplicity, incoming, stop))
            incomings_set.insert(incoming);
}


/**
 * @brief Finds the cycles in the graph
 * 
 * @param sdbg The succinct de Bruijn graph.
 * @param length_bound The maximum length of the path to be considered.
 * @param minimal_length The minimum length of the path to be considered.
 * @param genome_name The name of the genome.
 */
CycleFinder::CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name)
    : maximal_length(length_bound), minimal_length(minimal_length), sdbg(sdbg), genome_name(genome_name), cluster_bounds(500) {
    this->FindCycles();
}


/**
 * @brief Relaxes the lock values for nodes during backtracking to allow further exploration
 */
void CycleFinder::RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, size_t start_node, bool& stop) {  
    vector<pair<int, int>> relax_stack;
    relax_stack.push_back(make_pair(backtrack_length, v));

    while (!relax_stack.empty()) {
        int bl = relax_stack.back().first;
        int u = relax_stack.back().second;
        relax_stack.pop_back();
        if (lock.try_emplace(u, this->maximal_length).first->second < this->maximal_length - bl + 1) {
            lock[u] = this->maximal_length - bl + 1;
            unordered_set<uint64_t> incomings;
            this->_GetIncomings(u, incomings, sdbg.EdgeMultiplicity(start_node),stop);
            for (auto w : incomings)
                if (find(path.begin(), path.end(), w) == path.end())
                    relax_stack.push_back(make_pair(bl + 1, w));
        }
    }
}

/**
 * @brief Processes the neighbors of the current node and updates the path and stack accordingly.
 */
bool CycleFinder::ProcessNeighbors(unordered_set<uint64_t>& neighbors, int& current_node, uint64_t start_node, vector<uint64_t>& path, map<uint64_t, int>& lock, 
                                    vector<unordered_set<uint64_t>>& stack, vector<int>& backtrack_lengths, 
                                    vector<vector<uint64_t>>& cycles, int& counter,bool& stop) {
    bool flag = true;
    for (auto neighbor : neighbors) {
        current_node = neighbor;
        if (current_node == start_node ) {
            backtrack_lengths.back() = 1;
            if (path.size() > this->minimal_length) {
                cycles.push_back(path);
                counter += 1;
                if (counter >= this->cluster_bounds) {
                    cycles.clear();
                    return false;
                }
            }
        } else if (path.size() < lock.try_emplace(neighbor, this->maximal_length).first->second) {
            neighbors.erase(neighbor);
            path.push_back(neighbor);
            backtrack_lengths.push_back(this->maximal_length);
            lock[neighbor] = path.size();
            stack.back().erase(neighbor);
            unordered_set<uint64_t> outgoings;
            this->_GetOutgoings(neighbor, outgoings, sdbg.EdgeMultiplicity(start_node),stop);
            stack.push_back(outgoings);
            flag = false;
            break;
        }
    }
    return flag;
}
/**
 * @brief Performs backtracking during the cycle finding process to explore alternative paths
 */
void CycleFinder::Backtrack(vector<unordered_set<uint64_t>>& stack, vector<uint64_t>& path, vector<int>& backtrack_lengths, map<uint64_t, int>& lock, uint64_t start_node,bool& stop) {
    stack.pop_back();
    uint64_t v = path.back();
    path.pop_back();
    int backtrack_length = backtrack_lengths.back();
    backtrack_lengths.pop_back();

    if (!backtrack_lengths.empty()) {
        backtrack_lengths.back() = min(backtrack_lengths.back(), backtrack_length);
    }
    if (backtrack_length < this->maximal_length) {
        this->RelaxLock(backtrack_length, v, path, lock, start_node,stop);
    }
}

/**
 * @brief Finds cycles starting from a given node using a depth-first search approach with backtracking.
 */
vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, 
                                        vector<int> backtrack_lengths,bool& stop) {
    int counter = 0;
    int current_node = start_node;
    vector<vector<uint64_t>> cycles;
    int steps_counter = 0;

    while (!stack.empty()) {
        steps_counter += 1;
        if (stop) return {};
        if (steps_counter > 10000000) {
            return cycles;
        }
        
        unordered_set<uint64_t> neighbors = stack.back();
        bool flag = ProcessNeighbors(neighbors, current_node, start_node, path, lock, stack, backtrack_lengths, cycles, counter,stop);
        
        if (flag) {
            Backtrack(stack, path, backtrack_lengths, lock, start_node,stop);
        }
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
    bool stop = false;
    this->_GetOutgoings(start_node, outgoings, sdbg.EdgeMultiplicity(start_node),stop);
    stack.push_back(outgoings);
    backtrack_lengths.push_back(maximal_length);
    return FindCycle(start_node, path, lock, stack, backtrack_lengths,stop);
}
/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    int loaded = 0;
    #pragma omp parallel for num_threads(24)
    for (uint64_t node = 0; node < this->sdbg.size(); node++) {
        size_t edge_indegree = this->sdbg.EdgeIndegree(node);
        size_t edge_outdegree = this->sdbg.EdgeOutdegree(node);
        if (edge_indegree >= 2 && this->sdbg.EdgeMultiplicity(node) > 19)
        {
            loaded+=1; 
            if(loaded % 1000000 == 0) std::cout << "Loaded " << loaded << " nodes\n";
            if(this->_IncomingNotEqualToCurrentNode(node,edge_indegree)) continue;
            int reached_depth = 0;
            bool dls = this->DepthLevelSearch(node, node, 100, reached_depth);
            //if (!not_dls) std::cout << "Node " << node<< ", reached depth " << reached_depth<< endl;
            if(!dls) continue; //-> the last version!
            double log2_mult = ceil(log2(double(this->sdbg.EdgeMultiplicity(node))));
            start_nodes_chunked[log2_mult].push_back(node);
        }
    }
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "log2_mult: " << key << " Number of nodes: " << value.size() << endl;
        sum_of_all_quantities_in_all_chunks += value.size();
    }
    return sum_of_all_quantities_in_all_chunks;
}
/**
 * @brief Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.
 */
bool CycleFinder::DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) {
    std::stack<std::pair<uint64_t, int>> dls_stack;
    std::unordered_set<uint64_t> dls_visited;
    dls_stack.emplace(start, 0); // Use emplace for stack operations
    bool stop = false;

    while (!dls_stack.empty()) {
        auto [v, depth] = dls_stack.top();
        dls_stack.pop();

        if (v == target && depth > 1) {
            return true;
        }

        reached_depth = depth;

        if (depth < limit) {
            std::unordered_set<uint64_t> adj;
            adj.reserve(kAlphabetSize + 1); // Preallocate memory for adjacency list
            this->_GetOutgoings(v, adj, sdbg.EdgeMultiplicity(start), stop);
            if (stop) return false;
            for (uint64_t neighbor : adj) {
                if (dls_visited.find(neighbor) == dls_visited.end() || (neighbor == start && depth > 0)) {
                    dls_stack.emplace(neighbor, depth + 1); // Use emplace for stack operations
                    dls_visited.insert(neighbor);
                }
            }
        }
    }
    return false;
}

/**
 * @brief Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */
int CycleFinder::FindCycles() {
    struct mallinfo mem_info = mallinfo();
    size_t graph_mem_info = mem_info.uordblks;
    int cumulative = 0;
    printf("Number of nodes in a graph: %lu\n", this->sdbg.size());
    string mode = "fastq";
    this->visited.resize(this->sdbg.size(), false);
    std::cout << "Starting the cycle search:\n\n";
    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount = ChunkStartNodes(start_nodes_chunked);
    std::cout << "Number of start_nodes: " << start_nodes_amount << endl;
    size_t counter = 0;
    size_t n_th_counter = 0;
   //std::cin.get();
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        auto thread_count = 24;
        if (nodes_iterator->second.size() < thread_count)
            thread_count = nodes_iterator->second.size();
        #pragma omp parallel for num_threads(thread_count) reduction(+:cumulative) shared(nodes_iterator, sdbg, visited)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {


            uint64_t start_node = nodes_iterator->second[start_node_index];
            if (this->visited[start_node]) continue;
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            cumulative += cycles.size();
            #pragma omp critical
            {
                for (const auto& cycle : cycles) {
                    PathWriter("a", this->sdbg, cycle, this->genome_name, "fastq", this->visited);
                }
            }
            ++counter;
            if (counter % 100000 == 0) {
                n_th_counter += 1;
                std::cout << n_th_counter << " 100k nodes\n";
            }
        }
        
        malloc_trim(0);
        std::cout << "processed " << nodes_iterator->second.size() << " with ";
        //nodes_iterator = start_nodes_chunked.erase(nodes_iterator);
        printf("number of Cycles: %d\n", cumulative);
    }
    std::cout << endl;
    std::cout << "Number of Cycles: " << cumulative << endl;
    return cumulative;
}

CycleFinder::~CycleFinder() {}