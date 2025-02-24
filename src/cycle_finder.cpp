#include "cycle_finder.h"
#include "filters.h"
/**
 * @file cycle_finder.cpp
 * @brief Implementation of functions for cycle detection and analysis in a sequence graph.
 *
 * This file contains the implementation of the CycleFinder class, which includes methods for:
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

// @brief Writes the start nodes to a file - developer function
void CycleFinder::_WriteStartNodesToFile(const map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& [key, value] : start_nodes_chunked) {
        outfile << "log2_mult: " << key << " Number of nodes: " << value.size() << endl;
        for (const auto& node : value) {
            outfile << node << endl;
        }
    }

    outfile.close();
}
// @brief Reads the start nodes to a file - developer function
void CycleFinder::_ReadStartNodesFromFile(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    
    std::string line;
    int currentKey = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line.back() == ':') {
            line.pop_back(); // Remove the colon
            currentKey = std::stoi(line);
            start_nodes_chunked[currentKey] = {}; // Initialize the vector
        } else {
            uint64_t value = std::stoull(line);
            start_nodes_chunked[currentKey].push_back(value);
        }
    }
    
    file.close();
    
}
// @brief Writes the cycles to a file - developer function
void CycleFinder::_WriteMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& pair : cycles) {
        uint64_t key = pair.first;
        const std::vector<std::vector<uint64_t>>& nested_vectors = pair.second;

        // Write the key
        outfile << key << ":";

        // Write the nested vectors
        for (size_t i = 0; i < nested_vectors.size(); ++i) {
            const std::vector<uint64_t>& vec = nested_vectors[i];
            
            outfile << "[";
            for (size_t j = 0; j < vec.size(); ++j) {
                outfile << vec[j];
                if (j < vec.size() - 1) {
                    outfile << ",";
                }
            }
            outfile << "]";
            
            if (i < nested_vectors.size() - 1) {
                outfile << ";"; // Separate nested vectors with a semicolon
            }
        }

        outfile << std::endl; // Newline to separate entries
    }

    outfile.close();
}

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
bool CycleFinder::_BackgroundCheck(auto original_node, size_t repeat_multiplicity, auto neighbor_node) {
    auto neighbor_node_multiplicity = sdbg.EdgeMultiplicity(neighbor_node);
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
void CycleFinder::_GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity) {
   
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
     uint64_t outgoings[edge_outdegree];
    int flag =sdbg.OutgoingEdges(node, outgoings);
    if(flag!=-1)    
        for (const auto& outgoing : outgoings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(auto node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity) {
  
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag =sdbg.IncomingEdges(node, incomings);
    if (flag==-1)
        for (const auto& incoming : incomings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, incoming))
                incomings_set.insert(incoming);
}
// ## START: HELPER FUNCTIONS FOR DLS ##
/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set) {
   
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t outgoings[edge_outdegree];
    int flag = sdbg.OutgoingEdges(node, outgoings);
    if(flag!=-1)
        for (const auto& outgoing : outgoings)
            if (sdbg.EdgeMultiplicity(outgoing) > 1)
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(auto node, unordered_set<uint64_t>& incomings_set) {
  
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = sdbg.IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (sdbg.EdgeMultiplicity(incoming) > 1)
                incomings_set.insert(incoming);
}

// ## END: HELPER FUNCTIONS FOR DLS ##

/**
 * @brief Finds the cycles in the graph
 * 
 * @param sdbg The succinct de Bruijn graph.
 * @param length_bound The maximum length of the path to be considered.
 * @param minimal_length The minimum length of the path to be considered.
 * @param genome_name The name of the genome.
 */
CycleFinder::CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name,int threads_count)
    : maximal_length(length_bound), minimal_length(minimal_length), sdbg(sdbg), genome_name(genome_name), cluster_bounds(500), threads_count(threads_count) {
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
                if (path.size() > this->minimal_length) {
                    cycles.push_back(path);
                    counter += 1;
                    if (counter >= this->cluster_bounds) {
                        cycles.clear();
                        flag=false;
                    }
                }
            } 
            else if (path.size() < lock.try_emplace(neighbor, this->maximal_length).first->second) {
                neighbors.erase(neighbor);
                path.push_back(neighbor);
                backtrack_lengths.push_back(this->maximal_length);
                lock[neighbor] = path.size();
                stack.back().erase(neighbor);
                unordered_set<uint64_t> outgoings;
                this->_GetOutgoings(neighbor, outgoings, sdbg.EdgeMultiplicity(start_node));
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
            if (backtrack_length < this->maximal_length) {
                vector<pair<int, int>> relax_stack;
                relax_stack.push_back(make_pair(backtrack_length, v));

                unordered_set<uint64_t> path_set(path.begin(), path.end());

                while (!relax_stack.empty()) {
                    int bl = relax_stack.back().first;
                    int u = relax_stack.back().second;
                    relax_stack.pop_back();
                    if (lock.try_emplace(u, this->maximal_length).first->second < this->maximal_length - bl + 1) {
                        lock[u] = this->maximal_length - bl + 1;
                        unordered_set<uint64_t> incomings;
                        this->_GetIncomings(u, incomings, sdbg.EdgeMultiplicity(start_node));
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
    bool stop = false;
    this->_GetOutgoings(start_node, outgoings, sdbg.EdgeMultiplicity(start_node));
    stack.push_back(outgoings);
    backtrack_lengths.push_back(maximal_length);
    return FindCycle(start_node, path, lock, stack, backtrack_lengths);
}
/**
 * @brief Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.
 */
bool CycleFinder::DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) {
    std::stack<std::pair<uint64_t, int>> dls_stack;
    std::unordered_set<uint64_t> dls_visited;
    dls_stack.emplace(start, 0); // Use emplace for stack operations
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
            this->_GetOutgoings(v, adj);
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

vector<uint64_t> CycleFinder::CollectTips() {
    
    unordered_set<uint64_t> tips;

    #pragma omp parallel for
    for (uint64_t node = 0; node < this->sdbg.size(); node++) 
        if(this->sdbg.EdgeOutdegree(node) == 0) {
            #pragma omp critical
            tips.insert(node);
        }
    return vector<uint64_t>(tips.begin(), tips.end());
}

void CycleFinder::RecursiveReduction(uint64_t tip) {
    if (this->sdbg.EdgeOutdegree(tip)> 0) 
        return;
    unordered_set<uint64_t> parents;
    this->_GetIncomings(tip, parents);
    this->sdbg.SetInvalidEdge(tip);
    for (uint64_t parent : parents) 
        if(this->sdbg.IsValidEdge(parent)) 
            this->RecursiveReduction(parent);
        else 
            continue;
    return;
}
void CycleFinder::InvalidateMultiplicityOneNodes() {
    #pragma omp parallel for
    for (uint64_t node = 0; node < this->sdbg.size(); node++) {
        if (this->sdbg.EdgeMultiplicity(node) == 1) {
            this->sdbg.SetInvalidEdge(node);
        }
    }
}
/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    uint64_t loaded = 0;
    const int chunk_size = 20000;
    this->InvalidateMultiplicityOneNodes();
    #pragma omp parallel num_threads(this->threads_count)
    {
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t node = 0; node < this->sdbg.size(); node++) {
                size_t edge_indegree = this->sdbg.EdgeIndegree(node);
                size_t edge_outdegree = this->sdbg.EdgeOutdegree(node);
                loaded+=1; 
                if(loaded % 10000000 == 0) std::cout << "Loaded " << loaded << " nodes\n";
                if (edge_indegree >= 2 && this->sdbg.EdgeMultiplicity(node) > 19)
                {
                    
                    if(this->_IncomingNotEqualToCurrentNode(node,edge_indegree)) continue;
                    int reached_depth = 0;
                    
                    bool dls = this->DepthLevelSearch(node, node, this->maximal_length, reached_depth);
                    if(!dls) continue; //|-> the last version!
            
                    double log2_mult = ceil(log2(double(this->sdbg.EdgeMultiplicity(node))));
                    start_nodes_chunked[log2_mult].push_back(node);
            }
        }
    }
   //writeStartNodesToFile(start_nodes_chunked, "start_nodes.txt");
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "log2_mult: " << key << " Number of nodes: " << value.size() << endl;
        sum_of_all_quantities_in_all_chunks += value.size();
    }
    return sum_of_all_quantities_in_all_chunks;
}


/**
 * @brief Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */
int CycleFinder::FindApproximateCRISPRArrays()
 {
    /*
    vector<uint64_t> tips = this->CollectTips();
    std::cout<<"BEFORE number of nodes: " << this->sdbg.size() << endl;
    std::cout<< "Number of tips: " << tips.size() << endl;

    for (uint64_t tip : tips) {
        this->RecursiveReduction(tip);
    } 
    int invalid_edges = 0;
    #pragma omp parallel for reduction(+:invalid_edges)
    for (uint64_t node = 0; node < this->sdbg.size(); node++) {
        if (node==-1) {
            invalid_edges += 1;
        }
    }
    */
    //pause
    std::cout << "Recursive Reduction Ignored\n";
    //std::cout<<" Number of invalid edges " << invalid_edges << endl;
    //cin.get();
    struct mallinfo mem_info = mallinfo();
    size_t graph_mem_info = mem_info.uordblks;
    int cumulative = 0;
    printf("Number of nodes in a graph: %lu\n", this->sdbg.size());
    string mode = "fastq";
    this->visited.resize(this->sdbg.size(), false);
    std::cout << "Starting the cycle search:\n\n";
        std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> all_cycles;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "Number of start_nodes: " << start_nodes_amount << endl;
    size_t counter = 0;
    size_t n_th_counter = 0;
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        auto thread_count = this->threads_count;
        if (nodes_iterator->second.size() < thread_count)
            thread_count = nodes_iterator->second.size();
        #pragma omp parallel for num_threads(thread_count) reduction(+:cumulative) shared(nodes_iterator, sdbg, visited)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {
            uint64_t start_node = nodes_iterator->second[start_node_index];
            if (this->visited[start_node]) continue;
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            cumulative += cycles.size();
            this->results[start_node] = cycles;
            ++counter;
            if (counter % 100000 == 0) {
                n_th_counter += 1;
                std::cout << n_th_counter << " 100k nodes\n";
            }
        }
        
        malloc_trim(0);
        std::cout << "processed " << nodes_iterator->second.size() << " with ";
        //nodes_iterator = start_nodes_chunked.erase(nodes_iterator);
        printf("Cycles: %d\n", cumulative);
    }
    std::cout << "Number of cycles: " << cumulative << endl;
    // std::cout number of nodes in results
    std::cout << "Number of nodes in results: " << this->results.size() << endl;
    //writeMapToFile(this->results, "results.txt");

    std::cout << endl;
    std::cout << "Number of Cycles: " << cumulative << endl;
    return cumulative;
}

CycleFinder::~CycleFinder() {}