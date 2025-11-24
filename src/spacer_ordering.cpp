#include "spacer_ordering.h"

void find_strongly_connected_components_dfs(
    const uint64_t node,
    const SDBG& sdbg,
    unordered_map<uint64_t, int>& index_map,
    unordered_map<uint64_t, int>& lowlink_map,
    unordered_set<uint64_t>& on_stack,
    std::stack<uint64_t>& stack,
    vector<vector<uint64_t>>& components,
    int& index_counter
) {
    index_map[node] = index_counter;
    lowlink_map[node] = index_counter;
    index_counter++;
    stack.push(node);
    on_stack.insert(node);

    int outdegree = sdbg.EdgeOutdegree(node);
    if (outdegree > 0) {
        uint64_t* outgoings = new uint64_t[outdegree];
        if (sdbg.OutgoingEdges(node, outgoings) != -1) {
            for (int i = 0; i < outdegree; ++i) {
                const uint64_t neighbor = outgoings[i];
                if (!sdbg.IsValidEdge(neighbor)) continue;
                
                if (index_map.find(neighbor) == index_map.end()) {
                    find_strongly_connected_components_dfs(neighbor, sdbg, index_map, lowlink_map, on_stack, stack, components, index_counter);
                    lowlink_map[node] = std::min(lowlink_map[node], lowlink_map[neighbor]);
                } else if (on_stack.find(neighbor) != on_stack.end()) {
                    lowlink_map[node] = std::min(lowlink_map[node], index_map[neighbor]);
                }
            }
        }
        delete[] outgoings;
    }

    if (lowlink_map[node] == index_map[node]) {
        vector<uint64_t> component;
        uint64_t w;
        do {
            w = stack.top();
            stack.pop();
            on_stack.erase(w);
            component.push_back(w);
        } while (w != node);

        if (component.size() > 1) {
            components.push_back(component);
        }
    }
}

vector<vector<uint64_t>> find_strongly_connected_components(const SDBG& sdbg) {
    unordered_map<uint64_t, int> index_map;
    unordered_map<uint64_t, int> lowlink_map;
    unordered_set<uint64_t> on_stack;
    std::stack<uint64_t> stack;
    vector<vector<uint64_t>> components;
    int index_counter = 0;

    vector<uint64_t> valid_nodes;
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        if (sdbg.IsValidEdge(node)) {
            valid_nodes.push_back(node);
        }
    }

    for (const uint64_t node : valid_nodes) {
        if (index_map.find(node) == index_map.end()) {
            find_strongly_connected_components_dfs(node, sdbg, index_map, lowlink_map, on_stack, stack, components, index_counter);
        }
    }

    return components;
}

void keep_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const size_t& k,
    const vector<vector<uint64_t>>& cycles
) {
    unordered_set<uint64_t> cycle_edges;
    for (const auto& cycle : cycles) {
        for (uint64_t edge : cycle) {
            cycle_edges.insert(edge);
        }
    }

    // Extend cycle_edges by k
    size_t loop_counter = 0;
    while (loop_counter < k) {
        ++loop_counter;
        vector<uint64_t> new_edges;

        for (const auto& edge : cycle_edges) {
            if (!sdbg.IsValidEdge(edge)) {
                continue;
            }

            int indegree = sdbg.EdgeIndegree(edge);
            if (indegree > 0) {
                uint64_t* incomings = new uint64_t[indegree];
                if (sdbg.IncomingEdges(edge, incomings) != -1) {
                    for (int i = 0; i < indegree; ++i) {
                        const uint64_t neighbor = incomings[i];
                        new_edges.push_back(neighbor);
                    }
                }
                delete[] incomings;
            }

            int outdegree = sdbg.EdgeOutdegree(edge);
            if (outdegree > 0) {
                uint64_t* outgoings = new uint64_t[outdegree];
                if (sdbg.OutgoingEdges(edge, outgoings) != -1) {
                    for (int i = 0; i < outdegree; ++i) {
                        const uint64_t neighbor = outgoings[i];
                        new_edges.push_back(neighbor);
                    }
                }
                delete[] outgoings;
            }
        }

        for (const auto& edge : new_edges) {
            cycle_edges.insert(edge);
        }
    }

    for (uint64_t edge_id = 0; edge_id < sdbg.size(); ++edge_id) {
        if (sdbg.IsValidEdge(edge_id)) {
            // If this edge is not in any cycle, invalidate it
            if (cycle_edges.find(edge_id) == cycle_edges.end()) {
                sdbg.SetInvalidEdge(edge_id);
            }
        }
    }
}

vector<Graph> divide_graph_into_subgraphs(const SDBG& sdbg) {
    vector<Graph> subgraphs;
    
    auto components = find_strongly_connected_components(sdbg);
    
    for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
        const auto& component = components[comp_idx];
        Graph subgraph;
        
        for (const uint64_t edge_id : component) {
            if (!sdbg.IsValidEdge(edge_id)) continue;
            
            int outdegree = sdbg.EdgeOutdegree(edge_id);
            if (outdegree > 0) {
                uint64_t* outgoings = new uint64_t[outdegree];
                if (sdbg.OutgoingEdges(edge_id, outgoings) != -1) {
                    for (int i = 0; i < outdegree; ++i) {
                        const uint64_t neighbor = outgoings[i];
                        // Only add edge if neighbor is in the same component
                        if (std::find(component.begin(), component.end(), neighbor) != component.end()) {
                            subgraph.add_edge(edge_id, neighbor);
                        }
                    }
                }
                delete[] outgoings;
            }
        }
        
        if (subgraph.nodes.size() > 0) {
            subgraphs.push_back(std::move(subgraph));
        }
    }

    return subgraphs;
}

vector<Graph> get_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const size_t& k,
    const vector<vector<uint64_t>>& cycles
) {
    keep_crispr_regions_extended_by_k(sdbg, k, cycles);
    return divide_graph_into_subgraphs(sdbg);
}

vector<vector<uint64_t>> get_relevant_reads(
    const Graph& graph,
    const vector<vector<uint64_t>>& all_reads
) {
    vector<vector<uint64_t>> relevant_reads;

    for (const auto& read : all_reads) {
        if (graph.nodes.find(read.at(0)) != graph.nodes.end()
        || graph.nodes.find(read.at(read.size() - 1)) != graph.nodes.end()) {
            relevant_reads.push_back(read);
        }
    }

    return relevant_reads;
}

vector<vector<uint64_t>> get_relevant_cycles(
    const Graph& graph,
    const vector<vector<uint64_t>>& all_cycles
) {
    vector<vector<uint64_t>> relevant_cycles;

    for (const auto& cycle : all_cycles) {
        bool all_nodes_in_graph = true;
        for (uint64_t node : cycle) {
            if (graph.nodes.find(node) == graph.nodes.end()) {
                all_nodes_in_graph = false;
                break;
            }
        }
        if (all_nodes_in_graph) {
            relevant_cycles.push_back(cycle);
        }
    }

    return relevant_cycles;
}

void get_minimum_cycles_for_full_coverage(vector<vector<uint64_t>>& cycles) {
    if (cycles.size() == 0) {
        return;
    }
    
    unordered_set<uint32_t> universe;
    vector<vector<uint32_t>> sets;

    unordered_map<uint64_t, uint32_t> node_id_map;
    uint32_t problem_element_counter = 0;
    for (const auto& cycle : cycles) {
        vector<uint32_t> set;

        for (const auto& node : cycle) {
            uint32_t mapped_value;
            if (node_id_map.find(node) != node_id_map.end()) {
                mapped_value = node_id_map[node];
            } else {
                node_id_map[node] = problem_element_counter;
                mapped_value = problem_element_counter++;
            }

            set.push_back(mapped_value);
            universe.insert(mapped_value);
        }

        sets.push_back(set);
    }

    if (universe.size() == 0 || sets.size() == 0) {
        return;
    }

    vector<size_t> kept_indices = solve_min_cover_problem(universe, sets);
    std::sort(kept_indices.begin(), kept_indices.end(), std::greater<size_t>());
    
    for (size_t i = 0; i < cycles.size(); ++i) {
        size_t idx = cycles.size() - 1 - i;
        if (std::find(kept_indices.begin(), kept_indices.end(), idx) != kept_indices.end()) {
            continue;
        }

        cycles.erase(cycles.begin() + idx);
    }
}

vector<size_t> solve_min_cover_problem(
    const unordered_set<uint32_t>& universe,
    const vector<vector<uint32_t>>& sets
) {
    if (universe.empty() || sets.empty()) {
        std::cout << "Error: Unable to find min cover as the universe or sets are empty" << std::endl;
        return {};
    }
    for (const auto& element : universe) {
        if (element >= universe.size()) {
            std::cout << "Error: Unable to find min cover as the universe elements are invalid" << std::endl;
            return {};
        }
    }
    for (const auto& set : sets) {
        for (const auto& element : set) {
            if (element >= universe.size()) {
                std::cout << "Error: Unable to find min cover as the sets elements are invalid" << std::endl;
                return {};
            }
        }
    }

    auto set_cover_instance = cft::Instance();

    for (const auto& set : sets) {
        set_cover_instance.cols.push_back(set);
        set_cover_instance.costs.push_back(1.0);
    }

    cft::fill_rows_from_cols(set_cover_instance.cols, universe.size(), set_cover_instance.rows);
    
    auto env = cft::Environment();
    env.time_limit = 10.0;
    env.verbose = 5;
    env.timer.restart();
    env.verbose = false;

    vector<size_t> result;
    if (!sets.empty() && !universe.empty()) {
        const auto solution = cft::run(env, set_cover_instance);
        result = vector<size_t>(solution.sol.idxs.begin(), solution.sol.idxs.end());
    }
    return result;
}

unordered_map<uint64_t, uint32_t> get_node_to_unique_cycle_map(
    const vector<vector<uint64_t>>& cycles
) {
    vector<unordered_set<uint64_t>> cycle_node_sets;
    for (const auto& cycle : cycles) {
        cycle_node_sets.emplace_back(cycle.begin(), cycle.end());
    }

    unordered_map<uint64_t, uint32_t> node_to_cycle_map;
    for (uint32_t i = 0; i < cycle_node_sets.size(); ++i) {
        // Build union of all other sets
        unordered_set<uint64_t> other_nodes;
        for (uint32_t j = 0; j < cycle_node_sets.size(); ++j) {
            if (j == i) continue;
            other_nodes.insert(cycle_node_sets[j].begin(), cycle_node_sets[j].end());
        }
        // Find unique nodes in cycle i
        for (const uint64_t node : cycle_node_sets[i]) {
            if (other_nodes.find(node) == other_nodes.end()) {
                node_to_cycle_map[node] = i;
            }
        }
    }
    return node_to_cycle_map;
}

vector<uint32_t> get_all_cycle_indices(
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<uint32_t> cycle_indices;

    for (const auto& [node, cycle_idx] : node_to_cycle_map) {
        if (std::find(cycle_indices.begin(), cycle_indices.end(), cycle_idx) == cycle_indices.end()) {
            cycle_indices.push_back(cycle_idx);
        }
    }

    return cycle_indices;
}

vector<tuple<uint32_t, uint32_t>> every_possible_combination(const vector<uint32_t>& v) {
    vector<tuple<uint32_t, uint32_t>> possible_combination;
    
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = i + 1; j < v.size(); ++j) {
            const uint32_t element_i = v[i];
            const uint32_t element_j = v[j];

            if (element_i != element_j) {
                tuple<uint32_t, uint32_t> new_tuple = std::make_tuple(element_i, element_j);
                possible_combination.push_back(new_tuple);
            }
        }
    }

    return possible_combination;
}

vector<tuple<uint32_t, uint32_t>> generate_constraints_from_read(
    const Graph& graph,
    const vector<uint64_t>& read,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<uint32_t> cycle_indices_in_order;
    for (const auto& node_id : read) {
        const auto it = node_to_cycle_map.find(node_id);
        if (it != node_to_cycle_map.end()) {
            const uint32_t cycle_index = it->second;
            cycle_indices_in_order.push_back(cycle_index);
        }
    }

    // Merge common neighbors: A,A,B,C,C,C -> A,B,C
    vector<uint32_t> cycle_indices_merged;
    optional<uint32_t> last_cycle_index;
    for (size_t i = 0; i < cycle_indices_in_order.size(); ++i) {
        const uint32_t cycle_index = cycle_indices_in_order.at(i);

        if (i == 0 || !last_cycle_index || cycle_index != *last_cycle_index) {
            cycle_indices_merged.push_back(cycle_index);
            last_cycle_index = std::make_optional(cycle_index);
        }
    }

    return every_possible_combination(cycle_indices_in_order);

    // Only have direct neighbors (through benchmarking it led to worse results)
    // vector<tuple<uint32_t, uint32_t>> derived_constraints;
    // for (int i = 0; i < cycle_indices_merged.size() - 1; ++i) {
    //     const uint32_t cycle_index = cycle_indices_merged.at(i);
    //     const uint32_t cycle_index_next = cycle_indices_merged.at(i + 1);

    //     derived_constraints.push_back(std::make_tuple(cycle_index, cycle_index_next));
    // }

    // return derived_constraints;
}

vector<tuple<uint32_t, uint32_t>> generate_out_of_cycles_constraints_from_read(
    const Graph& graph,
    const vector<uint64_t>& read,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    if (node_to_cycle_map.find(read.at(0)) == node_to_cycle_map.end()
    || node_to_cycle_map.find(read.at(read.size() - 1)) == node_to_cycle_map.end()) {
        return {};
    }

    vector<uint32_t> cycle_indices_in_order;
    for (const auto& node_id : read) {
        const auto it = node_to_cycle_map.find(node_id);

        uint32_t cycle_index;
        if (it == node_to_cycle_map.end()) {
            cycle_index = NOT_IN_ANY_CYCLE_INDEX;
        } else {
            cycle_index = it->second;
        }

        cycle_indices_in_order.push_back(cycle_index);
    }

    // Merge common neighbors: A,A,B,C,C,C -> A,B,C
    vector<uint32_t> cycle_indices_merged;
    optional<uint32_t> last_cycle_index;
    for (size_t i = 0; i < cycle_indices_in_order.size(); ++i) {
        const uint32_t cycle_index = cycle_indices_in_order.at(i);

        if (i == 0 || !last_cycle_index || cycle_index != *last_cycle_index) {
            cycle_indices_merged.push_back(cycle_index);
            last_cycle_index = std::make_optional(cycle_index);
        }
    }

    vector<tuple<uint32_t, uint32_t>> derived_constraints;
    if (cycle_indices_merged.size() > 1) {
        const uint32_t cycle_index = cycle_indices_merged.at(0);
        const uint32_t cycle_index_next = cycle_indices_merged.at(1);

        derived_constraints.push_back(std::make_tuple(cycle_index, cycle_index_next));
    }

    return derived_constraints;
}

vector<tuple<uint32_t, uint32_t>> generate_constraints(
    const Graph& graph,
    const vector<vector<uint64_t>>& reads,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<tuple<uint32_t, uint32_t>> constraints;

    for (const auto& read : reads) {
        const auto read_constraints = generate_constraints_from_read(
            graph,
            read,
            node_to_cycle_map
        );
        const auto extra_constraints = generate_out_of_cycles_constraints_from_read(
            graph,
            read,
            node_to_cycle_map
        );

        for (const auto& constraint : read_constraints) {
            constraints.push_back(constraint);
        }
        for (const auto& constraint : extra_constraints) {
            constraints.push_back(constraint);
        }
    }

    return constraints;
}

vector<tuple<uint32_t, uint32_t>> get_maximal_spanning_tree(
    const vector<tuple<uint32_t, uint32_t>>& edges
) {
    unordered_map<uint32_t, uint32_t> parent;
    unordered_map<uint32_t, int> rank;

    function<uint32_t(uint32_t)> find = [&](uint32_t x) {
        if (parent.find(x) == parent.end()) {
            parent[x] = x;
            rank[x] = 0;
        }
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    };

    auto union_sets = [&](uint32_t x, uint32_t y) -> bool {
        uint32_t root_x = find(x);
        uint32_t root_y = find(y);
        
        if (root_x == root_y) {
            return false; // Would create cycle
        }
        
        if (rank[root_x] < rank[root_y]) {
            parent[root_x] = root_y;
        } else if (rank[root_x] > rank[root_y]) {
            parent[root_y] = root_x;
        } else {
            parent[root_y] = root_x;
            rank[root_x]++;
        }
        return true;
    };


    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash> edges_with_weights;
    for (const auto& edge : edges) {
        edges_with_weights[edge]++;
    }

    vector<pair<int, tuple<uint32_t, uint32_t>>> sorted_edges;
    for (const auto& [edge, weight] : edges_with_weights) {
        sorted_edges.push_back({weight, edge});
    }
    sort(sorted_edges.begin(), sorted_edges.end(), greater<pair<int, tuple<uint32_t, uint32_t>>>());
    
    vector<tuple<uint32_t, uint32_t>> mst_edges;
    for (const auto& [weight, edge] : sorted_edges) {
        const uint32_t from = std::get<0>(edge);
        const uint32_t to = std::get<1>(edge);
        
        if (union_sets(from, to)) {
            mst_edges.push_back(edge);
        }
    }

    return mst_edges;
}

void resolve_cycles_greedy(
    vector<tuple<uint32_t, uint32_t>>& constraints,
    unordered_map<uint32_t, int>& heuristic_node_values
) {
    const auto& mst_edges_vec = get_maximal_spanning_tree(constraints);
    const unordered_set<tuple<uint32_t, uint32_t>, TupleHash> mst_edges(mst_edges_vec.begin(), mst_edges_vec.end());
    
    vector<tuple<uint32_t, uint32_t>> filtered_constraints;
    for (const auto& constraint : constraints) {
        const uint32_t from = std::get<0>(constraint);
        const uint32_t to = std::get<1>(constraint);

        const bool is_removed_edge = mst_edges.find(constraint) == mst_edges.end();
        if (is_removed_edge && from != NOT_IN_ANY_CYCLE_INDEX && to != NOT_IN_ANY_CYCLE_INDEX) {
            heuristic_node_values[to] -= 1;
        } else {
            filtered_constraints.push_back(constraint);
        }
    }

    constraints = std::move(filtered_constraints);
}

void apply_topological_sort(
    vector<uint32_t>& possible_start_nodes,
    const unordered_map<uint32_t, int>& node_affection_to_start,
    unordered_map<uint32_t, int>& heuristic_node_values,
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash>& edges,
    vector<uint32_t>& total_order,
    float& confidence
) {
    if (possible_start_nodes.empty()) {
        return;
    }
    
    // Choose start_node by the heuristic
    int best_start_node = 0;
    float best_heuristic_value = std::numeric_limits<float>::lowest();
    float total_heuristic_value_summed = 0.0;
    
    for (size_t i = 0; i < possible_start_nodes.size(); ++i) {
        const uint32_t node = possible_start_nodes[i];
        const float current_affection = static_cast<float>(node_affection_to_start.at(node));
        const float current_other_heuristic_value = static_cast<float>(heuristic_node_values.at(node));

        // Apply good heuristic weight here
        const float current_heuristic_value = current_affection * 1.0 + current_other_heuristic_value;
        if (current_heuristic_value >= best_heuristic_value) {
            best_heuristic_value = current_heuristic_value;
            best_start_node = i;
        }

        total_heuristic_value_summed += std::abs(current_heuristic_value);
    }

    confidence += (std::abs(best_heuristic_value) / total_heuristic_value_summed);

    uint32_t start_node = possible_start_nodes[best_start_node];
    total_order.push_back(start_node);

    // Remove the start_node as choosable (by index, not value)
    possible_start_nodes.erase(possible_start_nodes.begin() + best_start_node);

    // Removes every edge containing the start_node
    vector<uint32_t> start_node_candidates;
    vector<tuple<uint32_t, uint32_t>> edges_to_remove;
    for (const auto& [edge, weight] : edges) {
        const auto from = std::get<0>(edge);
        const auto to = std::get<1>(edge);

        if (from == start_node) {
            start_node_candidates.push_back(to);
            heuristic_node_values[to] += weight;
            edges_to_remove.push_back(edge);
        }
    }
    for (const auto& edge : edges_to_remove) {
        edges.erase(edge);
    }

    // Finding new possible start nodes
    for (const auto& start_node_candidate : start_node_candidates) {
        bool has_other_incoming = false;
        for (const auto& [edge, _weight] : edges) {
            const auto to = std::get<1>(edge);

            if (to == start_node_candidate) {
                has_other_incoming = true;
                break;
            }
        }
        if (!has_other_incoming) {
            possible_start_nodes.push_back(start_node_candidate);
        }
    }

    apply_topological_sort(
        possible_start_nodes,
        node_affection_to_start,
        heuristic_node_values,
        edges,
        total_order,
        confidence
    );
}

vector<uint32_t> solve_constraints_with_topological_sort(
    const vector<tuple<uint32_t, uint32_t>>& constraints,
    unordered_map<uint32_t, int>& heuristic_node_values,
    const vector<uint32_t>& nodes,
    float& confidence
) {
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash> edges;
    for (const auto& constraint : constraints) {
        const uint32_t source = std::get<0>(constraint);
        const uint32_t dest = std::get<1>(constraint);

        if (source == NOT_IN_ANY_CYCLE_INDEX || dest == NOT_IN_ANY_CYCLE_INDEX) {
            continue;
        }

        edges[constraint]++;
    }

    // Find possible start nodes (nodes with no incoming edges)
    vector<uint32_t> possible_start_nodes;
    for (uint32_t node : nodes) {
        bool has_incoming = false;
        for (const auto& constraint : constraints) {
            const uint32_t source = std::get<0>(constraint);
            const uint32_t dest = std::get<1>(constraint);

            if (source != NOT_IN_ANY_CYCLE_INDEX && dest == node) {
                has_incoming = true;
                break;
            }
        }
        if (!has_incoming) {
            possible_start_nodes.push_back(node);
        }
    }

    // Used for the heuristic, which start node to choose from
    // Big positive number indicates a strong affection towards being the first choosen start node
    // Big negative, affection towards the end
    unordered_map<uint32_t, int> node_affection_to_start;
    for (uint32_t node : nodes) {
        node_affection_to_start[node] = 0;
    }

    for (const auto& constraint : constraints) {
        const uint32_t source = std::get<0>(constraint);
        const uint32_t dest = std::get<1>(constraint);

        if (source != NOT_IN_ANY_CYCLE_INDEX && dest != NOT_IN_ANY_CYCLE_INDEX) {
            continue;
        }

        if (source == NOT_IN_ANY_CYCLE_INDEX) {
            node_affection_to_start[dest]++;
        } else {
            node_affection_to_start[source]--;
        }
    }

    vector<uint32_t> total_order;
    confidence = 0.0;
    apply_topological_sort(
        possible_start_nodes,
        node_affection_to_start,
        heuristic_node_values,
        edges,
        total_order,
        confidence
    );

    confidence /= total_order.size();

    return total_order;
}

vector<uint32_t> order_cycles(
    const Graph& graph,
    const vector<vector<uint64_t>>& reads,
    const vector<vector<uint64_t>>& cycles,
    float& confidence_cycle_resolution,
    float& confidence_topological_sort
) {
    const auto node_to_cycle_map = get_node_to_unique_cycle_map(cycles);
    const auto all_cycle_indices = get_all_cycle_indices(node_to_cycle_map);
    auto constraints = generate_constraints(graph, reads, node_to_cycle_map);

    std::cout << "      ▸ " << constraints.size() << " constraints derived" << std::endl;

    unordered_map<uint32_t, int> heuristic_node_values;
    for (const auto& node : all_cycle_indices) {
        heuristic_node_values[node] = 0;
    }
    const int constraint_count_before = constraints.size();
    resolve_cycles_greedy(constraints, heuristic_node_values);
    confidence_cycle_resolution = static_cast<float>(constraints.size())
        / static_cast<float>(constraint_count_before);

    std::cout << "      ▸ " << constraints.size();
    std::cout << " constraints remain after resolving cycles (confidence = ";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << (confidence_cycle_resolution * 100) << "%)" << std::endl;

    return solve_constraints_with_topological_sort(
        constraints,
        heuristic_node_values,
        all_cycle_indices,
        confidence_topological_sort
    );
}

vector<vector<uint64_t>> get_ordered_cycles(
    const vector<uint32_t>& cycle_order,
    const vector<vector<uint64_t>>& cycles
) {
    vector<vector<uint64_t>> ordered_cycles;

    for (const auto& cycle_index : cycle_order) {
        if (cycle_index < cycles.size()) {
            ordered_cycles.push_back(cycles[cycle_index]);
        }
    }

    return ordered_cycles;
}
