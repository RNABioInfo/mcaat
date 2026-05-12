#ifndef NODE_COUNTER_H
#define NODE_COUNTER_H

#include <vector>
#include <utility>
#include <unordered_map>
struct NodeCounter {
    // The struct is used to count all the nodes in the cycle: unordered_map<std::uint64_t, std::vector<std::uint64_t>> systems;
    // INPUT PARAMETERS: input cycles: unordered_map<std::uint64_t, std::vector<std::uint64_t>> systems;
    // OUTPUT PARAMETERS: counts of each node in following format: unordered_map<std::pair<std::uint64_t, std::int>, int> counts;
    // where uint64_t is the node id and first int is the cycle order, second int is the number of occurrences in all the cycles
    #ifdef DEVELOP
    unordered_map<uint64_t, std::vector<uint64_t>> input_cycles;
    unordered_map<int, std::pair<uint64_t, int>> counts;
    
    void count_nodes() {
        for (const auto& cycle : this->input_cycles) {
            const auto& nodes = cycle.second;
            for (size_t i = 0; i < nodes.size(); ++i) {
                uint64_t node_id = nodes[i];
                int order = static_cast<int>(i);
                this->counts[order].first += 1; // Increment occurrence count
                this->counts[order].second = node_id; // Store the order of the node in the cycle
            }
        }
    }
    void filter_only_unique_nodes() {
        unordered_map<int, std::pair<uint64_t, int>> unique_counts;
        for (const auto& count : counts) {
            if (count.second.first == 1) { // Only keep nodes that appear once
                unique_counts[count.first] = count.second;
            }
        }
        counts = unique_counts;
    }
    //constructor to initialize the input cycles
    NodeCounter(const unordered_map<uint64_t, std::vector<uint64_t>>& input_cycles) : input_cycles(input_cycles) {
        count_nodes();
    }
    
    // Function to get the counts of nodes
    const unordered_map<int, std::pair<uint64_t, int>>& get_counts() const {
        return counts;
    }
    #endif
    
};

#endif // NODE_COUNTER_H