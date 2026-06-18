#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H
#include <sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <array>
#include <stack>
#include <omp.h>
#include "core/settings.h"
#include <unordered_set>
//#include "progressbar.hpp"
#include <list>
#include <malloc.h>
#include <bitset>
#include <atomic> 
// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

using namespace std;

struct SmallNodeSet {
    uint64_t nodes[4];
    uint8_t  count = 0;
    void insert(uint64_t n)  { nodes[count++] = n; }
    bool empty()       const { return count == 0; }
    void erase(uint64_t n) {
        for (uint8_t i = 0; i < count; ++i)
            if (nodes[i] == n) { nodes[i] = nodes[--count]; return; }
    }
    uint64_t*       begin()       { return nodes; }
    uint64_t*       end()         { return nodes + count; }
    const uint64_t* begin() const { return nodes; }
    const uint64_t* end()   const { return nodes + count; }
};

class CycleFinder {
    private:
        // Use Settings to configure CycleFinder globally
        Settings& settings;
        static constexpr size_t MAX_EDGE_COUNT = 4;
        // Use SDBG pointer from settings everywhere instead of storing a separate reference
        //SDBG& sdbg;
        uint16_t cluster_bounds;
        // visited bitset stored as 64-bit words (1 bit per node). Use atomic builtins on the words to avoid non-copyable std::atomic in vectors.
        vector<vector<uint64_t>> per_thread_visited;
        // dirty-list: indices of touched 64-bit words in per_thread_visited, per thread.
        // Cleared at the start of each BFS/DLS call; avoids O(N) memset on every call.
        vector<vector<uint32_t>> per_thread_dirty;
        vector<bool> look_up_table;

        // thread count obtained from settings

        //#### DEVELOPER FUNCTIONS ####
        void _WriteStartNodesToFile(const map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _ReadStartNodesFromFile(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _WriteMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename);
        //### DEVELOPER FUNCTIONS ####
        
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####
        bool _IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree);
        bool _BackgroundCheck(uint64_t original_node, size_t repeat_multiplicity, uint64_t current_node);
        void _GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set, size_t repeat_multiplicity);
        void _GetIncomings(uint64_t node, SmallNodeSet& incomings_set, size_t repeat_multiplicity);
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####


        //#### HELPER FUNCTIONS FOR DLS ###
        void _GetOutgoings(uint64_t node, SmallNodeSet& outgoings_set);
        void _GetIncomings(uint64_t node, SmallNodeSet& incomings_set);
        //#### HELPER FUNCTIONS FOR DLS ####

    public:
        // genome/cycles folder available via settings
        CycleFinder(Settings& settings);
        //write a getter for results
        unordered_map<uint64_t, vector<vector<uint64_t>>> results;

        // #### RECURSIVE REDUCTION ####
        vector<uint64_t> CollectTips();
        void RecursiveReduction(uint64_t tip);
        void InvalidateMultiplicityOneNodes();  
        int MultiCycleDepthLevelSearch(uint64_t start, uint64_t target, int limit, int min_cycles);
        void writeMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename);
        //#### DLS ####
        size_t ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked);
        bool DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth);
        bool BidirectionalBFS(uint64_t start, int limit);
        //#### DLS ####
        string CreateFolder();
        //#### CYCLE ENUMERATION ####
        vector<vector<uint64_t>> FindCycle(uint64_t start_node, vector<uint64_t>& path, phmap::flat_hash_map<uint64_t,int>& lock, vector<SmallNodeSet>& stack, vector<int>& backtrack_lengths);
        vector<vector<uint64_t>> FindCycleUtil(uint64_t startnode);
        //#### CYCLE ENUMERATION ####
        void MarkCycleNodesUpTo100();
        // 1. Call ChunkStartNodes to chunk the start nodes based on their multiplicity for parallel processing
        // 1.1 ChunkStartNodes will call DepthLevelSearch to determine if there is a cycle in a certain depth
        // 2. Call FindCycleUtil to find the cycles in the graph
        // 3. Write cycle labels to output_folder/cycles.txt
        // 4. Return the number of cycles found
        int FindApproximateCRISPRArrays();
        ~CycleFinder();
};  
#endif