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
#include "path_writer.h"
#include "settings.h"
#include <unordered_set>
//#include "progressbar.hpp"
#include <list>
#include <malloc.h>
#include <bitset>
#include <atomic> 
// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

using namespace std;


class CycleFinder {
    private:
        // Use Settings to configure CycleFinder globally
        Settings& settings;
        static constexpr size_t MAX_EDGE_COUNT = 4;
        // Use SDBG pointer from settings everywhere instead of storing a separate reference
        //SDBG& sdbg;
        uint16_t cluster_bounds;
        // visited bitset stored as 64-bit words (1 bit per node). Use atomic builtins on the words to avoid non-copyable std::atomic in vectors.
        std::vector<uint64_t> visited_words;
        vector<bool> look_up_table;

        // helpers for visited bitset (lock-free via gcc/clang atomic builtins)
        inline void InitializeVisited(size_t n) {
            size_t words = (n + 63) / 64;
            this->visited_words.resize(words);
            for (auto &w : this->visited_words) w = 0;
        }
        inline bool IsVisited(uint64_t node) const {
            size_t idx = node >> 6;
            uint64_t mask = 1ULL << (node & 63);
            return (__atomic_load_n(&this->visited_words[idx], __ATOMIC_RELAXED) & mask) != 0;
        }
        inline void MarkVisited(uint64_t node) {
            size_t idx = node >> 6;
            uint64_t mask = 1ULL << (node & 63);
            __atomic_fetch_or(&this->visited_words[idx], mask, __ATOMIC_RELAXED);
        }
        // thread count obtained from settings

        //#### DEVELOPER FUNCTIONS ####
        void _WriteStartNodesToFile(const map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _ReadStartNodesFromFile(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _WriteMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename);
        //### DEVELOPER FUNCTIONS ####
        
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####
        bool _IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree);
        bool _BackgroundCheck(uint64_t original_node, size_t repeat_multiplicity, uint64_t current_node);
        void _GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity);
        void _GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity);
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####


        //#### HELPER FUNCTIONS FOR DLS ###
        void _GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set);
        void _GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set);
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
        bool DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) ;
        //#### DLS ####
        string CreateFolder();
        //#### CYCLE ENUMERATION ####
        vector<vector<uint64_t>> FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, vector<int> backtrack_lengths);
        vector<vector<uint64_t>> FindCycleUtil(uint64_t startnode);
        //#### CYCLE ENUMERATION ####
        void MarkCycleNodesUpTo100();
        // 1. Call ChunkStartNodes to chunk the start nodes based on their multiplicity for parallel processing
        // 1.1 ChunkStartNodes will call DepthLevelSearch to determine if there is a cycle in a certain depth
        // 2. Call FindCycleUtil to find the cycles in the graph
        // 3. Call PathWriter to write the cycles to a file
        // 4. Return the number of cycles found
        int FindApproximateCRISPRArrays();
        ~CycleFinder();
};  
#endif