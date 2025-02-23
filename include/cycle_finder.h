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
#include <unordered_set>
//#include "progressbar.hpp"
#include <list>
#include <malloc.h>
#include <bitset>
using namespace std;


class CycleFinder {
    private:
        uint8_t maximal_length;
        uint8_t minimal_length;
        static constexpr size_t MAX_EDGE_COUNT = 4;
        SDBG& sdbg;
        string genome_name;
        uint16_t cluster_bounds;
        vector<bool> visited;
        vector<bool> look_up_table;

        //#### DEVELOPER FUNCTIONS ####
        void _WriteStartNodesToFile(const map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _ReadStartNodesFromFile(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked, const std::string& filename);
        void _WriteMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename);
        //### DEVELOPER FUNCTIONS ####
        
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####
        bool _IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree);
        bool _BackgroundCheck(auto original_node, size_t repeat_multiplicity, auto current_node);
        void _GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity);
        void _GetIncomings(auto node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity);
        //#### HELPER FUNCTIONS FOR CYCLE ENUMERATION ####


        //#### HELPER FUNCTIONS FOR DLS ###
        void _GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set);
        void _GetIncomings(auto node, unordered_set<uint64_t>& incomings_set);
        //#### HELPER FUNCTIONS FOR DLS ####

    public:
        CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name);
        //write a getter for results
        unordered_map<uint64_t, vector<vector<uint64_t>>> results;

        // #### RECURSIVE REDUCTION ####
        vector<uint64_t> CollectTips();
        void RecursiveReduction(uint64_t tip);
        void InvalidateMultiplicityOneNodes();  
        void writeMapToFile(const std::unordered_map<uint64_t,  std::vector<std::vector<uint64_t>>>& cycles, const std::string& filename);
        //#### DLS ####
        size_t ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked);
        bool DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth);
        //#### DLS ####

        //#### CYCLE ENUMERATION ####
        vector<vector<uint64_t>> FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, vector<int> backtrack_lengths);
        vector<vector<uint64_t>> FindCycleUtil(uint64_t startnode);
        //#### CYCLE ENUMERATION ####
  
        // 1. Call ChunkStartNodes to chunk the start nodes based on their multiplicity for parallel processing
        // 1.1 ChunkStartNodes will call DepthLevelSearch to determine if there is a cycle in a certain depth
        // 2. Call FindCycleUtil to find the cycles in the graph
        // 3. Call PathWriter to write the cycles to a file
        // 4. Return the number of cycles found
        int FindApproximateCRISPRArrays();
        ~CycleFinder();
};  
#endif